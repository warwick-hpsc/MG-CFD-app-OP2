// This code is to run thermo-mechanical test cases with non-linearities
#include "Structural.h"
#include "TransientThermal.h"
#include <math.h>
#include <boost/program_options.hpp>
#include <dolfinx.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/graph/partitioners.h>
#include <dolfinx/io/XDMFFile.h>
#include "structures.h"
#include "coupler_config.h"
#include "const_op.h"

namespace po = boost::program_options;
using namespace dolfinx;

class NLProblem
{
public:
  NLProblem(
      std::shared_ptr<fem::Form<PetscScalar>> L,
      std::shared_ptr<fem::Form<PetscScalar>> J,
      std::vector<std::shared_ptr<const fem::DirichletBC<PetscScalar>>> bcs)
      : _l(L), _j(J), _bcs(bcs),
        _b(L->function_spaces()[0]->dofmap()->index_map,
           L->function_spaces()[0]->dofmap()->index_map_bs()),
        //_matA(la::petsc::Matrix(fem::petsc::create_matrix(*J, "baij"), false))
        ////hyperelasticity uses this form
        _matA(la::petsc::Matrix(fem::petsc::create_matrix(*J, "aij"), false))
  {
    auto map = L->function_spaces()[0]->dofmap()->index_map;
    const int bs = L->function_spaces()[0]->dofmap()->index_map_bs();
    std::int32_t size_local = bs * map->size_local();

    std::vector<PetscInt> ghosts(map->ghosts().begin(), map->ghosts().end());
    std::int64_t size_global = bs * map->size_global();
    VecCreateGhostBlockWithArray(map->comm(), bs, size_local, size_global,
                                 ghosts.size(), ghosts.data(),
                                 _b.array().data(), &_b_petsc);
  }

  /// Destructor
  virtual ~NLProblem()
  {
    if (_b_petsc)
      VecDestroy(&_b_petsc);
  }

  auto form()
  {
    return [](Vec x)
    {
      VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
      VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    };
  }

  /// Compute F at current point x
  auto F()
  {
    return [&](const Vec x, Vec)
    {
      // Assemble b and update ghosts
      xtl::span<PetscScalar> b(_b.mutable_array());
      std::fill(b.begin(), b.end(), 0.0);
      fem::assemble_vector<PetscScalar>(b, *_l);
      VecGhostUpdateBegin(_b_petsc, ADD_VALUES, SCATTER_REVERSE);
      VecGhostUpdateEnd(_b_petsc, ADD_VALUES, SCATTER_REVERSE);

      // Set bcs
      Vec x_local;
      VecGhostGetLocalForm(x, &x_local);
      PetscInt n = 0;
      VecGetSize(x_local, &n);
      const PetscScalar* array = nullptr;
      VecGetArrayRead(x_local, &array);
      fem::set_bc<PetscScalar>(b, _bcs, xtl::span<const PetscScalar>(array, n),
                               -1.0);
      VecRestoreArrayRead(x, &array);
    };
  }

  /// Compute J = F' at current point x
  auto J()
  {
    return [&](const Vec, Mat A)
    {
      MatZeroEntries(A);
      fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A, ADD_VALUES), *_j,
                           _bcs);
      MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
      fem::set_diagonal(la::petsc::Matrix::set_fn(A, INSERT_VALUES),
                        *_j->function_spaces()[0], _bcs);
      MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    };
  }

  Vec vector() { return _b_petsc; }

  Mat matrix() { return _matA.mat(); }

private:
  std::shared_ptr<fem::Form<PetscScalar>> _l, _j;
  std::vector<std::shared_ptr<const fem::DirichletBC<PetscScalar>>> _bcs;
  la::Vector<PetscScalar> _b;
  Vec _b_petsc = nullptr;
  la::petsc::Matrix _matA;
};

// Function to compute the near nullspace for elasticity - it is made up
// of the six rigid body modes
MatNullSpace build_near_nullspace(const dolfinx::fem::FunctionSpace& V)
{
  // Create vectors for nullspace basis
  auto map = V.dofmap()->index_map;
  int bs = V.dofmap()->index_map_bs();
  std::vector<dolfinx::la::Vector<PetscScalar>> basis(
      6, dolfinx::la::Vector<PetscScalar>(map, bs));

  // x0, x1, x2 translations
  std::int32_t length_block = map->size_local() + map->num_ghosts();
  for (int k = 0; k < 3; ++k)
  {
    xtl::span<PetscScalar> x = basis[k].mutable_array();
    for (std::int32_t i = 0; i < length_block; ++i)
      x[bs * i + k] = 1.0;
  }

  // Rotations
  auto x3 = basis[3].mutable_array();
  auto x4 = basis[4].mutable_array();
  auto x5 = basis[5].mutable_array();

  const xt::xtensor<double, 2> x = V.tabulate_dof_coordinates(false);
  auto& dofs = V.dofmap()->list().array();
  for (std::size_t i = 0; i < dofs.size(); ++i)
  {
    x3[bs * dofs[i] + 0] = -x(dofs[i], 1);
    x3[bs * dofs[i] + 1] = x(dofs[i], 0);

    x4[bs * dofs[i] + 0] = x(dofs[i], 2);
    x4[bs * dofs[i] + 2] = -x(dofs[i], 0);

    x5[bs * dofs[i] + 2] = x(dofs[i], 1);
    x5[bs * dofs[i] + 1] = -x(dofs[i], 2);
  }

  // Orthonormalize basis
  dolfinx::la::orthonormalize(tcb::make_span(basis));
  if (!dolfinx::la::is_orthonormal(
          tcb::span<const decltype(basis)::value_type>(basis)))
  {
    throw std::runtime_error("Space not orthonormal");
  }

  // Build PETSc nullspace object
  std::int32_t length = bs * map->size_local();
  std::vector<xtl::span<const PetscScalar>> basis_local;
  std::transform(basis.cbegin(), basis.cend(), std::back_inserter(basis_local),
                 [length](auto& x)
                 { return xtl::span(x.array().data(), length); });
  MPI_Comm comm = V.mesh()->comm();
  std::vector<Vec> v = dolfinx::la::petsc::create_vectors(comm, basis_local);
  MatNullSpace ns = dolfinx::la::petsc::create_nullspace(comm, v);
  std::for_each(v.begin(), v.end(), [](auto v) { VecDestroy(&v); });
  return ns;
}

int main_dolfinx(int argc, char *argv[], MPI_Fint comm_int, int instance_number, struct unit units[], struct locators relative_positions[])
{
  // Initialization and Filename definitions
  MPI_Comm fenics_comm = MPI_Comm_f2c(comm_int);
  int internal_rank;
  int internal_size;
  int world_rank;
  Vec send;
  VecScatter scat;
  MPI_Comm_rank(fenics_comm, &internal_rank);
  MPI_Comm_size(fenics_comm, &internal_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PETSC_COMM_WORLD = fenics_comm;
  char filename[2];
  char default_name[24] = "FENICS_output_instance_";
  sprintf(filename, "%d", instance_number);
  strcat(default_name, filename);
  FILE *fp = fopen(default_name, "w");

  //set up the command line arguments for FENICSX
  FILE *input = fopen("Fenics_input", "r");
  int max_bufsize = 200;
  int numinputs = 20;
  char buff[max_bufsize];
  char **argv_fenics = (char **) malloc(numinputs*sizeof(char*));
  for(int i = 0; i < numinputs; i++){
    argv_fenics[i] = (char *) malloc(max_bufsize*sizeof(char));
  }
  int argc_fenics = 0; //first variable should be file name
  if(input == NULL){
    fprintf(stderr, "Can't open input file Fenics_input\n");
    exit(1);
  }else{
    while(fgets(buff, max_bufsize, input) != NULL){
      buff[strcspn(buff, "\n")] = 0;
      char temp[100];
      strcpy(temp, buff);
      char * tmp = strtok(temp, " \n");
      if(tmp != NULL){
        while(tmp != NULL){
          if(strcmp(tmp, " ")){
            argc_fenics++;
            strcpy(argv_fenics[argc_fenics], tmp);
          }
          tmp = strtok(NULL, " \n");
        }
      }else{
        argc_fenics++;
        strcpy(argv_fenics[argc_fenics], buff);
      }
    }
  }
  argc_fenics++;
  strcpy(argv_fenics[0], argv[0]);

  // Remove boost args
  std::vector<char*> argv_petsc;
  for (int i = 0; i < argc_fenics; i++)
  {
    std::string arg_str(argv_fenics[i]);
    if (arg_str.find("Thermal_iterations") == std::string::npos and
        arg_str.find("Structural_iterations") == std::string::npos and
        arg_str.find("Thermal_outer_its") == std::string::npos and
        arg_str.find("Structural_outer_its") == std::string::npos and
        arg_str.find("Thermomech_simulation") == std::string::npos and
        arg_str.find("num_vertices") == std::string::npos)
      argv_petsc.push_back(argv_fenics[i]);
  }
  int argc_ = argv_petsc.size();
  char** argv_ = argv_petsc.data();

  // Initialization of logging and petsc
  dolfinx::init_logging(argc_, argv_);
  PetscInitialize(&argc_, &argv_, nullptr, nullptr);

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "print usage message")(
                     "Thermal_iterations", po::value<std::string>(),
                     "The number of inner thermal Krylov iterations")(
                     "Structural_iterations", po::value<std::string>(),
                     "The number of inner structural Krylov iterations")(
                     "Thermal_outer_its", po::value<std::string>(),
                     "The number of outer thermal Newton iterations")(
                     "Structural_outer_its", po::value<std::string>(),
                     "The number of outer structural Newton iterations")(
                     "num_vertices", po::value<std::string>(),
                     "The target number of vertices in the mesh")(
                     "Thermomech_simulation", po::value<std::string>(),
                     "Set to 1 for a Thermomech simulation and 0 for a Thermal");

  po::variables_map vm;
  po::store(po::command_line_parser(argc_fenics, argv_fenics)
                .options(desc)
                .allow_unregistered()
                .run(),
            vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    return 0;
  }

  //loguru::g_stderr_verbosity = loguru::Verbosity_INFO;
  {

    const std::string st_thermal_its = vm["Thermal_iterations"].as<std::string>();
    const std::string st_structural_its = vm["Structural_iterations"].as<std::string>();
    const std::string st_thermal_outer_its = vm["Thermal_outer_its"].as<std::string>();
    const std::string st_structural_outer_its = vm["Structural_outer_its"].as<std::string>();
    const std::string st_num_vertices = vm["num_vertices"].as<std::string>();
    const std::string st_thermomech = vm["Thermomech_simulation"].as<std::string>();

    const int thermal_its = std::stoi(st_thermal_its);
    const int structural_its = std::stoi(st_structural_its);
    const int thermal_outer_its = std::stoi(st_thermal_outer_its);
    const int structural_outer_its = std::stoi(st_structural_outer_its);
    const int num_vertices = std::stoi(st_num_vertices);
    const int strut_flag = std::stoi(st_thermomech);

    if (internal_rank == 0)
      fprintf(fp, "Creating Mesh ...");

      // Set graph partitioner, e.g. ParMETIS, PT-SCOTCH or Kahip
#ifdef HAS_PARMETIS
    auto graph_part = dolfinx::graph::parmetis::partitioner();
#elif HAS_PTSCOTCH
    auto graph_part = dolfinx::graph::scotch::partitioner(
        dolfinx::graph::scotch::strategy::scalability);
#elif HAS_KAHIP
    auto graph_part = dolfinx::graph::kahip::partitioner();
#else
#error "No mesh partitioner has been selected"
#endif

    auto cell_part = dolfinx::mesh::create_cell_partitioner(graph_part);

    float tmp = round(std::cbrt(num_vertices));
    long unsigned int mesh_vertices = int(tmp) - 1;
    double left_corner = 10.0;

    std::shared_ptr<mesh::Mesh> mesh = std::make_shared<mesh::Mesh>(mesh::create_box(
                                        fenics_comm, {{{0.0, 0.0, 0.0},
                                        {left_corner, left_corner, left_corner}}},
                                        {mesh_vertices, mesh_vertices, mesh_vertices},
                                        mesh::CellType::tetrahedron, mesh::GhostMode::none,
                                        cell_part));

    mesh->topology_mutable().create_entities(2);
    mesh->topology_mutable().create_connectivity(2, 3);

    int tdim = mesh->topology().dim();
    auto domain2_indices2 = locate_entities(*mesh, tdim,
                           [](auto&& x) { return xt::isfinite(xt::row(x, 0)); });
    int num_domain_indices = domain2_indices2.size();
    std::vector<int> domain2_indices(num_domain_indices);
    std::iota(std::begin(domain2_indices), std::end(domain2_indices), 0);
    std::vector<int> domain2_values(domain2_indices.size(), 2);
    auto domain2 = mesh::MeshTags<int>(mesh, tdim, domain2_indices, domain2_values);

    int fdim = mesh->topology().dim()-1;
    auto facet2_indices2 = locate_entities(*mesh, fdim,
                           [](auto&& x) { return xt::isfinite(xt::row(x, 0)); });
    int num_facet_indices = facet2_indices2.size();
    std::vector<int> facet2_indices(num_facet_indices);
    std::iota(std::begin(facet2_indices), std::end(facet2_indices), 0);
    std::vector<int> facet2_values(facet2_indices.size(), 512);
    auto facet2 = mesh::MeshTags<int>(mesh, fdim, facet2_indices, facet2_values);

    //============================= Transient Thermal
    //==============================//

    // Create function space
    common::Timer t_11th("03 Transient Thermal: Setup");
    auto V_th = std::make_shared<fem::FunctionSpace>(fem::create_functionspace(
        functionspace_form_TransientThermal_F, "T", mesh));

    // Create solution function
    auto u_th = std::make_shared<fem::Function<PetscScalar>>(V_th);
    u_th->x()->set(293.15);
    // Create initial solution function

    auto u0_th = std::make_shared<fem::Function<PetscScalar>>(V_th);
    u0_th->x()->set(293.15);

    auto t = std::make_shared<fem::Constant<PetscScalar>>(0.0);
    auto dt = std::make_shared<fem::Constant<PetscScalar>>(100);

    // Define Variational Problem
    if (internal_rank == 0)
      fprintf(fp, "Building Thermal Bilinear forms ...\n");

    auto a_th = std::make_shared<fem::Form<PetscScalar>>(
        fem::create_form<PetscScalar>(
            *form_TransientThermal_J, {V_th, V_th},
            {{"T", u_th}, {"T0", u0_th}}, {{"time", t}, {"dt", dt}},
            {{fem::IntegralType::cell, &domain2},
             {fem::IntegralType::exterior_facet, &facet2}}));

    if (internal_rank == 0)
      fprintf(fp, "Building Thermal Linear forms ... \n");
    auto L_th = std::make_shared<fem::Form<PetscScalar>>(
        fem::create_form<PetscScalar>(
            *form_TransientThermal_F, {V_th}, {{"T", u_th}, {"T0", u0_th}},
            {{"time", t}, {"dt", dt}},
            {{fem::IntegralType::cell, &domain2},
             {fem::IntegralType::exterior_facet, &facet2}}));

    if (internal_rank == 0)
      fprintf(fp, "Create matrix ... \n");

    // Initialize stiffness matrix & vector
    la::petsc::Matrix A_th
        = la::petsc::Matrix(fem::petsc::create_matrix(*a_th), false);

    if (internal_rank == 0)
      fprintf(fp, "Create vector ... \n");

    la::petsc::Vector b_th(
        *L_th->function_spaces()[0]->dofmap()->index_map,
        L_th->function_spaces()[0]->dofmap()->index_map_bs());

    if (internal_rank == 0)
      fprintf(fp, "Create nonlinear problem ... \n");

    // Configure Solver
    NLProblem problem_th(L_th, a_th, {});

    nls::petsc::NewtonSolver nlth_solver(mesh->comm());
    nlth_solver.atol = -1.0;
    nlth_solver.rtol = -1.0;
    nlth_solver.error_on_nonconvergence = false;
    nlth_solver.max_it = thermal_outer_its;
    std::string th_its = std::to_string(thermal_its);

    // To change newtonsolver krylov solver settings
    const la::petsc::KrylovSolver& kspsolver = nlth_solver.get_krylov_solver();
    std::string prefix_t = kspsolver.get_options_prefix();
    la::petsc::options::set(prefix_t + "ksp_type", "cg");
    la::petsc::options::set(prefix_t + "ksp_max_it", th_its);
    la::petsc::options::set(prefix_t + "ksp_convergence_test", "skip");
    la::petsc::options::set(prefix_t + "pc_type", "hypre");
    la::petsc::options::set(prefix_t + "pc_hypre_type", "boomeramg");
    la::petsc::options::set(prefix_t + "pc_hypre_boomeramg_strong_threshold", 0.7);
    la::petsc::options::set(prefix_t + "pc_hypre_boomeramg_agg_nl", 4);
    la::petsc::options::set(prefix_t + "pc_hypre_boomeramg_agg_num_paths", 2);
    la::petsc::options::set(prefix_t + "ksp_monitor", "");
    la::petsc::options::set(prefix_t + "ksp_view", "");
    kspsolver.set_from_options();

    dolfinx::io::XDMFFile xdmf_file_th(fenics_comm, "Temperatures.xdmf",
                                       "w");
    xdmf_file_th.write_mesh(*mesh);

    t_11th.stop();

    //============================= Steady Structural
    //==============================//

    common::Timer t_12th("04 Structral: Setup");

    // Create function space and Define variational problem
    auto V_st = std::make_shared<fem::FunctionSpace>(
        dolfinx::fem::create_functionspace(functionspace_form_Structural_F, "u",
                                           mesh));

    // Create solution function
    auto u_st = std::make_shared<fem::Function<PetscScalar>>(V_st);

    // Define variational problem
    if (internal_rank == 0)
      fprintf(fp, "Building Structural Bilinear forms ... \n");

    // auto a_st =
    auto a_st = std::make_shared<fem::Form<PetscScalar>>(
        fem::create_form<PetscScalar>(
            *form_Structural_J, {V_st, V_st}, {{"T", u_th}, {"u", u_st}},
            {{"time", t}},
            {{fem::IntegralType::cell, &domain2},
             {fem::IntegralType::exterior_facet, &facet2}}));

    if (internal_rank == 0)
      fprintf(fp, "Building Structural Linear forms ... \n");

    auto L_st = std::make_shared<fem::Form<PetscScalar>>(
        fem::create_form<PetscScalar>(
            *form_Structural_F, {V_st}, {{"T", u_th}, {"u", u_st}},
            {{"time", t}},
            {{fem::IntegralType::cell, &domain2},
             {fem::IntegralType::exterior_facet, &facet2}}));

    //Set up boundary conditions
    std::vector<std::shared_ptr<const fem::DirichletBC<PetscScalar>>> bcs_st;
    auto VD = std::make_shared<fem::Function<PetscScalar>>(V_st);
    std::fill(VD->x()->mutable_array().begin(),
              VD->x()->mutable_array().end(), 0);
    auto bound_facets = dolfinx::mesh::locate_entities_boundary(
         *mesh, mesh->topology().dim() - 1,
         [](auto&& x) -> xt::xtensor<bool, 1>
         {
           auto x0 = xt::row(x, 0);
           return xt::isclose(x0, 0.0) or xt::isclose(x0, 10.0);
         });
    auto bc_1 = std::make_shared<const fem::DirichletBC<PetscScalar>>(VD, dolfinx::fem::locate_dofs_topological(*V_st, mesh->topology().dim() - 1, bound_facets));
    bcs_st.push_back(bc_1);

    // Initialize stiffness matrix & vector
    la::petsc::Matrix A_st
        = la::petsc::Matrix(fem::petsc::create_matrix(*a_st), false);
    la::petsc::Vector b_st(
        *L_st->function_spaces()[0]->dofmap()->index_map,
        L_st->function_spaces()[0]->dofmap()->index_map_bs());

    // Create near null space basis (required for smoothed aggregation AMG).
    MatNullSpace ns = build_near_nullspace(*V_st);
    // Attach near nullspace to matrix
    MatSetNearNullSpace(A_st.mat(), ns);
    MatNullSpaceDestroy(&ns);

    // Configure Solver
    NLProblem problem_st(L_st, a_st, bcs_st);
    nls::petsc::NewtonSolver nlst_solver(mesh->comm());
    nlst_solver.atol = -1.0;
    nlst_solver.rtol = -1.0;
    nlst_solver.error_on_nonconvergence = false;
    nlst_solver.max_it = structural_outer_its;
    std::string st_its = std::to_string(structural_its);

    // To change newtonsolver krylov solver settings
    const la::petsc::KrylovSolver& kspsolver2 = nlst_solver.get_krylov_solver();
    std::string prefix = kspsolver2.get_options_prefix();
    la::petsc::options::set(prefix + "ksp_type", "cg");
    la::petsc::options::set(prefix + "ksp_view", "");
    la::petsc::options::set(prefix + "ksp_monitor", "");
    la::petsc::options::set(prefix + "ksp_max_it", st_its);
    la::petsc::options::set(prefix + "ksp_convergence_test", "skip");
    la::petsc::options::set(prefix + "pc_type", "gamg");
    la::petsc::options::set(prefix + "pc_gamg_coarse_eq_limit", 1000);
    la::petsc::options::set(prefix + "mg_levels_ksp_type", "chebyshev");
    la::petsc::options::set(prefix + "mg_levels_pc_type", "jacobi");
    la::petsc::options::set(prefix + "pc_gamg_esteig_ksp_type", "cg");
    la::petsc::options::set("-matptap_via", "scalable");
    kspsolver2.set_from_options();

    dolfinx::io::XDMFFile xdmf_file_st(fenics_comm, "Displacements.xdmf",
                                       "w");
    xdmf_file_st.write_mesh(*mesh);

    t_12th.stop();

    //Find own structure
    int fenics_unit_num = relative_positions[world_rank].placelocator;
    int unit_count = 0;
    int work_count = 1; //since units start from 1
    bool found = false;
    while(!found){
      if(units[unit_count].type == 'F' && fenics_unit_num == work_count){
        found=true;
      }else{
        if(units[unit_count].type != 'C'){
          work_count++;
        }
        unit_count++;
      }
    }

    la::petsc::Vector _u_th(la::petsc::create_vector_wrap(*u_th->x()), false);
    la::petsc::Vector _u_st(la::petsc::create_vector_wrap(*u_st->x()), false);
    auto uth_span = u_th->x()->array();
    auto u0th_span = u0_th->x()->mutable_array();
    VecScatterCreateToZero(_u_th.vec(),&scat,&send);

    //get and send boundary sizes
    PetscInt size_u;
    VecGetSize(_u_th.vec(), &size_u);
    double nodes_size = round(size_u * 0.025); //TEMP boundary 2.5% of grid.
    if(internal_rank==0)
      fprintf(fp, "boudary size %f\n", nodes_size);

    double *p_variables_data = (double*) malloc(nodes_size * NVAR * sizeof(double));
    double *p_variables_recv = (double*) malloc(nodes_size * NVAR * sizeof(double));

    int total_coupler_unit_count = units[unit_count].coupler_ranks.size();
    int ranks_per_coupler;
    if (internal_rank == 0){
      for(int z = 0; z < total_coupler_unit_count; z++){
        ranks_per_coupler = units[unit_count].coupler_ranks[z].size();
        for(int z2 = 0; z2 < ranks_per_coupler; z2++){
          //Send the node sizes to each of the coupler ranks of each coupler unit.
          MPI_Send(&nodes_size, 1, MPI_DOUBLE, units[unit_count].coupler_ranks[z][z2], 0, MPI_COMM_WORLD);
        }
      }
    }


    if(internal_rank == 0)
      fprintf(fp, "Entering transient thermo-mechanical iteration loop... \n");

    //Calculate number of cycles
    int fenics_cycles = coupler_cycles*fenics_conversion_factor;
    common::Timer t_13th("05 Compute");

    for(int i = 0; i < fenics_cycles; i++){
      t->value[0] = t->value[0] + dt->value[0];
      if(internal_rank == 0){
        fprintf(fp, "Starting FEniCS X cycle %d of %d\n",i+1, fenics_cycles);
        fprintf(fp, "Solving for temperature at time: %f\n", t->value[0]);
      }

      common::Timer t_15th("07 pure Compute");
      // Solving for thermal
      nlth_solver.setF(problem_th.F(), problem_th.vector());
      nlth_solver.setJ(problem_th.J(), problem_th.matrix());
      nlth_solver.set_form(problem_th.form());

      nlth_solver.solve(_u_th.vec());
      xdmf_file_th.write_function(*u_th, t->value[0]);

      // Copy solution from this time step to previous time step
      std::copy(uth_span.cbegin(), uth_span.cend(), u0th_span.begin());

      if(strut_flag == 1){
        if(internal_rank == 0){
          fprintf(fp, "Solving displacement for time: %f\n", t->value[0]);
        }

        // Solving for structural
        nlst_solver.setF(problem_st.F(), problem_st.vector());
        nlst_solver.setJ(problem_st.J(), problem_st.matrix());
        nlst_solver.set_form(problem_st.form());
        nlst_solver.solve(_u_st.vec());
        xdmf_file_st.write_function(*u_st, t->value[0]);
      }
      t_15th.stop();
      //Send data
      if((i % fenics_conversion_factor == 0) || (hide_search == true && ((i % fenics_conversion_factor) == fenics_conversion_factor - 1))){
        VecScatterBegin(scat, _u_th.vec(), send, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scat, _u_th.vec(), send, INSERT_VALUES, SCATTER_FORWARD);
        if(internal_rank == 0){
          PetscScalar *send_array;
          VecGetArray(send, &send_array);
          for(int k = 0; k < NVAR; k++){
            for(int j = 0; j < nodes_size; j++){
              p_variables_data[(int) (nodes_size*k)+j] = send_array[j];
            }
          }
          VecRestoreArray(send, &send_array);
          printf("FEniCS X cycle %d comms starting\n", i+1);
          for(int j = 0; j < total_coupler_unit_count; j++){
            int coupler_rank = units[unit_count].coupler_ranks[j][0];
            int coupler_position = relative_positions[coupler_rank].placelocator;
            found = false;
            int unit_count_2 = 0;
            int coupler_unit_count = 1;
            while(!found){//this is used to find out the unit index of the coupler unit
              if(units[unit_count_2].type == 'C' && coupler_unit_count == coupler_position){
                found=true;
              }else{
                if(units[unit_count_2].type == 'C'){
                  coupler_unit_count++;
                }
                unit_count_2++;
              }
            }
            int coupler_vars = 0;
            if(units[unit_count_2].coupling_type == 'S'){
              coupler_vars = 5;
            }else if(units[unit_count_2].coupling_type == 'C'){
              coupler_vars = 1;
            }
			common::Timer t_15th("06 waiting");
            if(hide_search == true){
              if((i % (search_freq*fenics_conversion_factor)) == 0){
                MPI_Send(p_variables_data, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
				t_15th.stop();
              }else if((i % fenics_conversion_factor) == fenics_conversion_factor - 1){
				common::Timer t_14th("06 Coupling");
                MPI_Recv(p_variables_recv, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				t_14th.stop();
              }else{
                MPI_Send(p_variables_data, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
				t_15th.stop();
				common::Timer t_14th("06 Coupling");
                MPI_Recv(p_variables_recv, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                t_14th.stop();
			  }
            }else{
              MPI_Send(p_variables_data, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
			  t_15th.stop();
			  common::Timer t_14th("06 Coupling");
              MPI_Recv(p_variables_recv, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              t_14th.stop();
            }
          }
        }
      }
      if(internal_rank == 0){
        printf("FEniCS X cycle %d comms ending\n", i+1);
        fprintf(fp, "Ending FEniCS X cycle %d of %d\n", i+1, fenics_cycles);
      }
    }
    t_13th.stop();

    if(strut_flag == 1){
      xdmf_file_th.close();
    }
    xdmf_file_th.close();
    // Display timings (reduced using MPI::Max)
    Table times = timings({TimingType::wall});
    if(internal_rank == 0){
      fprintf(fp, times.str().c_str());
    }
  }
  MPI_Barrier(fenics_comm);
  PetscFinalize();
  MPI_Finalize();
  exit(0);
}
