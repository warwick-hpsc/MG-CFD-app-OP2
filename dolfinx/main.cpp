// This code is to run thermo-mechanical test cases with non-linearities
#include <boost/program_options.hpp>
#include <dolfinx.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/graph/partitioners.h>
#include "Structural.h"
#include "TransientThermal.h"
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
        //_matA(la::PETScMatrix(fem::create_matrix(*J, "baij"), false))  //hyperelasticity uses this form
        _matA(la::PETScMatrix(fem::create_matrix(*J, "aij"), false))
  {
    auto map = L->function_spaces()[0]->dofmap()->index_map;
    const int bs = L->function_spaces()[0]->dofmap()->index_map_bs();
    std::int32_t size_local = bs * map->size_local();

    std::vector<PetscInt> ghosts(map->ghosts().begin(), map->ghosts().end());
    std::int64_t size_global = bs * map->size_global();
    VecCreateGhostBlockWithArray(
        map->comm(common::IndexMap::Direction::forward), bs, size_local,
        size_global, ghosts.size(), ghosts.data(), _b.array().data(),
        &_b_petsc);
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
      fem::assemble_matrix(la::PETScMatrix::set_block_fn(A, ADD_VALUES), *_j,
                           _bcs);
      MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
      fem::set_diagonal(la::PETScMatrix::set_fn(A, INSERT_VALUES),
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
  la::PETScMatrix _matA;
};

dolfinx::la::VectorSpaceBasis build_near_nullspace(const fem::FunctionSpace& V)
{
  // Create vectors for nullspace basis
  auto map = V.dofmap()->index_map;
  int bs = V.dofmap()->index_map_bs();
  const std::int32_t length_block = map->size_local() + map->num_ghosts();
  const std::int32_t length = bs * length_block;
  xt::xtensor<PetscScalar, 2> basis = xt::zeros<PetscScalar>({length, 6});

  // x0, x1, x2 translations
  for (int k = 0; k < 3; ++k)
  {
    for (std::int32_t i = 0; i < length_block; ++i)
      basis(bs * i + k, k) = 1.0;
  }

  // Rotations
  const xt::xtensor<double, 2> x = V.tabulate_dof_coordinates(false);
  auto& dofs = V.dofmap()->list().array();
  for (size_t i = 0; i < dofs.size(); ++i)
  {
    basis(bs * dofs[i] + 0, 3) = -x(dofs[i], 1);
    basis(bs * dofs[i] + 1, 3) = x(dofs[i], 0);

    basis(bs * dofs[i] + 0, 4) = x(dofs[i], 2);
    basis(bs * dofs[i] + 2, 4) = -x(dofs[i], 0);

    basis(bs * dofs[i] + 2, 5) = x(dofs[i], 1);
    basis(bs * dofs[i] + 1, 5) = -x(dofs[i], 2);
  }

  const std::int32_t size = map->size_local() * bs;
  const std::int64_t size_global = map->size_global() * bs;
  std::vector<std::shared_ptr<la::PETScVector>> basis_vec;
  for (int i = 0; i < 6; ++i)
  {
    Vec vec0, vec1;
    xt::xarray<PetscScalar> basis_row = xt::col(basis, i);
    VecCreateMPIWithArray(V.mesh()->mpi_comm(), 3, size, size_global,
                          basis_row.data(), &vec0);
    VecDuplicate(vec0, &vec1);
    VecCopy(vec0, vec1);
    VecDestroy(&vec0);
    basis_vec.push_back(std::make_shared<la::PETScVector>(vec1, false));
  }

  // Create vector space and orthonormalize
  la::VectorSpaceBasis vector_space(basis_vec);
  vector_space.orthonormalize();

  return vector_space;
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

  int strut_flag = 0; //0 for thermal sim, 1 for thermomech sim;

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

  //common::subsystem::init_mpi();
  common::subsystem::init_petsc(argc_fenics, argv_fenics);
  common::subsystem::init_logging(argc_fenics, argv_fenics);

  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "print usage message")(
						      "mesh_file", po::value<std::string>(), "mesh filename");

  po::variables_map vm;
  po::store(po::command_line_parser(argc_fenics, argv_fenics).options(desc).allow_unregistered().run(),vm);
  po::notify(vm);

  if (vm.count("help"))
    {
      std::cout << desc << "\n";
      return 0;
    }

  const std::string mesh_file_name = vm["mesh_file"].as<std::string>();
  const std::string mesh_file = mesh_file_name;
 
  // Read mesh 
  std::shared_ptr<mesh::Mesh> mesh;
    
      if (internal_rank == 0)
        fprintf(fp, "Reading Mesh data ...");

      dolfinx::io::XDMFFile file(fenics_comm, mesh_file, "r");
      fem::CoordinateElement cmap
          = fem::CoordinateElement(mesh::CellType::tetrahedron, 1);

      // xt::xtensor<double, 2> geometry = file.read_geometry_data("geometry");
      // xt::xtensor<std::int64_t, 2> topology
      //     = file.read_topology_data("tetra marker");

      // For 'standard' XDMF files
      xt::xtensor<double, 2> x = file.read_geometry_data("mesh");
      xt::xtensor<std::int64_t, 2> topology = file.read_topology_data("mesh");

      auto [data, offset] = graph::create_adjacency_data(topology);
      graph::AdjacencyList<std::int64_t> cells(std::move(data),
                                               std::move(offset));

      if (internal_rank == 0)
        fprintf(fp, "Creating Mesh ...");

      // Set graph partitioner, e.g. ParMETIS or PT-SCOTCH

      // SCOTCH
      auto graph_part = graph::scotch::partitioner(graph::scotch::strategy::scalability);
      // ParMETIS
      //auto graph_part = graph::parmetis::partitioner(1.01);

      auto cell_part
          = [graph_part](MPI_Comm comm, int nparts, int tdim,
                         const graph::AdjacencyList<std::int64_t>& cells,
                         mesh::GhostMode ghost_mode)
      {
        return mesh::partition_cells_graph(comm, nparts, tdim, cells,
                                           ghost_mode, graph_part);
      };

      // Create distributed mesh
      mesh = std::make_shared<mesh::Mesh>(mesh::create_mesh(
          fenics_comm, cells, cmap, x, mesh::GhostMode::none, cell_part));
       
  mesh->topology_mutable().create_entities(2);
  mesh->topology_mutable().create_connectivity(2, 3);
     
  // Read domain meshtags
  if(internal_rank  == 0)
    fprintf(fp, "Reading domain MeshTags ...\n");
  auto domain1 = file.read_meshtags(mesh, "mesh");

  // Read facet meshtags
  if(internal_rank == 0)
    fprintf(fp, "Reading facet MeshTags ...\n");
  auto facet1 = file.read_meshtags(mesh, "facets");

  //============================= Transient Thermal ==============================//

  // Create function space
  common::Timer t_11th("03 Transient Thermal: Setup");
  auto V_th = fem::create_functionspace(functionspace_form_TransientThermal_F, "T", mesh);

  // Create solution function
  auto u_th  = std::make_shared<fem::Function<PetscScalar>>(V_th);
  VecSet(u_th->vector(), 293.15);
  VecGhostUpdateBegin(u_th->vector(), INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(u_th->vector(), INSERT_VALUES, SCATTER_FORWARD);
  VecScatterCreateToZero(u_th->vector(),&scat,&send);

  // Create initial solution function
  auto u0_th = std::make_shared<fem::Function<PetscScalar>>(V_th);
  VecSet(u0_th->vector(), 293.15);
  VecGhostUpdateBegin(u0_th->vector(), INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(u0_th->vector(), INSERT_VALUES, SCATTER_FORWARD);
 

  auto t = std::make_shared<fem::Constant<PetscScalar>>(0.0);
  auto dt = std::make_shared<fem::Constant<PetscScalar>>(100);
  //const int t_idx_max_th = 1;   //final time point
            
  // Define Variational Problem
  if(internal_rank == 0)
    fprintf(fp, "Building Thermal Bilinear forms ...\n");

  auto a_th =  std::make_shared<fem::Form<PetscScalar>>(fem::create_form<PetscScalar>(*form_TransientThermal_J, {V_th, V_th},
    {{"T", u_th},{"T0", u0_th}},
    {{"time", t},{"dt", dt}},
    {{fem::IntegralType::cell, &domain1},{fem::IntegralType::exterior_facet, &facet1}}));

  if(internal_rank == 0)
    fprintf(fp, "Building Thermal Linear forms ... \n");

  auto L_th = std::make_shared<fem::Form<PetscScalar>>(fem::create_form<PetscScalar>(*form_TransientThermal_F, {V_th},
    {{"T", u_th},{"T0", u0_th}},
    {{"time", t},{"dt", dt}},
    {{fem::IntegralType::cell, &domain1},{fem::IntegralType::exterior_facet, &facet1}}));
                    
  // Initialize stiffness matrix & vector            
  la::PETScMatrix A_th = la::PETScMatrix(fem::create_matrix(*a_th), false);
  la::PETScVector b_th(*L_th->function_spaces()[0]->dofmap()->index_map,
		       L_th->function_spaces()[0]->dofmap()->index_map_bs()); 
            
  // Configure Solver
  NLProblem problem_th(L_th, a_th, {});
  nls::NewtonSolver nlth_solver(mesh->mpi_comm());
  nlth_solver.atol = 1e-4;
 
  // To change newtonsolver krylov solver settings
  const la::PETScKrylovSolver& kspsolver = nlth_solver.get_krylov_solver();
  std::string prefix_t = kspsolver.get_options_prefix();
  la::PETScOptions::set(prefix_t + "ksp_type", "cg");
  la::PETScOptions::set(prefix_t + "ksp_rtol", 1.0e-8);
  //new options
  la::PETScOptions::set(prefix_t + "pc_type", "gamg");
  la::PETScOptions::set(prefix_t + "pc_gamg_coarse_eq_limit", 1000);
  la::PETScOptions::set(prefix_t + "mg_levels_ksp_type", "chebyshev");
  la::PETScOptions::set(prefix_t + "mg_levels_pc_type", "jacobi");
  la::PETScOptions::set(prefix_t + "mg_levels_esteig_ksp_type", "cg");
  /*la::PETScOptions::set(prefix_t + "pc_type", "hypre");
  la::PETScOptions::set(prefix_t + "pc_hypre_type", "boomeramg");
  la::PETScOptions::set(prefix_t + "hypre_boomeramg_strong_threshold", 0.7);
  la::PETScOptions::set(prefix_t + "hypre_boomeramg_agg_nl", 4);
  la::PETScOptions::set(prefix_t + "hypre_boomeramg_agg_num_paths", 2);*/
  kspsolver.set_from_options();
  
  
  dolfinx::io::XDMFFile xdmf_file_th(fenics_comm, "Temperatures.xdmf", "w");                    
  xdmf_file_th.write_mesh(*mesh);
  t_11th.stop();
  //============================= Steady Structural ==============================//

  common::Timer t_12th("04 Strutral: Setup");
  // Create function space and Define variational problem
  auto V_st = dolfinx::fem::create_functionspace(functionspace_form_Structural_F, "u", mesh);

  // Create solution function
  auto u_st = std::make_shared<fem::Function<PetscScalar>>(V_st);
  VecGhostUpdateBegin(u_st->vector(), INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(u_st->vector(), INSERT_VALUES, SCATTER_FORWARD);

  // Define variational problem
  if(internal_rank == 0)
    fprintf(fp, "Building Structural Bilinear forms ...\n");

  //auto a_st = std::make_shared<fem::create_form<PetscScalar>(create_form_Structural_J,{V_st,V_st},
  auto a_st = std::make_shared<fem::Form<PetscScalar>>(fem::create_form<PetscScalar>(*form_Structural_J,{V_st,V_st},
    {{"T", u_th},{"u", u_st}},
    {{"time", t}},
    {{fem::IntegralType::cell, &domain1},{fem::IntegralType::exterior_facet, &facet1}}));

  if(internal_rank == 0)
    fprintf(fp, "Building Structural Linear forms ...\n");

  auto L_st = std::make_shared<fem::Form<PetscScalar>>(fem::create_form<PetscScalar>(*form_Structural_F,{V_st},
    {{"T", u_th}, {"u", u_st}},
    {{"time", t}},
    {{fem::IntegralType::cell, &domain1},{fem::IntegralType::exterior_facet, &facet1}}));

  // Set up boundary condition values
  auto zero_0 = std::make_shared<fem::Function<PetscScalar>>(V_st->sub({0})->collapse().first);
  std::fill(zero_0->x()->mutable_array().begin(),zero_0->x()->mutable_array().end(),0);
  auto zero_1 = std::make_shared<fem::Function<PetscScalar>>(V_st->sub({1})->collapse().first);
  std::fill(zero_1->x()->mutable_array().begin(),zero_1->x()->mutable_array().end(),0);
  auto zero_2 = std::make_shared<fem::Function<PetscScalar>>(V_st->sub({2})->collapse().first);
  std::fill(zero_2->x()->mutable_array().begin(),zero_2->x()->mutable_array().end(),0);

  // Set up boundary condition tags
  const auto facet_array = facet1.values();
  const auto facet_indices = facet1.indices();
  std::vector<std::int32_t> tag0{459, 460};
  std::vector<std::int32_t> tag1{803, 804, 805, 813, 814, 815, 1153, 1154, 1157, 1158};
  std::vector<std::shared_ptr<const fem::DirichletBC<PetscScalar>>> bcs_st;

  for (auto i : tag0)
    {
      auto facets_dim0 = facet1.find(i);
      auto dofs0 = dolfinx::fem::locate_dofs_topological({*V_st->sub({0}),*zero_0->function_space()},mesh->topology().dim()-1,facets_dim0);
      auto x_bc = std::make_shared<const fem::DirichletBC<PetscScalar>>(zero_0, std::move(dofs0), V_st);
      bcs_st.push_back(x_bc);
    }
  for (auto i : tag1)
    {
      auto facets_dim1 = facet1.find(i);
      auto dofs1 = dolfinx::fem::locate_dofs_topological({*V_st->sub({1}),*zero_1->function_space()},mesh->topology().dim()-1,facets_dim1);
      auto y_bc = std::make_shared<const fem::DirichletBC<PetscScalar>>(zero_1, std::move(dofs1), V_st);
      bcs_st.push_back(y_bc);
      auto dofs2 = dolfinx::fem::locate_dofs_topological({*V_st->sub({2}),*zero_2->function_space()},mesh->topology().dim()-1,facets_dim1);
      auto z_bc = std::make_shared<const fem::DirichletBC<PetscScalar>>(zero_2, std::move(dofs2), V_st);
      bcs_st.push_back(z_bc);
    }

  // Initialize stiffness matrix & vector            
  la::PETScMatrix A_st = la::PETScMatrix(fem::create_matrix(*a_st), false);
  la::PETScVector b_st(*L_st->function_spaces()[0]->dofmap()->index_map,
	  	          L_st->function_spaces()[0]->dofmap()->index_map_bs()); 

  // Create near null space basis (required for smoothed aggregation AMG).
  la::VectorSpaceBasis near_null_space = build_near_nullspace(*V_st);

  // Attach near nullspace to matrix (TODO: how to use near nullspace with newton solver?
  A_st.set_near_nullspace(near_null_space);
            
  // Configure Solver
  NLProblem problem_st(L_st, a_st, bcs_st);
  nls::NewtonSolver nlst_solver(mesh->mpi_comm());
  nlst_solver.atol = 1e-4;

  // To change newtonsolver krylov solver settings
  const la::PETScKrylovSolver& kspsolver2 = nlst_solver.get_krylov_solver();
  std::string prefix = kspsolver2.get_options_prefix();
  la::PETScOptions::set(prefix + "ksp_type", "cg");
  la::PETScOptions::set(prefix + "pc_type", "gamg");
  la::PETScOptions::set(prefix + "pc_gamg_coarse_eq_limit", 1000);
  la::PETScOptions::set(prefix + "ksp_rtol", 1e-8);
  la::PETScOptions::set(prefix + "mg_levels_ksp_type", "chebyshev");
  la::PETScOptions::set(prefix + "mg_levels_pc_type", "jacobi");
  la::PETScOptions::set(prefix + "mg_levels_esteig_ksp_type", "cg");
  la::PETScOptions::set(prefix + "matptap_via", "scalable");
  kspsolver2.set_from_options();

  dolfinx::io::XDMFFile xdmf_file_st(fenics_comm, "Displacements.xdmf", "w");
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

  //get and send boundary sizes
  PetscInt size_u;
  VecGetSize(u_th->vector(), &size_u);
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
	fprintf(fp, "Entering transient thermo-mechanical iteration loop...\n");

  //Calculate number of cycles
  int fenics_cycles = coupler_cycles*fenics_conversion_factor; 
  common::Timer t_13th("05 Compute");

  for(int i = 0; i < fenics_cycles; i++){ 
    t->value[0] = t->value[0] + dt->value[0];
    if(internal_rank == 0){
      fprintf(fp, "Starting FEniCS X cycle %d of %d\n",i+1,fenics_cycles);
	  fprintf(fp, "Solving for temperature at time: %f\n", t->value[0]);
    }
    common::Timer t_15th("07 pure Compute");
    // Solving for thermal
    nlth_solver.setF(problem_th.F(), problem_th.vector());
    nlth_solver.setJ(problem_th.J(), problem_th.matrix());
    nlth_solver.set_form(problem_th.form());
    nlth_solver.solve(u_th->vector());
    xdmf_file_th.write_function(*u_th, t->value[0]);
    VecCopy(u_th->vector(), u0_th->vector());

	if(internal_rank == 0){
	  fprintf(fp, "Solving displacement for time: %f\n", t->value[0]);
	}
    if(strut_flag == 1){
      // Solving for structural
	  nlst_solver.setF(problem_st.F(), problem_st.vector());
	  nlst_solver.setJ(problem_st.J(), problem_st.matrix());
	  nlst_solver.set_form(problem_st.form());
	  nlst_solver.solve(u_st->vector());
      xdmf_file_st.write_function(*u_st, t->value[0]);
    }
    t_15th.stop();

    //Send data
    if((i % fenics_conversion_factor == 0) || (hide_search == true && ((i % fenics_conversion_factor) == fenics_conversion_factor - 1))){
      VecScatterBegin(scat, u_th->vector(), send, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(scat, u_th->vector(), send, INSERT_VALUES, SCATTER_FORWARD);
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
          common::Timer t_14th("06 Coupling");
		  if(hide_search == true){
		    if((i % (search_freq*fenics_conversion_factor)) == 0){
			  MPI_Send(p_variables_data, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
		    }else if((i % fenics_conversion_factor) == fenics_conversion_factor - 1){
			  MPI_Recv(p_variables_recv, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }else{
			  MPI_Send(p_variables_data, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
			  MPI_Recv(p_variables_recv, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		  }else{
		    MPI_Send(p_variables_data, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
		    MPI_Recv(p_variables_recv, nodes_size * coupler_vars, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  }
          t_14th.stop();
		}
	  }
    }
    if(internal_rank == 0){
      printf("FEniCS X cycle %d comms ending\n", i+1);
      fprintf(fp, "Ending FEniCS X cycle %d of %d\n", i+1, fenics_cycles);
    }
  }
  t_13th.stop();

  // Display timings
  Table times = timings({TimingType::wall});
  if(internal_rank == 0){
    fprintf(fp, times.str().c_str());
  }
  MPI_Barrier(fenics_comm);
  MPI_Finalize();
  exit(0);
}
