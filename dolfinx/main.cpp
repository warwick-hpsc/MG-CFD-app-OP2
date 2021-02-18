#include "poisson.h"
#include "structures.h"
#include "const_op.h"
#include "petscvec.h"
#include <stdio.h>
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/function/Constant.h>

using namespace dolfinx;

int main_dolfinx(int argc, char* argv[], MPI_Fint comm_int, int instance_number, struct unit units[], struct locators relative_positions[])
{
  MPI_Comm fenics_comm = MPI_Comm_f2c(comm_int);
  int internal_rank;
  int internal_size;
  int worldrank;
  Vec send;
  VecScatter scat;
  MPI_Comm_rank(fenics_comm, &internal_rank);
  MPI_Comm_size(fenics_comm, &internal_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);
  PETSC_COMM_WORLD = fenics_comm;
  common::SubSystemsManager::init_logging(argc, argv);
  common::SubSystemsManager::init_petsc(argc, argv);

  {
    // Create mesh and function space
    auto cmap = fem::create_coordinate_map(create_coordinate_map_poisson);
    auto mesh = std::make_shared<mesh::Mesh>(generation::RectangleMesh::create(
        fenics_comm,
        {Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 1.0, 0.0)},
        {100, 100}, cmap, mesh::GhostMode::shared_facet));

    auto V = fem::create_functionspace(create_functionspace_form_poisson_a, "u",
                                       mesh);

    // Prepare and set Constants for the bilinear form
    auto kappa = std::make_shared<function::Constant<PetscScalar>>(2.0);
    auto f = std::make_shared<function::Function<PetscScalar>>(V);
    auto g = std::make_shared<function::Function<PetscScalar>>(V);

    // Define variational forms
    auto a = fem::create_form<PetscScalar>(create_form_poisson_a, {V, V}, {},
                                           {{"kappa", kappa}}, {});
    auto L = fem::create_form<PetscScalar>(create_form_poisson_L, {V},
                                           {{"f", f}, {"g", g}}, {}, {});


    // FIXME: zero function and make sure ghosts are updated
    // Define boundary condition
    auto u0 = std::make_shared<function::Function<PetscScalar>>(V);

    const auto bdofs = fem::locate_dofs_geometrical({*V}, [](auto& x) {
      static const double epsilon = std::numeric_limits<double>::epsilon();
      return (x.row(0).abs() < 10.0 * epsilon
              or (x.row(0) - 1.0).abs() < 10.0 * epsilon);
    });

    std::vector bc{
        std::make_shared<const fem::DirichletBC<PetscScalar>>(u0, bdofs)};

    f->interpolate([](auto& x) {
      auto dx = Eigen::square(x - 0.5);
      return 10.0 * Eigen::exp(-(dx.row(0) + dx.row(1)) / 0.02);
    });
    g->interpolate([](auto& x) { return Eigen::sin(5 * x.row(0)); });
        
    #include "coupler_config.h"
    char filename[2];
    char default_name[24] = "FENICS_output_instance_";
    sprintf(filename,"%d",instance_number);
    strcat(default_name, filename);

    FILE *fp = fopen(default_name,"w");

    //Find own self.
    int fenics_unit_num = relative_positions[worldrank].placelocator;
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
    
    function::Function<PetscScalar> u(V);
    la::PETScMatrix A = fem::create_matrix(*a);
    la::PETScVector b(*L->function_spaces()[0]->dofmap()->index_map,
                      L->function_spaces()[0]->dofmap()->index_map_bs());

    VecScatterCreateToZero(u.vector(),&scat,&send);

    //get and send boundary conditions.
    int size_u;
    double nodes_sizes[4] = {0.0, 0.0, 0.0, 0.0};
    VecGetSize(u.vector(), &size_u);
    nodes_sizes[0] = round(size_u * 0.025); //TEMP boundary 2.5% of grid.
    
    double *p_variables_data;
    double *p_variables_recv_l0, *p_variables_recv_l1, *p_variables_recv_l2, *p_variables_recv_l3;    

    p_variables_data = (double*) malloc(nodes_sizes[0] * NVAR * sizeof(double));

    p_variables_recv_l0 = (double*) malloc(nodes_sizes[0] * NVAR * sizeof(double));
    p_variables_recv_l1 = (double*) malloc(nodes_sizes[1] * NVAR * sizeof(double));
    p_variables_recv_l2 = (double*) malloc(nodes_sizes[2] * NVAR * sizeof(double));
    p_variables_recv_l3 = (double*) malloc(nodes_sizes[3] * NVAR * sizeof(double));
    
    
    int total_coupler_unit_count = units[unit_count].coupler_ranks.size();
    if(internal_rank == 0){
        for(int z = 0; z < total_coupler_unit_count; z++){
            int ranks_per_coupler = units[unit_count].coupler_ranks[z].size();
            for(int z2 = 0; z2 < ranks_per_coupler; z2++){
                MPI_Send(nodes_sizes, 4, MPI_DOUBLE, units[unit_count].coupler_ranks[z][z2], 0, MPI_COMM_WORLD);
            }
        }
    }

    

    for(int i =0 ; i < mgcycles ; i++){
        //compute sol
        if(internal_rank == 0){
            fprintf(fp,"performing fenics cycle %d of %d\n",i+1,mgcycles); 
        }
        MatZeroEntries(A.mat());
        fem::assemble_matrix(la::PETScMatrix::add_fn(A.mat()), *a, bc);
        fem::add_diagonal(la::PETScMatrix::add_fn(A.mat()), *V, bc);
        MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);

        VecSet(b.vec(), 0.0);
        VecGhostUpdateBegin(b.vec(), INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(b.vec(), INSERT_VALUES, SCATTER_FORWARD);
        fem::assemble_vector_petsc(b.vec(), *L);
        fem::apply_lifting_petsc(b.vec(), {a}, {{bc}}, {}, 1.0);
        VecGhostUpdateBegin(b.vec(), ADD_VALUES, SCATTER_REVERSE);
        VecGhostUpdateEnd(b.vec(), ADD_VALUES, SCATTER_REVERSE);
        fem::set_bc_petsc(b.vec(), bc, nullptr);

        la::PETScKrylovSolver lu(fenics_comm);
        la::PETScOptions::set("ksp_type", "preonly");
        la::PETScOptions::set("pc_type", "lu");
        lu.set_from_options();

        lu.set_operator(A.mat());
        lu.solve(u.vector(), b.vec());
        if(i % conversion_factor == 0){
            VecScatterBegin(scat, u.vector(), send, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scat, u.vector(), send, INSERT_VALUES, SCATTER_FORWARD);            
            if(internal_rank == 0){
                int size1;
                PetscScalar *array;
                VecGetSize(send, &size1);
                VecGetArray(send, &array);
                for(int j=0; j<nodes_sizes[0];j++){
                    for(int k=0; k < NVAR; k++){
                        if(k == 0){
                            p_variables_data[NVAR*j+k] = array[j];
                        }else{
                            p_variables_data[NVAR*j+k] = 0;
                        }
                    }
                }
                VecRestoreArray(send, &array);
                printf("starting fenics comms on cycle %d of %d\n",i,mgcycles);
                for(int j = 0; j < total_coupler_unit_count; j++){
                    int coupler_rank = units[unit_count].coupler_ranks[j][0];
                    MPI_Send(p_variables_data, nodes_sizes[0] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
                    MPI_Send(p_variables_data, nodes_sizes[1] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
                    MPI_Send(p_variables_data, nodes_sizes[2] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
                    MPI_Send(p_variables_data, nodes_sizes[3] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
                    MPI_Recv(p_variables_recv_l0, nodes_sizes[0] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(p_variables_recv_l1, nodes_sizes[1] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(p_variables_recv_l2, nodes_sizes[2] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(p_variables_recv_l3, nodes_sizes[3] * NVAR, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                printf("ending fenics comms on cycle %d of %d\n",i,mgcycles);
            }
            if(i+1 == mgcycles){
                io::VTKFile file("u.pvd");
                file.write(u);
            }    
        }
    }
  }

  common::SubSystemsManager::finalize_petsc();
  return 0;
}
