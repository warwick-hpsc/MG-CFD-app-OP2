#include "poisson.h"
#include "structures.h"
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/function/Constant.h>

using namespace dolfinx;

int main_dolfinx(int argc, char* argv[], MPI_Fint comm_int, int instance_number, struct unit units[], struct locators relative_positions[])
{
  MPI_Comm fenics_comm = MPI_Comm_f2c(comm_int);
  PETSC_COMM_WORLD = fenics_comm;
  common::SubSystemsManager::init_logging(argc, argv);
  common::SubSystemsManager::init_petsc(argc, argv);

  {
    // Create mesh and function space
    auto cmap = fem::create_coordinate_map(create_coordinate_map_poisson);
    auto mesh = std::make_shared<mesh::Mesh>(generation::RectangleMesh::create(
        fenics_comm,
        {Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 1.0, 0.0)},
        {32, 32}, cmap, mesh::GhostMode::shared_facet));

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

    // Compute solution
    function::Function<PetscScalar> u(V);
    la::PETScMatrix A = fem::create_matrix(*a);
    la::PETScVector b(*L->function_spaces()[0]->dofmap()->index_map,
                      L->function_spaces()[0]->dofmap()->index_map_bs());

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

    // Save solution in VTK format
    io::VTKFile file("u.pvd");
    file.write(u);
  }

  common::SubSystemsManager::finalize_petsc();
  return 0;
}
