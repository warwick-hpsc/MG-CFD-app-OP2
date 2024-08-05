#include <cmath>
#include <dolfinx.h>
#include <dolfinx/common/sort.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/refinement/plaza.h>
#include <dolfinx/refinement/utils.h>
//#include <xtensor/xadapt.hpp>
//#include <xtensor/xarray.hpp>
//#include <xtensor/xio.hpp>
//#include <xtensor/xview.hpp>

namespace
{
/// Read a mesh
/// @param[in] filename The file name
/// @return The tuple (mesh, domain tags, facet tags)
auto read_mesh(const std::string& filename)
{
  // Read and create mesh
  dolfinx::io::XDMFFile file(MPI_COMM_WORLD, filename, "r");
  dolfinx::fem::CoordinateElement cmap = dolfinx::fem::CoordinateElement(
      dolfinx::mesh::CellType::tetrahedron, 1);
  auto geometry = file.read_geometry_data("geometry");
  auto topology = file.read_topology_data("volume markers");

  dolfinx::graph::AdjacencyList<std::int64_t> cells_adj
      = dolfinx::graph::regular_adjacency_list(topology.first,
                                               topology.second.back());
  auto mesh = std::make_shared<dolfinx::mesh::Mesh>(dolfinx::mesh::create_mesh(
      MPI_COMM_WORLD, cells_adj, cmap, geometry.first, geometry.second,
      dolfinx::mesh::GhostMode::none));
      
  mesh->topology_mutable().create_entities(2);
  mesh->topology_mutable().create_connectivity(2, 3);

  // Create entity-vertex connectivity
  constexpr int tdim = 3;
  mesh->topology_mutable().create_entities(tdim - 1);
  mesh->topology_mutable().create_connectivity(tdim - 1, tdim);
  // Read domain meshtags
  if (dolfinx::MPI::rank(mesh->comm()) == 0)
    std::cout << "Reading domain MeshTags ..." << std::endl;
  auto domain1 = file.read_meshtags(mesh, "volume markers");

  // Read facet meshtags
  if (dolfinx::MPI::rank(mesh->comm()) == 0)
    std::cout << "Reading facet MeshTags ..." << std::endl;
  auto facet1 = file.read_meshtags(mesh, "facet markers");

  return std::make_tuple(mesh, domain1, facet1);
}
} // namespace


using namespace dolfinx;
using T = PetscScalar;

int main(int argc, char* argv[])
{
  dolfinx::init_logging(argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);

  {
    // Get 'my' MPI rank
    const int rank = dolfinx::MPI::rank(MPI_COMM_WORLD);

    
    // --- Read mesh
    if (rank == 0)
      std::cout << "Reading Mesh ..." << std::endl;
    std::string mesh_file;
    int target;
    if(argc > 1){
      mesh_file = argv[1];
    }else{
      printf("Please specify a mesh and target number of cells to refine\n");
      exit(1);
    }
    const auto [mesh, domain1, facet1] = read_mesh(mesh_file);

    io::XDMFFile xdmf(mesh->comm(), "cell_marker.xdmf", "w");
    mesh->name="geometry";
    xdmf.write_mesh(*mesh);
    xdmf.write_meshtags(domain1, "/Xdmf/Domain/Grid[@Name='geometry']/Geometry");
    xdmf.write_meshtags(facet1, "/Xdmf/Domain/Grid[@Name='geometry']/Geometry");
    xdmf.close();

    // Refine mesh and transfer a mesh-tag
    mesh->topology_mutable().create_entities(1);
    auto [fine_mesh, parent_cell, parent_facet]
        = dolfinx::refinement::plaza::refine(
            *mesh, false,
            dolfinx::refinement::plaza::RefinementOptions::
                parent_cell_and_facet);

    mesh::MeshTags<std::int32_t> new_meshtag
        = dolfinx::refinement::transfer_cell_meshtag(domain1,
            std::make_shared<const dolfinx::mesh::Mesh>(fine_mesh),
                                                         parent_cell);


    mesh->topology_mutable().create_connectivity(mesh->topology().dim(),
                                                 mesh->topology().dim() - 1);
    fine_mesh.topology_mutable().create_connectivity(
        fine_mesh.topology().dim(), fine_mesh.topology().dim() - 1);
    new_meshtag.name = "volume markers";


    mesh::MeshTags<std::int32_t> new_meshtag_facet
         = dolfinx::refinement::transfer_facet_meshtag(facet1, 
                                                      std::make_shared<const dolfinx::mesh::Mesh>(fine_mesh),
                                                      parent_cell, parent_facet);
    new_meshtag_facet.mesh()->topology_mutable().create_connectivity(
        fine_mesh.topology().dim() - 1, fine_mesh.topology().dim());													  
    new_meshtag_facet.name = "facet markers";


    io::XDMFFile xdmf_fine(mesh->comm(), "refined_cell_marker.xdmf", "w");
    fine_mesh.name="geometry";
    xdmf_fine.write_mesh(fine_mesh);
    xdmf_fine.write_meshtags(new_meshtag, "/Xdmf/Domain/Grid/Geometry");
    xdmf_fine.write_meshtags(new_meshtag_facet, "/Xdmf/Domain/Grid/Geometry");
    xdmf_fine.close();


    
  }

  PetscFinalize();

  return 0;
}


