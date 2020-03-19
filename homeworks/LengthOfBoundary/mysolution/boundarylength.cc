/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>


namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double volume = 0.0;
  //====================

  lf::mesh::Mesh *mesh = mesh_p.get();

  // Loop through cells
  for ( const lf::mesh::Entity *cell : mesh->Entities(0) ) {
    volume += lf::geometry::Volume( *cell->Geometry() );
  }

  //====================

  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double length = 0.0;
  //====================

  lf::mesh::Mesh *mesh = mesh_p.get();

  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags
          = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1);

  // Loop through all edges
  for ( const lf::mesh::Entity *edge : mesh_p->Entities(1) ) {
    if ( bd_flags(*edge) )
      length += lf::geometry::Volume( *edge->Geometry() );
  }

  //====================

  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string filename) {
  double volume, length;

  //====================
  // Your code goes here
  //====================

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
