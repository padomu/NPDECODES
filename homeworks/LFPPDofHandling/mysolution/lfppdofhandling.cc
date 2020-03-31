/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "lfppdofhandling.h"

#include <array>
#include <memory>

#include <Eigen/Dense>

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/utils/utils.h"

namespace LFPPDofHandling {

/* SAM_LISTING_BEGIN_1 */
std::array<std::size_t, 3>
countEntityDofs(const lf::assemble::DofHandler &dofhandler) {
  std::array<std::size_t, 3> entityDofs;
  //====================

  // Get the mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh(); 

  // Loop through codimensions
  for( std::size_t codim = 0; codim <=2; ++codim) {
    entityDofs[codim] = 0;
    // Loop through all entities of given codim
    for( const auto *el : mesh->Entities(codim) ) {
      if( el->RefEl() == lf::base::RefEl::kQuad() ) {
        throw "only triangular meshes are allowed!";
      }
      entityDofs[codim] += dofhandler.NumInteriorDofs(*el);
    }
  }

  //====================
  return entityDofs;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::size_t countBoundaryDofs(const lf::assemble::DofHandler &dofhandler) {
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  // given an entity, bd\_flags(entity) == true, if the entity is on the
  // boundary
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh));
  std::size_t no_dofs_on_bd = 0;
  //====================

  // Edges on boundary
  for( const auto *edge : mesh->Entities(1) ) {
    if( bd_flags(*edge) )
      no_dofs_on_bd += dofhandler.NumInteriorDofs(*edge);
  }

  // Nodes on boundary
  for( const auto *node : mesh->Entities(2) ) {
    if( bd_flags(*node) )
      no_dofs_on_bd += dofhandler.NumInteriorDofs(*node);
  }

  //====================
  return no_dofs_on_bd;
}
/* SAM_LISTING_END_2 */

// clang-format off
/* SAM_LISTING_BEGIN_3 */
double integrateLinearFEFunction(
    const lf::assemble::DofHandler& dofhandler,
    const Eigen::VectorXd& mu) {
  double I = 0.0;
  //====================

  // Loop over all cells
  // - Get global node indices
  // - I += u1+u2+u3
  // - I *= K/3.0

  // Get mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();

  // Loop over all cells
  for( const auto *cell : mesh->Entities(0) ) {
    // Check if the FE space is really S_1^0
    if( dofhandler.NumLocalDofs(*cell) != 3 ) {
      throw "Not a S_1^0 space!";
    }

    // Get all dofs for nodes of cell
    // Note: We are in S_0^1 so the "interior" dofs
    // happen to match the nodes.
    auto int_dofs = dofhandler.GlobalDofIndices(*cell);

    // Loop over all those (three) interior dofs
    for( auto dof_idx_p = int_dofs.begin();
         dof_idx_p < int_dofs.end();
         ++dof_idx_p ) {

      const double K = lf::geometry::Volume(*(cell->Geometry())) / 3.0;

      I += K*mu(*dof_idx_p);
    }
    // 
  }


  //====================
  return I;
}
/* SAM_LISTING_END_3 */
// clang-format on

/* SAM_LISTING_BEGIN_4 */
double integrateQuadraticFEFunction(const lf::assemble::DofHandler &dofhandler,
                                    const Eigen::VectorXd &mu) {
  double I = 0;
  //====================

  // Get the mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh(); 

  // Loop over all cells
  for( const auto *cell : mesh->Entities(0) ) {
    // Check if the FE space really is S_2^0
    if( dofhandler.NumLocalDofs(*cell) != 6 ) {
      throw "Not a S_2^0 FE space!";
    }

    // Get area of cell
    const double K = lf::geometry::Volume(*(cell->Geometry())) / 3.0;

    // Get internal dofs
    auto int_dofs = dofhandler.GlobalDofIndices(*cell);

    // See solution and calculate it with 2.7.5.5!
    // => int_k bk^j dx = 0 for j=1,2,3 => SKIP

    I += K * ( mu(int_dofs[3]) + mu(int_dofs[4]) + mu(int_dofs[5]) );
  }

  //====================
  return I;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd
convertDOFsLinearQuadratic(const lf::assemble::DofHandler &dofh_Linear_FE,
                           const lf::assemble::DofHandler &dofh_Quadratic_FE,
                           const Eigen::VectorXd &mu) {
  if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
    throw "Underlying meshes must be the same for both DOF handlers!";
  }
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      dofh_Linear_FE.Mesh();                         // get the mesh
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NumDofs()); // initialise empty zeta
  // safety guard: always set zero if you're not sure to set every entry later
  // on for us this shouldn't be a problem, but just to be sure
  zeta.setZero();

  for (const auto *cell : mesh->Entities(0)) {
    //====================

    //====================
  }
  return zeta;
}
/* SAM_LISTING_END_5 */

} // namespace LFPPDofHandling
