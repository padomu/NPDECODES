/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "lfppdofhandling.h"

#include <Eigen/Dense>
#include <array>
#include <memory>

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/utils/utils.h"

namespace LFPPDofHandling {

/* SAM_LISTING_BEGIN_1 */
std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler &dofhandler) {
  std::array<std::size_t, 3> entityDofs;
  //====================
  
  // Idea: iterate over entities in the mesh and get interior number of dofs for
  // each
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (std::size_t codim = 0; codim <= 2; ++codim) {
    entityDofs[codim] = 0;
    for (const auto *el : mesh->Entities(codim)) {
      if (el->RefEl() == lf::base::RefEl::kQuad()) {
        throw "Only triangular meshes are allowed!";
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
  
  const lf::mesh::Mesh::size_type N_dofs = dofhandler.NumDofs();

  for (lf::assemble::gdof_idx_t i=0; i < N_dofs; ++i) {
    // Get i-th entity
    const lf::mesh::Entity &el = dofhandler.Entity(i);
    if (el.RefEl() == lf::base::RefEl::kQuad()) {
      throw "Only triangular meshes are allowed!";
    }
    if(bd_flags(el)) {
      no_dofs_on_bd++;
    }
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

  const lf::mesh::Mesh::size_type N_dofs = dofhandler.NumDofs();
  const lf::mesh::Mesh& mesh = *dofhandler.Mesh();

  // Loop over triangles
  for (const lf::mesh::Entity* triangle : mesh.Entities(0)) {
    if (dofhandler.NumLocalDofs(*triangle) != 3) {
      throw "Not a S_1^0 FE Space!";
    }
    // Get area of triangle
    double K = lf::geometry::Volume(*triangle->Geometry());

    // Get nodes of triangle
    nonstd::span<const lf::mesh::Entity *const> nodes = triangle->SubEntities(2);
    // Loop over nodes of a i-th triangle
    for (const lf::mesh::Entity *node : nodes) {
      int id = mesh.Index(*node);
      I += (1.0/3.0)*mu(id)*K;
    }
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

    std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  for (const auto *cell : mesh->Entities(0)) {
    // check if we the FE space is really $\Cs_2^0$
    if (dofhandler.NumLocalDofs(*cell) != 6) {
      throw "Not a S_2^0 FE space!";
    }
    const double weight = 1.0 / 3.0 * lf::geometry::Volume(*(cell->Geometry()));
    // iterate over dofs
    auto int_dofs = dofhandler.GlobalDofIndices(*cell);
    // The integrated basis functions associated with the nodes are 0:
    //  $\int_K b_K^j dx = 0$ for $j = 1,2,3$. Skip!
    for (int l = 3; l < 6; ++l) {
      // The integrated basis functions associated with the edges are:
      // $\int_K b_K^j dx = |K|/3$ for $j = 4,5,6$
      // multiply by the value at the dof to get local contribution
      I += (weight * mu(int_dofs[l]));
    }
  }
  
  //====================
  return I;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler &dofh_Linear_FE,
    const lf::assemble::DofHandler &dofh_Quadratic_FE,
    const Eigen::VectorXd &mu) {
  if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
    throw "Underlying meshes must be the same for both DOF handlers!";
  }
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      dofh_Linear_FE.Mesh();                          // get the mesh
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NumDofs());  // initialise empty zeta
  // safety guard: always set zero if you're not sure to set every entry later
  // on for us this shouldn't be a problem, but just to be sure
  zeta.setZero();

  for (const auto *cell : mesh->Entities(0)) {
    // check if the spaces are actually linear and quadratic
    //====================
    // Your code goes here
    //====================
    // get the global dof indices of the linear and quadratic FE spaces, note
    // that the vectors obey the LehrFEM++ numbering, which we will make use of
    // lin\_dofs will have size 3 for the 3 dofs on the nodes and
    // quad\_dofs will have size 6, the first 3 entries being the nodes and
    // the last 3 the edges
    //====================
    // Your code goes here
    // assign the coefficients of mu to the correct entries of zeta, use
    // the previous subproblem 2-9.a
    //====================
  }
  return zeta;
}
/* SAM_LISTING_END_5 */

}  // namespace LFPPDofHandling
