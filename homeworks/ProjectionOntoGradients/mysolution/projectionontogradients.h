/**
 * @file
 * @brief NPDE homework ProjectionOntoGradients code
 * @author ?, Philippe Peter
 * @date December 2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <utility>

namespace ProjectionOntoGradients {

/* SAM_LISTING_BEGIN_1 */
class ElementMatrixProvider {
 public:
  Eigen::Matrix3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity & /*entity*/) const { return true; }
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix3d ElementMatrixProvider::Eval(const lf::mesh::Entity &entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;

  //====================

  // Get area of entity
  const double K = lf::geometry::Volume(*geo_ptr);

  // Get coordiantes of corners of entity/triangle
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);

  // We know A_K=a(b_i,b_j) for i,j=1,2,3
  // We plug our a from 2.10.a int that and get
  // A_K=2*int_K grad(b_h^i)grad(b_h^j)dx fpr i,j=1,2,3
  // = 2int_k grad(lambda_i)grad(lambda_j)dx
  // We can use 2.4.5.10 to get the element matrix for this.
  // See code: 2.4.5.11 and 2.4.5.13

  // First Matrix (with a's in it)
  Eigen::Matrix3d gradients;
  gradients.col(0) = Eigen::Vector3d::Ones();
  gradients.rightCols(2) = corners.transpose();

  // Second Matrix (with the betas in it)
  // Solve it:
  const Eigen::MatrixXd betas = gradients.inverse().bottomRows(2);

  // Now we use code 2.4.5.13
  loc_mat = K * betas.transpose() * betas;

  //====================
  return loc_mat;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
class GradProjRhsProvider {
 public:
  explicit GradProjRhsProvider(FUNCTOR f) : f_(f) {}

  Eigen::Vector3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity & /*entity*/) const { return true; }

 private:
  FUNCTOR f_;
};
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename FUNCTOR>
Eigen::Vector3d GradProjRhsProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Vector3d loc_vec;

  //====================
  // Your code goes here
  //====================
  return loc_vec;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
template <typename FUNCTOR>
Eigen::VectorXd projectOntoGradients(const lf::assemble::DofHandler &dofh,
                                     FUNCTOR f) {
  const lf::assemble::size_type N_dofs = dofh.NumDofs();
  Eigen::VectorXd sol_vec;

  // I. Build the (full) Galerkin matrix for the linear system.
  //====================
  // Your code goes here
  //====================

  // II. Build the (full) right hand side vector
  //====================
  // Your code goes here
  //====================

  // III. Enforce homogeneous dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE using an Eigen solver
  //====================
  // Your code goes here
  //====================
  return sol_vec;
}
/* SAM_LISTING_END_5 */

}  // namespace ProjectionOntoGradients
