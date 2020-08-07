/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearfeelementmatrix.h"

#include <Eigen/Core>

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <iostream>

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 4, 4>
MyLinearFEElementMatrix::Eval(const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());
  // Matrix for returning element matrix
  Eigen::Matrix<double, 4, 4> elem_mat;

  //====================

  // Get element matrix for negative laplacian of given cell
  lf::uscalfe::LinearFELaplaceElementMatrix prov;
  Eigen::MatrixXd A = prov.Eval(cell);

  // First get the area
  double K = lf::geometry::Volume(*geo_ptr);

  // Initialize elem. stiffness matrix
  Eigen::Matrix<double, 4, 4> M;

  if (ref_el == lf::base::RefEl::kTria()) {
    M << 2.0, 1.0, 1.0, 0.0,
         1.0, 2.0, 1.0, 0.0,
         1.0, 1.0, 2.0, 0.0,
         0.0, 0.0, 0.0, 0.0;
    M *= K/12.0;
  } else if (ref_el == lf::base::RefEl::kQuad()) {
    M << 4.0, 2.0, 1.0, 2.0,
         2.0, 4.0, 2.0, 1.0,
         1.0, 2.0, 4.0, 2.0,
         2.0, 1.0, 2.0, 4.0;
    M *= K/36.0;
  }

  elem_mat = A + M;

  //====================

  return elem_mat;
}
/* SAM_LISTING_END_1 */
} // namespace ElementMatrixComputation
