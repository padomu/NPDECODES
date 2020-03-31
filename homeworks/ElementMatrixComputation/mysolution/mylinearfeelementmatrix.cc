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

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_1 */
<<<<<<< HEAD
Eigen::Matrix<double, 4, 4> MyLinearFEElementMatrix::Eval(
  const lf::mesh::Entity &cell) {
=======
Eigen::Matrix<double, 4, 4>
MyLinearFEElementMatrix::Eval(const lf::mesh::Entity &cell) {
>>>>>>> 3a2043cea445a57ecb3d77f1b8ed46cd255a6ade
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

  // a(u,v) = int grad(u)grad(v) + uv dx can be split into
  // the first term [Laplacian] and the second term [Mass].

  // Laplacian-Term - we use the built in Elemental Matrix Builder
  lf::uscalfe::LinearFELaplaceElementMatrix laplace_elmat_builder;
  auto laplace_elem_mat = laplace_elmat_builder.Eval(cell);

  // Mass-Term - We use our own computated element matrix
  // Note that we have to distinguish between
  // rectangular and triangular cells
  
//std::cout << "\n\n\n\n" << vertices << "\n\n\n\n";

  Eigen::Matrix<double, 4, 4> mass_elem_mat;

  if( ref_el == lf::base::RefEl::kTria() ) {
    mass_elem_mat << 2.0, 1.0, 1.0, 0.0,
                     1.0, 2.0, 1.0, 0.0,
                     1.0, 1.0, 2.0, 0.0,
                     0.0, 0.0, 0.0, 0.0;
    /*double K = 0.5
     * ( (vertices(0,1) - vertices(0,0)) * (vertices(1,2) - vertices(1,0))
     - (vertices(1,1) - vertices(1,0)) * (vertices(0,2) - vertices(0,0)) );
    */
    double K = lf::geometry::Volume(*geo_ptr);
    mass_elem_mat *= (K / 12.0 );
  } else if ( ref_el == lf::base::RefEl::kQuad() ) {
    mass_elem_mat << 4.0, 2.0, 1.0, 2.0,
                     2.0, 4.0, 2.0, 1.0,
                     1.0, 2.0, 4.0, 2.0,
                     2.0, 1.0, 2.0, 4.0;
    /*double K = (vertices(0,1) - vertices(0,0))*(vertices(1,3)-vertices(1,0));*/
    double K = lf::geometry::Volume(*geo_ptr);
    mass_elem_mat *= ( K / 36.0 );
  } else {
    LF_ASSERT_MSG(false, "Illegal cell type");
  }

  elem_mat = laplace_elem_mat + mass_elem_mat;

  //====================

  return elem_mat;
}
/* SAM_LISTING_END_1 */
} // namespace ElementMatrixComputation
