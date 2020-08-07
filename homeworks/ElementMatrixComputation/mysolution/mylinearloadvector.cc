/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearloadvector.h"

#include <functional>

#include <Eigen/Core>

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

namespace ElementMatrixComputation {

namespace {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector4d
computeLoadVector(const Eigen::MatrixXd &vertices,
                  std::function<double(const Eigen::Vector2d &)> f) {
  // Number of nodes of the element: triangles = 3, rectangles = 4
  const int num_nodes = vertices.cols();
  // Vector for returning element vector
  Eigen::Vector4d elem_vec = Eigen::Vector4d::Zero();

  //====================

  // first two midpoints of edges of cell
  Eigen::Vector2d m1 = 0.5*(vertices.col(0) + vertices.col(1));
  Eigen::Vector2d m2 = 0.5*(vertices.col(1) + vertices.col(2));

  // == TRIANGLE ==
  if (num_nodes == 3) {
    // Third midpoint of edge
    Eigen::Vector2d m3 = 0.5*(vertices.col(2) + vertices.col(0));
    
    // Area
    double K = 0.5*( (vertices(0,1) - vertices(0,0)) * (vertices(1,2) - vertices(1,0)) -
                     (vertices(1,1) - vertices(1,0)) * (vertices(0,2) - vertices(0,0)));
    K /= num_nodes;

    elem_vec(0) = f(m1) + f(m3);
    elem_vec(1) = f(m1) + f(m2);
    elem_vec(2) = f(m2) + f(m3);
    elem_vec *= K/2.0;

    // == QUADRILITERA ==
    } else if (num_nodes == 4) {
    // Third and forth midpoints of edge
    Eigen::Vector2d m3 = 0.5*(vertices.col(2) + vertices.col(3));
    Eigen::Vector2d m4 = 0.5*(vertices.col(3) + vertices.col(0));
    
    // Area
    double K = ( vertices(0,1) - vertices(0,0) )*
               ( vertices(1,3) - vertices(1,0) );
    K /= num_nodes;

    elem_vec(0) = f(m1) + f(m4);
    elem_vec(1) = f(m1) + f(m2);
    elem_vec(2) = f(m2) + f(m3);
    elem_vec(3) = f(m3) + f(m4);
    elem_vec *= K/2.0;
    }

  //====================

  return elem_vec;
}
/* SAM_LISTING_END_1 */

} // namespace

Eigen::Vector4d MyLinearLoadVector::Eval(const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();

  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());

  return computeLoadVector(vertices, f_);
}

} // namespace ElementMatrixComputation
