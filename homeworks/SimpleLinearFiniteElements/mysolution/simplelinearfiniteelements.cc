/**
 * @file
 * @brief NPDE homework SimpleLinearFiniteElements
 * @author Amélie Loher
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include "simplelinearfiniteelements.h"
#include <iostream>

namespace SimpleLinearFiniteElements {

double getArea(const Eigen::Matrix<double, 2, 3> &triangle) {
  return std::abs(
      0.5 *
      ((triangle(0, 1) - triangle(0, 0)) * (triangle(1, 2) - triangle(1, 1)) -
       (triangle(0, 2) - triangle(0, 1)) * (triangle(1, 1) - triangle(1, 0))));
}

Eigen::Matrix<double, 2, 3>
gradbarycoordinates(const Eigen::Matrix<double, 2, 3> &triangle) {
  Eigen::Matrix3d X;

  // solve for the coefficients of the barycentric coordinate functions, see
  // \eqref{eq:lambdalse}
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = triangle.transpose();
  return X.inverse().block<2, 3>(1, 0);
}

/**
 *  @brief Computation of Element Matrix for the Laplacian
 */
Eigen::Matrix3d
ElementMatrix_Lapl_LFE(const Eigen::Matrix<double, 2, 3> &triangle) {
  Eigen::Matrix<double, 2, 3> X = gradbarycoordinates(triangle);
  // compute inner products of gradients through matrix multiplication
  return getArea(triangle) * X.transpose() * X;
}

/**
 *  @brief Computation of full Galerkin Matrix
 */
Eigen::Matrix3d
ElementMatrix_LaplMass_LFE(const Eigen::Matrix<double, 2, 3> &triangle) {
  return ElementMatrix_Lapl_LFE(triangle) + ElementMatrix_Mass_LFE(triangle);
}

/**
 *  @brief Computation of element mass matrix on planar triangle
 *  @param triangle 2x3 matrix of vertex coordinates
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix3d
ElementMatrix_Mass_LFE(const Eigen::Matrix<double, 2, 3> &triangle) {
  Eigen::Matrix3d element_matrix;
  //====================
  double K = getArea(triangle);
  element_matrix << 2.0, 1.0, 1.0,
                    1.0, 2.0, 1.0,
                    1.0, 1.0, 2.0;

  element_matrix *= K/12.0;

  //====================
  return element_matrix;
}
/* SAM_LISTING_END_1 */

/**
 * @brief L2Error Computes the L2 error between the approximate solution and
 *                the exact solution
 * @param mesh the mesh to use
 * @param uFEM the solution approximated through FEM
 * @param exact the exact solution
 * @return the L2 difference
 */
/* SAM_LISTING_BEGIN_2 */
double L2Error(const TriaMesh2D &mesh, const Eigen::VectorXd &uFEM,
               const std::function<double(const Eigen::Vector2d &)> exact) {
  double l2error_squared = 0.0;
  //====================
  
  int N_dim = mesh.elements.rows();

  for(int i = 0; i<N_dim; ++i) {
    // Get triangle
    Eigen::Vector3i triangle_node_indexes = mesh.elements.row(i);

    // Get vertices
    Eigen::Vector2d v1 = mesh.vertices.row(triangle_node_indexes(0));
    Eigen::Vector2d v2 = mesh.vertices.row(triangle_node_indexes(1));
    Eigen::Vector2d v3 = mesh.vertices.row(triangle_node_indexes(2));

    // Save triangle
    Eigen::Matrix<double, 2, 3> triangle;
    //triangle.col(0) = v1;
    //triangle.col(1) = v2;
    //riangle.col(2) = v3;
    triangle << v1, v2, v3;
    // get area
    double K = getArea(triangle);

    double tmp = 0.0;
    tmp = ( uFEM(triangle_node_indexes(0)) - exact(v1) ) * ( uFEM(triangle_node_indexes(0)) - exact(v1) );
    tmp = ( uFEM(triangle_node_indexes(1)) - exact(v2) ) * ( uFEM(triangle_node_indexes(1)) - exact(v2) );
    tmp = ( uFEM(triangle_node_indexes(2)) - exact(v3) ) * ( uFEM(triangle_node_indexes(2)) - exact(v3) );
    tmp *= K/3.0;

    l2error_squared += tmp;
  }

  //====================

  return std::sqrt(l2error_squared);
}
/* SAM_LISTING_END_2 */

/**
 * @brief H1Serror Computes the H^1 error between the approximate solution and
 *                the exact solution
 * @param mesh the mesh to use
 * @param uFEM the solution approximated through FEM
 * @param exact the exact gradient of the solution
 * @return the H^1 difference
 *
 * @note This implementation seems to be flawed!
 */
/* SAM_LISTING_BEGIN_3 */
double
H1Serror(const TriaMesh2D &mesh, const Eigen::VectorXd &uFEM,
         const std::function<Eigen::Vector2d(const Eigen::Vector2d &)> exact) {
  double H1Serror_squared = 0.0;
  //====================
  
  int N_dim = mesh.elements.rows();

  for(int i = 0; i<N_dim; ++i) {
    // Get triangle
    Eigen::Vector3i triangle_node_indexes = mesh.elements.row(i);

    // Get vertices
    Eigen::Vector2d v1 = mesh.vertices.row(triangle_node_indexes(0));
    Eigen::Vector2d v2 = mesh.vertices.row(triangle_node_indexes(1));
    Eigen::Vector2d v3 = mesh.vertices.row(triangle_node_indexes(2));

    // Save triangle
    Eigen::Matrix<double, 2, 3> triangle;
    //triangle.col(0) = v1;
    //triangle.col(1) = v2;
    //riangle.col(2) = v3;
    triangle << v1, v2, v3;
    // get area
    double K = getArea(triangle);

    Eigen::Matrix<double, 2, 3> gradTriangle = gradbarycoordinates(triangle);

    // Get gradients of barycentric coord. functions (the lambdas) for the given triangle
    // We need those to compute grad_uh
    Eigen::Matrix<double, 2, 3> gradBary = gradbarycoordinates(triangle);

    // Get grad(u_h(a^i))
    // Note that grad(u_h(a^1))=grad(u_h(a^2))=grad(u_h(a^3)) since we are on a plane.
    // A plane has only one gradient.. since... you know ... it's a plane
    Eigen::Vector2d grad_uh = uFEM(triangle_node_indexes(0))*gradBary.col(0)
                            + uFEM(triangle_node_indexes(1))*gradBary.col(1)
                            + uFEM(triangle_node_indexes(2))*gradBary.col(2);

    double tmp = 0.0;

    tmp = (grad_uh - exact(v1)).squaredNorm();
    tmp += (grad_uh - exact(v2)).squaredNorm();
    tmp += (grad_uh - exact(v2)).squaredNorm();

    //tmp = ( grad_uh - exact(triangle_node_indexes(0)) ).squaredNorm();

    tmp *= K/3.0;

    H1Serror_squared += tmp;
  }
  //====================

  return std::sqrt(H1Serror_squared);
}
/* SAM_LISTING_END_3 */

/**
 * @brief assemLoad_LFE Assembles the Load Vector
 * @param mesh the mesh to use
 * @param getElementVector
 * @param f function handle for f
 * @return assembled load vector
 */
Eigen::VectorXd
assemLoad_LFE(const TriaMesh2D &mesh,
              const std::function<double(const Eigen::Vector2d &)> &f) {
  // obtain the number of triangles
  int M = mesh.elements.rows();

  // obtain the number of vertices
  int N = mesh.vertices.rows();
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N);

  // loop over all triangles
  for (int i = 0; i < M; i++) {
    Eigen::Matrix<double, 2, 3> triangle = mesh[i];

    // loop over vertices of current triangle
    double factor = getArea(triangle) / 3.0;
    for (int j = 0; j < 3; ++j) {
      // from local to global load vector
      phi(mesh.elements(i, j)) += factor * f(triangle.col(j));
    }
  }

  return phi;
}

/**
 * @brief GalerkinAssembly Assembles the Galerkin Matrix
 * @param mesh the mesh to use
 * @param getElementMatrix Element Matrix
 * @return Galerkin Matrix
 */
Eigen::SparseMatrix<double> GalerkinAssembly(
    const TriaMesh2D &mesh,
    const std::function<Eigen::Matrix3d(const Eigen::Matrix<double, 2, 3> &)>
        &getElementMatrix) {

  // obtain the number of vertices
  int N = mesh.vertices.rows();
  // obtain the number of elements/cells
  int M = mesh.elements.rows();
  std::vector<Eigen::Triplet<double>> triplets;
  // loop over elements and add local contributions
  for (int i = 0; i < M; i++) {
    // get local$\to$global index mapping for current element, \emph{cf.}
    // \lref{eq:idxdef}
    Eigen::Vector3i element = mesh.elements.row(i);
    Eigen::Matrix<double, 2, 3> triangle;
    // extract vertices of current element
    for (int j = 0; j < 3; j++) {
      triangle.col(j) = mesh.vertices.row(element(j)).transpose();
    }
    // compute element contributions
    Eigen::Matrix3d Ak = getElementMatrix(triangle);
    // build triplets from contributions
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        triplets.push_back({element(j), element(k), Ak(j, k)});
      }
    }
  }
  // build sparse matrix from triplets
  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  return A;
}

/**
 * @brief solves system and prints H1-semierror, L2 error, the mesh and a
 * surface plot
 * @param mesh: discretisation of the computational domain
 */
/* SAM_LISTING_BEGIN_4 */
std::tuple<Eigen::VectorXd, double, double>
Solve(const SimpleLinearFiniteElements::TriaMesh2D &mesh) {
  const double pi = 3.1415926535897;

  // define the source function f
  auto f = [pi](const Eigen::Vector2d &x) {
    return (1.0 + 8.0 * pi * pi) * std::cos(2.0 * pi * x(0)) *
           std::cos(2.0 * pi * x(1));
  };
  // the exact solution of the linear variational problem
  auto uExact = [pi](const Eigen::Vector2d &x) {
    return std::cos(2 * pi * x(0)) * std::cos(2 * pi * x(1));
  };

  Eigen::VectorXd U;
  double l2error;
  double h1error;

  //====================
  // Your code goes here

  auto gradUExact = [pi] (const Eigen::Vector2d& x) {
    Eigen::Vector2d gradient;
    gradient << -2.0*pi*std::sin(2*pi*x(0))*std::cos(2.0*pi*x(1)),
                -2.0*pi*std::cos(2*pi*x(0))*std::sin(2.0*pi*x(1));

    return gradient;
  };

  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::GalerkinAssembly(mesh, SimpleLinearFiniteElements::ElementMatrix_LaplMass_LFE);
  Eigen::VectorXd rhs_vector = SimpleLinearFiniteElements::assemLoad_LFE(mesh, f);

  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
  //solver.compute(A);
  solver.analyzePattern(A);
  solver.factorize(A);

  U = solver.solve(rhs_vector);




  l2error = SimpleLinearFiniteElements::L2Error(mesh, U, uExact);
  h1error = SimpleLinearFiniteElements::H1Serror(mesh, U, gradUExact);

  //====================
  return std::make_tuple(U, l2error, h1error);
}
/* SAM_LISTING_END_4 */

} // namespace SimpleLinearFiniteElements
