/**
 * @file  maximumprinciple.cc
 * @brief NPDE homework "MaximumPrinciple" code
 * @author Oliver Rietmann
 * @date 25.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "maximumprinciple.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

namespace MaximumPrinciple {

/**
 * @brief Assembly on a tensor product mesh
 *
 * Compute the global Galerkin matrix from the local
 * element matrix.
 *
 * @param M Number of interior vertices in x and y direction.
 * @param B_K Local element matrix.
 * @return Global Galerkin matrix of size M^2 times M^2.
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> assemble(int M, const Eigen::Matrix3d &B_K) {
  int M2 = M * M;
  Eigen::SparseMatrix<double> A(M2, M2);
  //====================
  
  // allocate memory
  A.reserve(7);

  // Get |K| - note that each cell has the same shape, so we can
  // get it outside of the loop
  double h = 1.0/(M+1);
  double K = 0.5*h*h;

  // Looper over all nodes
  for(int i=0; i<M2; ++i) {

    A.coeffRef(i,j) = v;
  }

  //====================
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<double> computeGalerkinMatrix(int M, double c) {
  Eigen::Matrix3d B_K;
  //====================
  double h = 1./(M+1);
  double K = 0.5*h*h;

  Eigen::Matrix3d A_K;
  A_K << 2.0, -1.0, -1.0,
        -1.0,  2.0,  0.0
        -1.0,  0.0,  2.0;
  A_K *= 0.5;

  Eigen::Matrix3d M_K;
  M_K << 2.0, 1.0, 1.0,
         1.0, 2.0, 1.0,
         1.0, 1.0, 2.0;
  M_K *= K/12.0;

  B_K = (1.0 - c)*A_K + c*M_K;
  //====================
  return assemble(M, B_K);
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
Eigen::SparseMatrix<double> computeGalerkinMatrixTR(int M, double c) {
  Eigen::Matrix3d B_K;
  //====================
  // Your code goes here
  //====================
  return assemble(M, B_K);
}
/* SAM_LISTING_END_4 */

}  // namespace MaximumPrinciple
