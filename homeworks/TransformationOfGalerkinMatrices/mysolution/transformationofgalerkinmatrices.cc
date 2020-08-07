/**
 * @file
 * @brief NPDE homework TransformationOfGalerkinMatrices code
 * @author Erick Schulz
 * @date 01/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "transformationofgalerkinmatrices.h"
#include <cassert>

namespace TransformationOfGalerkinMatrices {

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Triplet<double>> transformCOOmatrix(
    const std::vector<Eigen::Triplet<double>> &A) {
  std::vector<Eigen::Triplet<double>> A_t{};  // return value

  // First step: find the size of the matrix by searching the maximal
  // indices. Depends on the assumption that no zero rows/columns occur.
  int rows_max_idx = 0, cols_max_idx = 0;
  for (const Eigen::Triplet<double> &triplet : A) {
    rows_max_idx =
        (triplet.row() > rows_max_idx) ? triplet.row() : rows_max_idx;
    cols_max_idx =
        (triplet.col() > cols_max_idx) ? triplet.col() : cols_max_idx;
  }
  int n_rows = rows_max_idx + 1;
  int n_cols = cols_max_idx + 1;

  // Make sure we deal with a square matrix
  assert(n_rows == n_cols);
  // The matrix size must have even parity
  assert(n_cols % 2 == 0);

  int N = n_cols;      // Size of (square) matrix
  int M = n_cols / 2;  // Half the size
  //====================
  Eigen::SparseMatrix<double> AA(n_rows,n_cols);
  AA.setFromTriplets(A.begin(), A.end());
  Eigen::MatrixXd K = Eigen::MatrixXd(AA);
  Eigen::Triplet<double> el; //(i,j,value)
  double val = 0;
  for(size_t i=0; i<N; ++i) {
    for(size_t j=0; j<N; ++j) {
      // 1 <= i,j <= M
      if (i<M && j<M) {
        val = K(2*i,2*j)   + K(2*i,2*j-1)
            + K(2*i-1,2*j) + K(2*i-1,2*j-1);
        el = Eigen::Triplet<double>(i,j, val);
        A_t.push_back(el);
      }

      // M+1 <= i,j <= N
      if (i >= M && j < N) {
        val = K(2*(i-M)-1, 2*(j-M)-1) - K(2*(i-M), 2*(j-M)-1)
            - K(2*(i-M)-1, 2*(j-M))   + K(2*(i-M), 2*(j-M));
        el = Eigen::Triplet<double>(i,j, val);
        A_t.push_back(el);
      }

      // M+1 <= i <= N, 1 <= j <= M
      if (i >= M && j < M) {
        val = K(2*(i-N)-1, 2*j-1) - K(2*(i-M), 2*j-1)
            + K(2*(i-M)-1, 2*j)   - K(2*(i-M), 2*j);
        el = Eigen::Triplet<double>(i,j, val);
        A_t.push_back(el);
      }
      // 1<=i,j<=M
      if (i < M && j >= M ) {
        val = K(2*i-1, 2*(j-M)-1) + K(2*i, 2*(j-M)-1)
            - K(2*i-1,2*(j-M)) - K(2*i, 2*(j-M));
        el = Eigen::Triplet<double>(i,j, val);
        A_t.push_back(el);
      }

    }
  }

  //====================
  return A_t;
}
/* SAM_LISTING_END_1 */

}  // namespace TransformationOfGalerkinMatrices
