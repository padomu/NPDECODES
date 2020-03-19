/**
 * @ file LinearFE1D.h
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher, Erick Schulz
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/SparseCholesky>

namespace LinearFE1D {

// Calculate the matrix entries corresponding to the Galerkin matrix
// of the Laplace operator with non-constant coefficient function alpha
// using composite midpoint integration rule. The entries are stored in
// the Eigen triplet data structure.
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR1>
std::vector<Eigen::Triplet<double>> computeA(const Eigen::VectorXd& mesh,
                                             FUNCTOR1&& alpha) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializing the vector of triplets whose size corresponds to the
  // number of entries in a (N+1) x (N+1) tridiagonal (band) matrix
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(3 * N + 1);

  //====================
  
  // NOTE: This is a general Galerking Matrix. For Problems with vanish-
  //       boundaries and non-vanishing bundaries. Which means, we
  //       include b_0 and b_M => so we have (N+1) basis-vectors.
  //    !! So A is a (N+1)x(N+1) Matrix.
//
  // A(0,0)
  triplets.push_back( Eigen::Triplet<double>(0,0, 0) );

  double value = -1.0/ ( mesh(1) - mesh(0) );
  // A(1,0)
  triplets.push_back( Eigen::Triplet<double>(1,0, value) );
  // A(0,1)
  triplets.push_back( Eigen::Triplet<double>(0,1, value) );

  // A(i,i) - Diagonal
  for(size_t i=1; i<N-1; ++i) {
    value = 1.0 / ( mesh(i+1) - mesh(i)   )
                 + 1.0 / ( mesh(i+2) - mesh(i+1) );
    triplets.push_back( Eigen::Triplet<double>(i,i, value) );
  }
  // A(i-1,i) & A(i,i-1) - Sub-Diagonals
  for(size_t i=1; i<N-1; ++i) {
    value = -1.0 / ( mesh(i+1) - mesh(i) );
    triplets.push_back( Eigen::Triplet<double>(i,i+1, value) );
    triplets.push_back( Eigen::Triplet<double>(i+1,i, value) );
  }
  // A(N,N)
  triplets.push_back( Eigen::Triplet<double>(0,0, 0) );

  value = -1.0/ ( mesh(N) - mesh(N-1) );
  // A(N,N-1)
  triplets.push_back( Eigen::Triplet<double>(N,N-1, value) );
  // A(N-1,N)
  triplets.push_back( Eigen::Triplet<double>(N-1,N, value) );

  //====================

  return triplets;
}  // computeA
/* SAM_LISTING_END_1 */

// Calculate the matrix entries corresponding to the mass matrix
// associated to the L2 pairing weighted by the variable coefficient
// function gamma using the trapezoidal integration rule.
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR1>
std::vector<Eigen::Triplet<double>> computeM(const Eigen::VectorXd& mesh,
                                             FUNCTOR1&& gamma) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;

  // The basis tent functions spanning the space of continuous piecewise
  // linear polynomial satisfy the cardinal basis property with respect
  // to the nodes of the mesh. Using  the trapezoidal rule to approximate
  // the integrals therefore lead to a diagonal (N+1) x (N+1) matrix
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(N + 1);

  //====================

  // M(0,0)
  double value = mesh(1) - mesh(0);
  value *= 0.5*gamma(mesh(0));
  triplets.push_back( Eigen::Triplet<double>(0,0, value) );
  // M(N,n)
  value = mesh(N) - mesh(N-1);
  value *= 0.5*gamma(mesh(N));
  triplets.push_back( Eigen::Triplet<double>(N,N, value) );

  // Diaognal entries
  for(size_t i=1; i<N; ++i) {
    value = mesh(i+1) - mesh(i-1);
    value *= 0.5*gamma(mesh(i));

    triplets.push_back( Eigen::Triplet<double>(i,i, value) );
  }

  //====================

  return triplets;
}  // computeM
/* SAM_LISTING_END_2 */

// Calculate the entries of the right hand side vector using
// the composite trapezoidal integration rule
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR1>
Eigen::VectorXd computeRHS(const Eigen::VectorXd& mesh, FUNCTOR1&& f) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializing right hand side vector
  Eigen::VectorXd rhs_vec = Eigen::VectorXd::Zero(N + 1);

  // Some tool variables
  double dx;

  /* Calculate the entries */
  // First entry (left boundary node)
  rhs_vec(0) = f(mesh(0)) * (mesh(1) - mesh(0)) * 0.5;
  // Last entry (right boundary node)
  rhs_vec(N) = f(mesh(N)) * (mesh(N) - mesh(N - 1)) * 0.5;
  // All other enties associated to interior nodes
  for (unsigned i = 1; i < N; ++i) {
    dx = mesh(i + 1) - mesh(i - 1);
    rhs_vec(i) = f(mesh(i)) * 0.5 * dx;
  }

  return rhs_vec;
}   // computeRHS
/* SAM_LISTING_END_3 */
  
// SOLVE THE LINEAR SYSTEM OF PROBLEM (A)
/* SAM_LISTING_BEGIN_A */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveA(const Eigen::VectorXd& mesh, FUNCTOR1&& gamma,
                       FUNCTOR2&& f) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializations (notice initialization with zeros here)
  Eigen::VectorXd u = Eigen::VectorXd::Zero(N + 1);  // solution vec
  Eigen::SparseMatrix<double> A(N + 1, N + 1);       // laplacian galerkin mat
  Eigen::SparseMatrix<double> M(N + 1, N + 1);       // mass galerkin mat
  Eigen::SparseMatrix<double> L(N + 1, N + 1);       // full galerkin mat

  // I. Build the (full) Galerkin matrix L for the lin. sys.
  // I.i Compute the entries of the Laplace Galerkin matrix A
  // (const. ceoff. func. alpha = 1.0)
  std::vector<Eigen::Triplet<double>> triplets_A =
      computeA(mesh, [](double) -> double { return 1.0; });
  // I.ii Compute the entries of the mass matrix M
  std::vector<Eigen::Triplet<double>> triplets_M = computeM(mesh, gamma);
  // I.iii Assemble the sparse matrices
  A.setFromTriplets(triplets_A.begin(), triplets_A.end());
  M.setFromTriplets(triplets_M.begin(), triplets_M.end());
  L = A + M;  // Full Galerkin matrix of the LSE

  // II. Build the right hand side source vector
  Eigen::VectorXd rhs_vec = computeRHS(mesh, f);

  // III. Enforce (zero) dirichlet boundary conditions
  // The B.C. are u(0) = u(N) = 0. The boundary nodes are indexed by 0
  // and N. We can therefore opt for the option of droping first and last rows
  // and columns of L, along with the first and last entries of the rhs_vec.
  Eigen::SparseMatrix<double> L_reduced = L.block(1, 1, N - 1, N - 1);
  Eigen::VectorXd rhs_vec_reduced = rhs_vec.segment(1, N - 1);

  // IV. Solve the LSE L*u = rhs_vec using an Eigen solver
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
  solver.compute(L_reduced);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  u.segment(1, N - 1) = solver.solve(rhs_vec_reduced);

  // The solution vector u was initialized with zeros, and therefore already
  // contains the zero Dirichlet boundary data in the first and last entry
  return u;
}  // solveA
/* SAM_LISTING_END_A */

// SOLVE THE LINEAR SYSTEM OF PROBLEM (B)
/* SAM_LISTING_BEGIN_B */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveB(const Eigen::VectorXd& mesh, FUNCTOR1&& alpha,
                       FUNCTOR2&& f, double u0, double u1) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializations
  Eigen::VectorXd u(N + 1);                     // solution vec
  Eigen::SparseMatrix<double> A(N + 1, N + 1);  // laplacian galerkin mat
  // Some tool variables
  double dx_left, dx_right;  // cell widths

  // I. Build the Galerkin matrix A
  //====================
  // Your code goes here
  //====================

  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  // III. Enforce dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE A*u = rhs_vec using an Eigen solver
  //====================
  // Your code goes here
  //====================

  // The solution vector still needs to be supplemented with the known
  // boundary values
  u(0) = u0;  // left boundary node
  u(N) = u1;  // right boundary node
  return u;
}  // solveB
/* SAM_LISTING_END_B */

// Build an sol!ve the LSE corresponding to (C)
/* SAM_LISTING_BEGIN_C */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveC(const Eigen::VectorXd& mesh, FUNCTOR1&& alpha,
                       FUNCTOR2&& gamma) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializations (notice initialization with zeros here)
  Eigen::VectorXd u(N + 1);                     // solution vec
  Eigen::SparseMatrix<double> A(N + 1, N + 1);  // laplacian galerkin mat
  Eigen::SparseMatrix<double> M(N + 1, N + 1);  // mass galerkin mat
  Eigen::SparseMatrix<double> L(N + 1, N + 1);  // full galerkin mat

  // I. Build the (full) Galerkin matrix L for the lin. sys.
  //====================
  // Your code goes here
  //====================

  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE A*u = rhs_vec using an Eigen solver
  //====================
  // Your code goes here
  //====================

  return u;
}  // solveC
/* SAM_LISTING_END_C */

}  // namespace LinearFE1D
