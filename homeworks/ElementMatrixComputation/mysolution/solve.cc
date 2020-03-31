/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 06.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "solve.h"

#include <iostream>

#include <Eigen/Core>

#include "mylinearfeelementmatrix.h"
#include "mylinearloadvector.h"

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solvePoissonBVP() {
  // Convert the globally defined function f to a LehrFEM++ mesh function object
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // The basis expansion coefficient vector for the finite-element solution
  Eigen::VectorXd solution = Eigen::VectorXd::Zero(1);

  // Elemental Matrix Builder
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder;

  // Elemental Vector Builder
  lf::uscalfe::LinearFELocalLoadVector<double, decltype(mf_f)> elvec_builder(mf_f);

  solution = solve(elmat_builder, elvec_builder);

  //====================

  return solution;
}
/* SAM_LISTING_END_2 */

Eigen::VectorXd solveNeumannEq() {
  // Define the solution vector
  Eigen::VectorXd solution;

  //====================
  // Your code goes here
  //====================

  return solution;
}

}  // namespace ElementMatrixComputation
