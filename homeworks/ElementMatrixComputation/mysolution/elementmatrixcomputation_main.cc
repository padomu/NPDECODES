/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "solve.h"

int main(int /*argc*/, const char ** /*argv*/) {
  ElementMatrixComputation::solvePoissonBVP();
  std::cout << "\n !!!!! \n";
  ElementMatrixComputation::solveNeumannEq();
  return 0;
}
