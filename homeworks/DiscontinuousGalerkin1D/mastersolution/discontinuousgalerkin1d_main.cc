/**
 * @file discontinuousgalerkin1d_main.cc
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "discontinuousgalerkin1d.h"

int main() {
  DiscontinuousGalerkin1D::Solution solution =
      DiscontinuousGalerkin1D::solveTrafficFlow();

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");

  std::ofstream file;
  file.open("solution.csv");
  file << solution.x_.transpose().format(CSVFormat) << std::endl;
  file << solution.u_.transpose().format(CSVFormat) << std::endl;
  file.close();

  std::cout << "Generated " CURRENT_BINARY_DIR "/solution.csv" << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR
              "/plot_solution.py " CURRENT_BINARY_DIR
              "/solution.csv " CURRENT_BINARY_DIR "/solution.eps");

  return 0;
}
