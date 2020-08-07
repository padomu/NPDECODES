/**
 * @ file LinearFE1D_main.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher, Erick Schulz
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

/* This problem does not rely on Lehrfempp. It is solely based on the
library Eigen. As the computational domain is simply the 1D interval
(0,1), the mesh is stored as a vector whose entries are the
coordinates of the nodes (i.e. the distance of the node from origin 0.0).
*/

#include "linearfe1d.h"
#include <cmath>

int main() {
  // Constant and variable parameters
  // Creating a 1D mesh of the interval (0,1)
  int N = 9;                  // nb. of cells
  auto f = [](double x) { return 0.0; };
  auto alpha = [](double x) { return 1.0; };
  double u0 = 0.0;
  double u1 = 0.0;


  Eigen::VectorXd mesh(N + 1);  // nb. of nodes
  // Nodes are equally spaced
  for (int i = 0; i < N + 1; i++) {
    mesh[i] = i * (1.0 / N);
  }
 


  // Solving the BVPs
  Eigen::VectorXd uA, uB, uC;
  //uA = LinearFE1D::solveA(mesh, identity, const_one);
  //uB = LinearFE1D::solveB(mesh, identity, const_one, 0.1, 0.5);
                               // alpha      f       u0    u1
  uB = LinearFE1D::solveB(mesh, alpha, f, u0, u1);

std::cout << "\n\n\nSolution u:\n" << uB << "\n\n\n\n\n";


  //uC = LinearFE1D::solveC(mesh, identity, identity);

  // PRINTING RESULTS TO.csv FILE
  // Defining CSV output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream ofs;
  // Printing results to file for problem (A)
  std::string filename = CURRENT_BINARY_DIR "/uA.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uA.format(CSVFormat);
    std::cout << "Components of uA were written to uA.csv" << std::endl;
  }
  ofs.close();
  if (ofs.is_open()) {
    std::cout << "File uA.csv was not properly closed." << std::endl;
  }
  // Printing results to file for problem (B)
  filename = CURRENT_BINARY_DIR "/uB.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uB.format(CSVFormat);
    std::cout << "Components of uB were written to uB.csv" << std::endl;
  }
  ofs.close();
  if (ofs.is_open()) {
    std::cout << "File uB.csv was not properly closed." << std::endl;
  }
  // Printing results to file for problem (C)
  filename = CURRENT_BINARY_DIR "/uC.csv";
  ofs.open(filename.c_str());
  if (ofs.is_open()) {
    ofs << uC.format(CSVFormat);
    std::cout << "Components of uC were written to uC.csv" << std::endl;
  }
  ofs.close();
}
