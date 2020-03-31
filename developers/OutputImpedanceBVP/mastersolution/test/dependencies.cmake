#if SOLUTION
# Dependencies of mastersolution tests:
#else
# Add your custom dependencies here:
#endif

# PROBLEM_NAME and DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/outputimpedancebvp_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.assemble
  LF::lf.io
  LF::lf.mesh.hybrid2d
  LF::lf.uscalfe
)
