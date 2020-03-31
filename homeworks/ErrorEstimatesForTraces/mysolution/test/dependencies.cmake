# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/test/errorestimatesfortraces_test.cc
)

set(LIBRARIES
  Eigen3::Eigen
  GTest::gtest_main
  LF::lf.mesh.utils
  LF::lf.mesh.hybrid2d
  LF::lf.io
)
