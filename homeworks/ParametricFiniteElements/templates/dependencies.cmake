# Add your custom dependencies here:

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/parametricfiniteelements.h
  ${DIR}/parametricfiniteelements_main.cc
)

set(LIBRARIES
  Eigen3::Eigen
)
