# Check if cmake has the required version
cmake_minimum_required(VERSION 2.6.4)

if(CMAKE_VERSION VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
endif()

project(SMPT)

set(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp")

# Load Macro to locate ROOT and needed libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
find_package(ROOT REQUIRED Core RIO Net Hist Graf Graf3d Gpad Tree Rint
Postscript Matrix Physics MathCore Thread Gui TreePlayer Minuit TMVA)

#Check the compiler and set the compile and link flags
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)


include_directories(${ROOT_INCLUDE_DIRS} ${SMPT_SOURCE_DIR}/include)

# Build one shared lib with all sources
add_library(cropsimpletools SHARED ${SMPT_SOURCE_DIR}/src/cropcutspace.cc
${SMPT_SOURCE_DIR}/src/cropdataset.cc
${SMPT_SOURCE_DIR}/src/cropvarensemble.cc
${SMPT_SOURCE_DIR}/src/cropoptimisationengine.cc
${SMPT_SOURCE_DIR}/src/cropcutensemble.cc
${SMPT_SOURCE_DIR}/src/cropdatastore.cc
)


set(EXELIST thresholdeffs bwdiv crop stackergen rangefinder sepper mergevars updatedatastore stacker varstocuts eff corr cuttester columnmaker tuplesampler tuplescrambler cutapplier)

foreach(exe ${EXELIST})
  add_executable(${exe} ${SMPT_SOURCE_DIR}/src/${exe}.cc)
  target_link_libraries (${exe} cropsimpletools ${ROOT_LIBRARIES}) 
endforeach()


install(TARGETS cropsimpletools ${EXELIST}
  RUNTIME DESTINATION ${SMPT_SOURCE_DIR}/bin
  LIBRARY DESTINATION ${SMPT_SOURCE_DIR}/lib
  ARCHIVE DESTINATION ${SMPT_SOURCE_DIR}/lib/static)
