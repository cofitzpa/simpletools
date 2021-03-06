# Check if cmake has the required version
cmake_minimum_required(VERSION 2.6.4)

if(CMAKE_VERSION VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
endif()

project(SMPT)

# Load Macro to locate ROOT and needed libraries
set(ROOT_REQUIRED_COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint
Postscript Matrix Physics MathCore Thread Gui TreePlayer Minuit TMVA)
#set(ROOT_OPTIONAL_COMPONENTS Proof)
find_package(ROOT REQUIRED ${ROOT_REQUIRED_COMPONENTS}
OPTIONAL_COMPONENTS ${ROOT_OPTIONAL_COMPONENTS})

# Set preprocessor flags for optional components
foreach(opt ${ROOT_OPTIONAL_COMPONENTS})
  if(ROOT_${opt}_LIBRARY)
    add_definitions(-Duse_${opt})
  endif()
endforeach()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

# Locate Boost
find_package(Boost REQUIRED system)

# Check the compiler and set the compile and link flags
#set(CMAKE_BUILD_TYPE Debug)
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
target_link_libraries(cropsimpletools ${ROOT_LIBRARIES})

# List executables and link each of them to the local and ROOT libraries
set(EXELIST thresholdeffs bwdiv crop stackergen rangefinder sepper mergevars updatedatastore stacker varstocuts eff corr cuttester columnmaker multicolumnmaker tuplesampler tuplescrambler cutapplier)

foreach(exe ${EXELIST})
  add_executable(${exe} ${SMPT_SOURCE_DIR}/src/${exe}.cc)
  target_link_libraries (${exe} cropsimpletools ${ROOT_LIBRARIES})
endforeach()

# Configure install destination
install(TARGETS cropsimpletools ${EXELIST}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static)

# Provide easy (but dumb) way to uninstall
add_custom_target(uninstall xargs rm -v < install_manifest.txt)

