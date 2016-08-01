# setting up of my project
cmake_minimum_required(VERSION 3.0)
if (POLICY CMP0042)
	cmake_policy (SET CMP0042 OLD)
endif()
project(G4Hodoscope_vega)

set(CMAKE_BUILD_TYPE debug)
set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_CXX_FLAGS "-w -fPIC -O0")
add_definitions(-std=c++11)

# Find Geant4 package, activating all available UI and Vis drivers
option(WITH_GEANT4_UIVIS "Build Tomography with UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()

# Setup geant4 include directories and compile definitions
# Setup include directory for this project
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable and link to geant4 Libraries
add_executable(muonHodo muonHodo.cc ${sources} ${headers})
target_link_libraries(muonHodo ${Geant4_LIBRARIES})



# Install the executable to bin directory under CMAKE_INSTALL)PREFIX
install(TARGETS muonHodo DESTINATION bin)
