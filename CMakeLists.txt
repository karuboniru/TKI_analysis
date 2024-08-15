cmake_minimum_required(VERSION 3.12)
project(MyProject)

set(CMAKE_CXX_STANDARD 23)
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
add_link_options(-Wl,--as-needed)

# Locate CERN ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Locate NuWro
if (DEFINED NUWRO)
  set(NUWRO_DIR ${NUWRO})
else()
  set(NUWRO_DIR $ENV{NUWRO})
endif()

if (NOT NUWRO_DIR)
  message(FATAL_ERROR "Please set NUWRO environment variable to the path to the NuWro installation")
endif()

add_library(EVENT1 SHARED IMPORTED)
set_property(TARGET EVENT1 PROPERTY IMPORTED_LOCATION ${NUWRO_DIR}/bin/event1.so)

# Include ROOT and NuWro headers
include_directories(${ROOT_INCLUDE_DIRS} ${NUWRO_DIR}/src include)

# Find all .cxx files in the src directory
file(GLOB SOURCES src/*.cxx)
file(GLOB HEADERS include/*.h)

#  ${SOURCES} -> libanalysis.so
add_library(analysis SHARED ${SOURCES})
target_link_libraries(analysis PRIVATE EVENT1 ROOT::Physics ROOT::ROOTDataFrame)
root_generate_dictionary(G__Dict ${HEADERS} MODULE analysis LINKDEF LinkDef.h)

# Add the main executable
add_executable(T2K T2K.cpp)

# Link all the object files to the main executable
target_link_libraries(T2K PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)

add_executable(MINERvA_pi0 MINERvA_pi0.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_pi0 PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)

add_executable(MINERvA_0pi MINERvA_0pi.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_0pi PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)


add_executable(rawtheta rawtheta.cpp)

# Link all the object files to the main executable
target_link_libraries(rawtheta PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)

add_executable(MINERvA_all MINERvA_all.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_all PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)


add_executable(MINERvA_all_gibuu MINERvA_all_gibuu.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_all_gibuu PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)


add_executable(gibuu2root gibuu2root.cpp)

# Link all the object files to the main executable
target_link_libraries(gibuu2root PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
)


add_executable(rawq2 rawq2.cpp)

# Link all the object files to the main executable
target_link_libraries(rawq2 PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)

add_executable(rawq2_gibuu rawq2_gibuu.cpp)

# Link all the object files to the main executable
target_link_libraries(rawq2_gibuu PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)


add_executable(minerva_all_gibuu_plot minerva_all_gibuu_plot.cpp)

# Link all the object files to the main executable
target_link_libraries(minerva_all_gibuu_plot PRIVATE 
  ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree ROOT::Hist ROOT::Gpad
  EVENT1
  analysis
)
