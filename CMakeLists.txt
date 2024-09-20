cmake_minimum_required(VERSION 3.12)
project(MyProject)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
add_link_options(-Wl,--as-needed)

# Locate Boost
find_package(Boost REQUIRED COMPONENTS system)

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})



# Locate CERN ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# enable O3
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native")

# Locate NuWro
if (DEFINED NUWRO)
  set(NUWRO_DIR ${NUWRO})
else()
  set(NUWRO_DIR $ENV{NUWRO})
endif()

if (NOT NUWRO_DIR)
  # message(FATAL_ERROR "Please set NUWRO environment variable to the path to the NuWro installation")
  set(NUWRO_DIR /var/home/yan/neutrino/nuwro/nuwro.nuwro_21.09.2)
endif()

add_library(EVENT1 SHARED IMPORTED)
set_property(TARGET EVENT1 PROPERTY IMPORTED_LOCATION ${NUWRO_DIR}/bin/event1.so)

# Include ROOT and NuWro headers
include_directories(${ROOT_INCLUDE_DIRS} ${NUWRO_DIR}/src include)

# Find all .cxx files in the src directory
file(GLOB SOURCES src/*.cxx)
file(GLOB HEADERS include/*.h)


#  ${SOURCES} -> libanalysis.so
add_library(analysis STATIC ${SOURCES})
target_link_libraries(analysis PRIVATE EVENT1 ROOT::Physics ROOT::ROOTDataFrame)
root_generate_dictionary(G__Dict ${HEADERS} MODULE analysis LINKDEF LinkDef.h)

set(GENERAL_LIBS ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree  ROOT::Hist ROOT::Gpad 
boost_program_options
)

# Add the main executable
add_executable(T2K T2K.cpp)

# Link all the object files to the main executable
target_link_libraries(T2K PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)

add_executable(MINERvA_pi0 MINERvA_pi0.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_pi0 PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)

add_executable(MINERvA_0pi MINERvA_0pi.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_0pi PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)


add_executable(rawtheta rawtheta.cpp)

# Link all the object files to the main executable
target_link_libraries(rawtheta PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)

add_executable(MINERvA_all MINERvA_all.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_all PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)


add_executable(MINERvA_all_gibuu MINERvA_all_gibuu.cpp)

# Link all the object files to the main executable
target_link_libraries(MINERvA_all_gibuu PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)


add_executable(gibuu2root gibuu2root.cpp)

# Link all the object files to the main executable
target_link_libraries(gibuu2root PRIVATE 
  ${GENERAL_LIBS}
)


add_executable(rawq2 rawq2.cpp)

# Link all the object files to the main executable
target_link_libraries(rawq2 PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)

add_executable(rawq2_gibuu rawq2_gibuu.cpp)

# Link all the object files to the main executable
target_link_libraries(rawq2_gibuu PRIVATE 
  ${GENERAL_LIBS}
  EVENT1
  analysis
)


add_executable(minerva_all_gibuu_plot minerva_all_gibuu_plot.cpp)

# Link all the object files to the main executable
target_link_libraries(minerva_all_gibuu_plot PRIVATE 
  ${GENERAL_LIBS}
  analysis
)


add_executable(t2k_stk_gibuu t2k_stk_gibuu.cpp)

# Link all the object files to the main executable
target_link_libraries(t2k_stk_gibuu PRIVATE 
  ${GENERAL_LIBS}
  analysis
)


add_executable(GiBUU2pi GiBUU2pi.cpp)

# Link all the object files to the main executable
target_link_libraries(GiBUU2pi PRIVATE 
  ${GENERAL_LIBS}
  analysis
)


add_executable(microboone_pi0_gibuu microboone_pi0_gibuu.cpp)

# Link all the object files to the main executable
target_link_libraries(microboone_pi0_gibuu PRIVATE 
  ${GENERAL_LIBS}
  analysis
)
