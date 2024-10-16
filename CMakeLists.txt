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

# if (NOT NUWRO_DIR)
#   # message(FATAL_ERROR "Please set NUWRO environment variable to the path to the NuWro installation")
#   set(NUWRO_DIR /var/home/yan/neutrino/nuwro/nuwro.nuwro_21.09.2)
# endif()
if (NUWRO_DIR)
  add_library(EVENT1 SHARED IMPORTED)
  set_property(TARGET EVENT1 PROPERTY IMPORTED_LOCATION ${NUWRO_DIR}/bin/event1.so)
  include_directories(${NUWRO_DIR}/src ${NUWRO_DIR}/include)
endif()

# Include ROOT and NuWro headers
include_directories(${ROOT_INCLUDE_DIRS} include)


set(GLOB Dict_HEADERS include/tkigeneral.h
                      include/tkievent.h)
set(COMMON_TOOL src/dochi2.cxx 
                src/expdata.cxx 
                src/MINERvA_tki_cut.cxx 
                src/plottools.cxx 
                src/T2K_tki_cut.cxx
                src/tkievent.cxx
                src/tkigeneral.cxx)
set(GIBUU_TOOL src/EvtTracker2event.cxx)
set(NuWro_TOOL src/NuWro2event.cxx )

add_library(analysis SHARED ${COMMON_TOOL} ${GIBUU_TOOL})
if(NUWRO_DIR)
  target_sources(analysis PRIVATE ${NuWro_TOOL})
  target_link_libraries(analysis PRIVATE EVENT1)
endif()
target_link_libraries(analysis PRIVATE ROOT::Physics ROOT::ROOTDataFrame)

target_compile_definitions(analysis PRIVATE 
  PIZERO_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data/PIZEROTKI_MINERvA.root"
  ZEROPI_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data/PhysRevD.101.092001.root"
  T2K_STK_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data/datResults.root"
)

root_generate_dictionary(G__Dict ${Dict_HEADERS} MODULE analysis LINKDEF LinkDef.h)


set(GENERAL_LIBS ROOT::Physics ROOT::ROOTDataFrame ROOT::Tree  ROOT::Hist ROOT::Gpad 
    boost_program_options
)

if(NUWRO_DIR)
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


endif()

add_executable(gibuu2root gibuu2root.cpp)

# Link all the object files to the main executable
target_link_libraries(gibuu2root PRIVATE 
  ${GENERAL_LIBS}
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
