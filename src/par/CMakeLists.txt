############################################################################
##
## BSD 3-Clause License
##
## Copyright (c) 2021, The Regents of the University of California
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
## * Redistributions of source code must retain the above copyright notice, this
##   list of conditions and the following disclaimer.
##
## * Redistributions in binary form must reproduce the above copyright notice,
##   this list of conditions and the following disclaimer in the documentation
##   and/or other materials provided with the distribution.
##
## * Neither the name of the copyright holder nor the names of its
##   contributors may be used to endorse or promote products derived from
##   this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
##
############################################################################

include("openroad")

option(LOAD_CPLEX "Load CPLEX" OFF)
if (LOAD_CPLEX)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_MODULE_PATH               "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  set(CMAKE_CXX_FLAGS                 "${CMAKE_CXX_FLAGS} -DIL_STD -m64 -Wall -fPIC")
  set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX   "-isystem ")
  set(LINKER_OPTIONS                  "-Wl,--export-dynamic")
  find_package(Cplex)
  include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})
endif (LOAD_CPLEX)

find_package(Threads REQUIRED)
find_package(ortools REQUIRED)
find_package(CapnProto REQUIRED)

capnp_generate_cpp(CAPNP_SRCS CAPNP_HDRS
  delay_graph.capnp
)

add_executable(text_graph_loader_test
  src/tests/text_loader_test.cpp
  src/delay_graph.cpp
)

add_executable(initial_loop_test
  src/tests/test_initial_loops.cpp
  src/find_initial_loops.cpp
  src/delay_graph.cpp
)

add_library(par_lib
  ${CAPNP_SRCS}
  src/PartitionMgr.cpp
  src/TritonPart.cpp
  src/Hypergraph.cpp
  src/Utilities.cpp
  src/delay_graph.cpp
  src/Coarsener.cpp
  src/Multilevel.cpp
  src/Refiner.cpp
  src/Partitioner.cpp
  src/Evaluator.cpp
  src/GreedyRefine.cpp
  src/ILPRefine.cpp
  src/KWayFMRefine.cpp
  src/KWayPMRefine.cpp
  src/PriorityQueue.cpp
  src/delay_graph.cpp
  src/delay_graph_to_hypergraph.cpp
  src/load_delay_graph_capnp.cpp
)

target_include_directories(par_lib
  PUBLIC
    include
    ${CMAKE_CURRENT_BINARY_DIR}

)

target_link_libraries(par_lib
  PRIVATE
    odb
    OpenSTA
    utl_lib
    dbSta_lib
    ortools::ortools
    CapnProto::capnp
)

if (LOAD_CPLEX)
  target_compile_definitions(par_lib
    PRIVATE
      LOAD_CPLEX
  )

  target_link_libraries(par_lib
    PRIVATE
      ${CPLEX_LIBRARIES}
      ${CMAKE_THREAD_LIBS_INIT}
      dl
  )
endif (LOAD_CPLEX)


swig_lib(NAME      par
         NAMESPACE par
         I_FILE    src/partitionmgr.i
         SCRIPTS   src/partitionmgr.tcl
)

target_sources(par
  PRIVATE
    src/MakePartitionMgr.cpp
)

messages(
  TARGET par_lib
)

target_include_directories(par
  PUBLIC
    include
)

target_link_libraries(par
  PRIVATE
    odb
    OpenSTA
    par_lib
)

target_link_libraries(initial_loop_test
  PRIVATE
    ortools::ortools
)

if (LOAD_CPLEX)
  target_compile_definitions(par
    PRIVATE
      LOAD_CPLEX
  )

  target_link_libraries(par
    PRIVATE
      ${CPLEX_LIBRARIES}
      ${CMAKE_THREAD_LIBS_INIT}
      dl
  )
endif (LOAD_CPLEX)


if (Python3_FOUND AND BUILD_PYTHON)
  swig_lib(NAME          par_py
           NAMESPACE     par
           LANGUAGE      python
           I_FILE        src/partitionmgr-py.i
           SWIG_INCLUDES ${PROJECT_SOURCE_DIR}/include/par
           SCRIPTS       ${CMAKE_CURRENT_BINARY_DIR}/par_py.py
  )

  target_include_directories(par_py
    PUBLIC
      include
  )

  target_link_libraries(par_py
    PUBLIC
      par_lib
  )

endif()
