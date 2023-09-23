project (ExclusiveAnalysis)
cmake_minimum_required (VERSION 2.8.12)

#sources for LIB_ExclusiveAnalysis
set (SRCS
   src/MonteCarloEfficiency.cxx
#  src/MatchFillPosition.cxx		
#  src/ReadFillPositionFile.cxx
  )

#name of LIB_ExclusiveAnalysis
set (LIB_ExclusiveAnalysis ExclusiveAnalysis)

#binary executables
set (BIN MonteCarloEfficiency)

#UPC picoDsts library
set (UPC_LIB ../../star-upc/build/libstar-upc.so)
set (UPC_INC ../star-upc/include)

#general compiler options
set (CMAKE_CXX_COMPILER /usr/bin/g++)
#option -Wl,--no-as-needed is necessary only with root5
#set (CMAKE_CXX_FLAGS "-Wall -g -Wl,--no-as-needed")
set (CMAKE_CXX_FLAGS "-Wall -g")
include_directories (include ${UPC_INC})

#ROOT section

#compiler flags
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_FLAGS_CMD)
string(REPLACE "\n" "" ROOT_FLAGS "${ROOT_FLAGS_CMD}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_FLAGS}")

#libraries
set(ROOT_DEPS Core EG Hist Physics RIO Tree MathCore)

#include dir
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOT_LIBDIR_CMD)
string(REPLACE "\n" "" ROOT_LIBDIR "${ROOT_LIBDIR_CMD}")
link_directories (${ROOT_LIBDIR})

#end of ROOT section

#create library LIC_ExclusiveAnalysis
add_library (${LIB_ExclusiveAnalysis} SHARED ${SRCS})
target_link_libraries(${LIB_ExclusiveAnalysis} ${ROOT_DEPS} )
#create executables in BIN
foreach(IBIN ${BIN})
  add_executable(${IBIN} src/${IBIN}.cxx)
  target_link_libraries(${IBIN} ${LIB_ExclusiveAnalysis} ${UPC_LIB} ${ROOT_DEPS} )
endforeach(IBIN ${BIN})


