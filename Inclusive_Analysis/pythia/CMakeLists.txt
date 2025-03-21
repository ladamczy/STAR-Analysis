
cmake_minimum_required (VERSION 3.0)
project (pythia)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#PYTHIA stuff
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
find_package(Pythia8 REQUIRED)

message(${PYTHIA8_INCLUDE_DIR})

#folder with shared code
set(SHARED_FOLDER ../../star-upc)

#sources for shared library LIB
file(GLOB_RECURSE SRCS "${SHARED_FOLDER}/src/*.cxx")
list(TRANSFORM SRCS REPLACE ${CMAKE_CURRENT_SOURCE_DIR}/ "")

#name of shared library LIB
set (LIB star-upc)

#binary executables
file(GLOB_RECURSE BIN "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx")

#sources for ROOT dictionary
set(SRCS_DICT ${SRCS})

#LinkDef file for ROOT dictionary
set (LINKDEF "${LIB}LinkDef.h")
set (LINKDEF_PATH "${SHARED_FOLDER}/include/${LINKDEF}")

#general compiler options
set (CMAKE_CXX_COMPILER /usr/bin/g++)
#option -Wl,--no-as-needed is necessary only with root5
set (CMAKE_CXX_FLAGS "-Wall -g -Wl,--no-as-needed")
include_directories ("${SHARED_FOLDER}/include")

#ROOT section

#compiler flags
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_FLAGS_CMD)
string(REPLACE "\n" "" ROOT_FLAGS "${ROOT_FLAGS_CMD}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_FLAGS}")

#libraries
set(ROOT_DEPS Core EG Hist Physics RIO Tree MathCore Gpad Thread Imt MultiProc Graf3d Graf Net Matrix TreePlayer)

#include dir
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOT_LIBDIR_CMD)
string(REPLACE "\n" "" ROOT_LIBDIR "${ROOT_LIBDIR_CMD}")
link_directories (${ROOT_LIBDIR})

#generate ROOT dictionary by CINT/CLING
set (CLING rootcint)
set (DICT G__${LIB}.cxx)

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS_DICT}")
string(REPLACE "${SHARED_FOLDER}/src/" "" HDRS_PLAIN "${HDRS}")
string(REPLACE "src/" "include/" HDRS_INC "${HDRS}")

# directory with headers for CLING
get_directory_property(incdirs INCLUDE_DIRECTORIES)

# command to invoke CLING, options -c -p are used by root5
# and deprecated in root6
message(${CMAKE_CURRENT_BINARY_DIR})

add_custom_command (
  OUTPUT ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${DICT}
  COMMAND ${CLING}
  ARGS -f ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${DICT} -c -p -I${incdirs}
  ${HDRS_PLAIN} ${LINKDEF}
  DEPENDS ${HDRS_INC} ${LINKDEF_PATH}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

#end of ROOT section

#create library LIB
add_library (${LIB} SHARED ${SRCS} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${DICT})
target_link_libraries(${LIB} ${ROOT_DEPS} )

#add custom (mine) headers
include_directories(include)
include_directories(../common_assets)

#add PYTHIA stuff
include_directories(${PYTHIA8_INCLUDE_DIR})

#create executables in BIN
foreach(IBIN ${BIN})
  #IBIN is full path, TEMP_BIN only the filename
  get_filename_component(TEMP_BIN ${IBIN} NAME_WLE)
  add_executable(${TEMP_BIN} ${IBIN})
  target_link_libraries(${TEMP_BIN} ${LIB} ${ROOT_DEPS} ${PYTHIA8_LIBRARIES})
endforeach(IBIN ${BIN})