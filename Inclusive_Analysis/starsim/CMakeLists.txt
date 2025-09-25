
cmake_minimum_required (VERSION 3.0)
project (starsim)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#because root_generate_dictionary uses
#some additional directories
include(GNUInstallDirs)

#outside libraries
find_package(ROOT REQUIRED COMPONENTS EG)
# set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
# find_package(Pythia8)

#folder with shared code
set(SHARED_FOLDER ../../star-upc-new)

#headers
include_directories(include)
include_directories(../common_assets)
include_directories(../star-upc/include)
include_directories("${SHARED_FOLDER}/include")
include_directories(${ROOT_INCLUDE_DIRS})
# include_directories(${PYTHIA8_INCLUDE_DIR})

#sources for shared library LIB
#both from star-upc and starsim
file(GLOB_RECURSE SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${SHARED_FOLDER}/src/*.cxx" "${CMAKE_CURRENT_SOURCE_DIR}/srcs/*.cxx")

#name of shared library LIB
set (LIB starsim)

#LinkDef file for ROOT dictionary
set (LINKDEF "${LIB}LinkDef.h")

#generate ROOT dictionary
set (DICT G__${LIB})
string(REPLACE "src/" "include/" INCLUDES_NAMES_ONLY "${SRCS}")
string(REPLACE "srcs/" "" INCLUDES_NAMES_ONLY "${INCLUDES_NAMES_ONLY}")
string(REPLACE ".cxx" ".h" INCLUDES_NAMES_ONLY "${INCLUDES_NAMES_ONLY}")
root_generate_dictionary(${DICT} ${INCLUDES_NAMES_ONLY} LINKDEF "include/${LINKDEF}")

#create library LIB
add_library (${LIB} SHARED ${SRCS} ${DICT})
target_link_libraries(${LIB} ${ROOT_LIBRARIES})
# target_link_libraries(${LIB} ${ROOT_LIBRARIES} ${PYTHIA8_LIBRARIES})

#create executables in BIN
file(GLOB_RECURSE BIN "${CMAKE_CURRENT_SOURCE_DIR}/execs/*.cxx")
foreach(IBIN ${BIN})
  #IBIN is full path, TEMP_BIN only the filename
  get_filename_component(TEMP_BIN ${IBIN} NAME_WE)
  add_executable(${TEMP_BIN} ${IBIN})
  target_link_libraries(${TEMP_BIN} ${LIB})
endforeach(IBIN ${BIN})