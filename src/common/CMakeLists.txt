#####################################################################################################################
# Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
#
# This program is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#####################################################################################################################

set(INC
    .  
    ${OPENVDB_INCLUDE_DIRS}
    ${NANOVDB_INCLUDE_DIRS}
    ${BOOST_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIRS}
    ${MPI_CXX_HEADER_DIR}
    ${EMBREE_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/src/common
    ${PROJECT_SOURCE_DIR}/src/utility
    ${PROJECT_SOURCE_DIR}/submodules/nanoflann/include
)

set(SRC_COMMON
    connection_tcp.cpp
    connection_tcp.h
    convert_vdb.cpp
    convert_vdb.h
    convert_common.cpp
    convert_common.h
    dense_common.cpp
    dense_common.h

    data_common.h

    ${PROJECT_SOURCE_DIR}/src/utility/dense_utility.h
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_NANOVDB)
    add_definitions(-DWITH_NANOVDB)

    list(APPEND INC
        ${NANOVDB_INCLUDE_DIRS}
    )    

    add_definitions(-DOPENVDB_VERSION=${OPENVDB_VERSION})
endif()

if(WITH_OPENVDB)
    add_definitions(-DWITH_OPENVDB)
    add_definitions(-DOPENVDB_VERSION=${OPENVDB_VERSION})
endif()

if(WITH_MPI)
    add_definitions(-DWITH_MPI)
endif()

if(WITH_TBB)
    add_definitions(-DWITH_TBB)
endif()

if(WITH_EMBREE)
    add_definitions(-DWITH_EMBREE)
endif()

if(WITH_GENERICIO)
    add_definitions(-DWITH_GENERICIO)
endif()

if(WITH_NANOFLANN)
    add_definitions(-DWITH_NANOFLANN)     
    list(APPEND SRC_COMMON
        ${PROJECT_SOURCE_DIR}/submodules/nanoflann/include/nanoflann.hpp
    )      
endif()

if(WITH_CUDAKDTREE)
    add_definitions(-DWITH_CUDAKDTREE)

    # if (WITH_CUDAKDTREE_CPU)
    #     add_definitions(-DWITH_CUDAKDTREE_CPU)
    # endif()    
endif()

if(WITH_MULTIRES)
    add_definitions(-DWITH_MULTIRES)
endif()

if(WITH_NO_DATA_TEMP)
    add_definitions(-DWITH_NO_DATA_TEMP)
endif()

include_directories(${INC})
########################################################
add_library(space_common ${SRC_COMMON})
target_link_libraries(
   space_common 
   ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${EMBREE_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
   ${OpenMP_omp_LIBRARY}
)

if(WITH_CUDAKDTREE)
   target_link_libraries(space_common cudakdtree_tool)
endif()

if(WITH_NANOFLANN)
   target_link_libraries(space_common nanoflann_tool)
endif()

install (TARGETS space_common DESTINATION lib)
