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

    ${PROJECT_SOURCE_DIR}/src/common
)

set(SRC
    csv_extract_iolib.cpp
    csv_extract_iolib.h
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_OPENVDB OR WITH_NANOVDB)
    list(APPEND SRC
        csv_convert_vdb.cpp
        csv_convert_vdb.h
    )
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

include_directories(${INC})

########################################################
add_library(csv_converter ${SRC})
target_link_libraries(
   csv_converter 
   debug ${OPENVDB_LIBRARIES_DEBUG}
   optimized ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
)

target_link_libraries(csv_converter space_common)
    
install (TARGETS csv_converter DESTINATION lib)
########################################################
