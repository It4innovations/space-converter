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
    ${PROJECT_SOURCE_DIR}/src/changa/utility
    ${PROJECT_SOURCE_DIR}/src/changa/utility/structures    
    ${PROJECT_SOURCE_DIR}/src/changa/utility/xdr
)

set(SRC_TIPSY
    changa_tipsy_extract_iolib.cpp
    changa_tipsy_extract_iolib.h

    utility/structures/TipsyReader.cpp
    utility/structures/TipsyFile.cpp    
)

set(SRC_NCHILADA
    changa_nchilada_extract_iolib.cpp
    changa_nchilada_extract_iolib.h

    utility/xdr/htonl.c
    utility/xdr/types.h
    utility/xdr/xdr.c
    utility/xdr/xdr.h
    utility/xdr/xdr_array.c
    utility/xdr/xdr_float.c
    utility/xdr/xdr_intXX_t.c
    utility/xdr/xdr_mem.c
    utility/xdr/xdr_sizeof.c
    utility/xdr/xdr_stdio.c    
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_OPENVDB OR WITH_NANOVDB)
    list(APPEND SRC_TIPSY
        changa_tipsy_convert_vdb.cpp
        changa_tipsy_convert_vdb.h
    )

    list(APPEND SRC_NCHILADA
        changa_nchilada_convert_vdb.cpp
        changa_nchilada_convert_vdb.h
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

add_definitions(-DWITH_IO_CHANGA)

include_directories(${INC})

########################TIPSY################################
add_library(changa_tipsy_converter ${SRC_TIPSY})
target_link_libraries(
   changa_tipsy_converter 
   debug ${OPENVDB_LIBRARIES_DEBUG}
   optimized ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
)

target_link_libraries(changa_tipsy_converter space_common)

install (TARGETS changa_tipsy_converter DESTINATION lib)
###########################NCHILADA#############################
add_library(changa_nchilada_converter ${SRC_NCHILADA})
target_link_libraries(
   changa_nchilada_converter 
   debug ${OPENVDB_LIBRARIES_DEBUG}
   optimized ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
)

target_link_libraries(changa_nchilada_converter space_common)
    
install (TARGETS changa_nchilada_converter DESTINATION lib)