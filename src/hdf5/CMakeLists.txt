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
    ${HDF5_INCLUDE_DIRS}

    ${PROJECT_SOURCE_DIR}/src/common
)

set(SRC
    hdf5_extract_iolib.cpp
    hdf5_extract_iolib.h
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_OPENVDB OR WITH_NANOVDB)
    list(APPEND SRC
        hdf5_convert_vdb.cpp
        hdf5_convert_vdb.h
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

find_package(HDF5 REQUIRED COMPONENTS C CXX)

include_directories(${INC})

########################################################
add_library(hdf5_converter ${SRC})
target_link_libraries(
   hdf5_converter
   ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
   ${HDF5_C_LIBRARIES}
   ${HDF5_CXX_LIBRARIES}

   space_common
)

#target_compile_definitions(hdf5_converter PUBLIC -DHAVE_HDF5 -DH5_BUILT_AS_DYNAMIC_LIB)
target_compile_definitions(hdf5_converter PRIVATE ${HDF5_DEFINITIONS})
    
install (TARGETS hdf5_converter DESTINATION lib)
########################################################
