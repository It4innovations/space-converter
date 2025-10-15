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
    utility/OpenGadget3/IO
)

set(SRC
    gadget_extract_iolib.cpp
    gadget_extract_iolib.h
)

set(SRC_SIMPLE
    gadget_simple_extract_iolib.cpp
    gadget_simple_extract_iolib.h 
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_OPENVDB OR WITH_NANOVDB)
    list(APPEND SRC
        gadget_convert_vdb.cpp
        gadget_convert_vdb.h
    )

    list(APPEND SRC_SIMPLE
        gadget_simple_convert_vdb.cpp
        gadget_simple_convert_vdb.h
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

if(WITH_HDF5)
    find_package(HDF5 REQUIRED)
endif()

if(GADGET_MAX_HSML)
    add_definitions(-DGADGET_MAX_HSML)
endif()

if(GADGET_READ_ID)
    add_definitions(-DGADGET_READ_ID)
endif()

if(GADGET_READ_ID64)
    add_definitions(-DGADGET_READ_ID64)
endif()

if(GADGET3_IO_LIB)
    add_definitions(-DGADGET3_IO_LIB)
endif()

if(GADGET3_IO_LIB_MPI)
    add_definitions(-DGADGET3_IO_LIB_MPI)
endif()

include_directories(${INC})

########################################################
add_subdirectory(utility)
########################################################
add_library(gadget_converter ${SRC})
target_link_libraries(
   gadget_converter 
   debug ${OPENVDB_LIBRARIES_DEBUG}
   optimized ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
)

target_link_libraries(gadget_converter space_common gadget3_io_lib_static)
    
install (TARGETS gadget_converter DESTINATION lib)
########################################################
add_library(gadget_simple_converter ${SRC_SIMPLE})
target_link_libraries(
   gadget_simple_converter
   debug ${OPENVDB_LIBRARIES_DEBUG}
   optimized ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
)

target_link_libraries(gadget_simple_converter space_common gadget3_io_lib_static)

install (TARGETS gadget_simple_converter DESTINATION lib)
########################################################