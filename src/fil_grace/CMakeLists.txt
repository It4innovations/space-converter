#####################################################################################################################
# Copyright(C) 2023-2026 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
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
    ${TBB_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIRS}
    ${MPI_CXX_HEADER_DIR}
    ${HDF5_INCLUDE_DIRS}

    ${PROJECT_SOURCE_DIR}/src/common
)

set(SRC_FIL_GRACE
    fil_grace_extract_iolib.cpp
    fil_grace_extract_iolib.h
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_OPENVDB OR WITH_NANOVDB)
    list(APPEND SRC_FIL_GRACE
        fil_grace_convert_vdb.cpp
        fil_grace_convert_vdb.h
    )
endif()

if(WITH_NANOVDB)
    add_definitions(-DWITH_NANOVDB)

    list(APPEND INC
        ${NANOVDB_INCLUDE_DIRS}
    )    
endif()

if(WITH_OPENVDB)
    add_definitions(-DWITH_OPENVDB)
endif()

if(WITH_MPI)
    add_definitions(-DWITH_MPI)
endif()

if(WITH_TBB)
    add_definitions(-DWITH_TBB)
endif()

if(WITH_HDF5)
    add_definitions(-DWITH_HDF5)
endif()


include_directories(${INC})

########################################################
add_library(fil_grace_converter ${SRC_FIL_GRACE})
target_link_libraries(
   fil_grace_converter 
   $<IF:$<CONFIG:Debug>,${OPENVDB_LIBRARIES_DEBUG},${OPENVDB_LIBRARIES}>
   $<IF:$<CONFIG:Debug>,${TBB_LIBRARIES_DEBUG},${TBB_LIBRARIES}>
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
   ${HDF5_LIBRARIES}
)

target_link_libraries(fil_grace_converter space_common)
    
install (TARGETS fil_grace_converter 
    EXPORT space_converter_targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)
