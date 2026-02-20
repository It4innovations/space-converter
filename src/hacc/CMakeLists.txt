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

    ${PROJECT_SOURCE_DIR}/src/common
    utility/GenericIO
)

set(SRC_GENERICIO
    hacc_genericio_extract_iolib.cpp
    hacc_genericio_extract_iolib.h
)

set(SRC_HACCBIN
    haccbin_extract_iolib.cpp
    haccbin_extract_iolib.h
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_OPENVDB OR WITH_NANOVDB)
    list(APPEND SRC_GENERICIO
        hacc_genericio_convert_vdb.cpp
        hacc_genericio_convert_vdb.h
    )

    list(APPEND SRC_HACCBIN
        haccbin_convert_vdb.cpp
        haccbin_convert_vdb.h
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

include_directories(${INC})

########################################################
add_subdirectory(utility)
########################################################
add_library(hacc_genericio_converter ${SRC_GENERICIO})
target_link_libraries(
   hacc_genericio_converter 
   $<IF:$<CONFIG:Debug>,${OPENVDB_LIBRARIES_DEBUG},${OPENVDB_LIBRARIES}>
   $<IF:$<CONFIG:Debug>,${TBB_LIBRARIES_DEBUG},${TBB_LIBRARIES}>
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
)

target_link_libraries(hacc_genericio_converter space_common genericio_lib_static sz_lib)
    
install (TARGETS hacc_genericio_converter DESTINATION lib)
########################################################
add_library(haccbin_converter ${SRC_HACCBIN})
target_link_libraries(
   haccbin_converter 
   $<IF:$<CONFIG:Debug>,${OPENVDB_LIBRARIES_DEBUG},${OPENVDB_LIBRARIES}>
   $<IF:$<CONFIG:Debug>,${TBB_LIBRARIES_DEBUG},${TBB_LIBRARIES}>
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
)

target_link_libraries(haccbin_converter space_common)
    
install (TARGETS haccbin_converter DESTINATION lib)
