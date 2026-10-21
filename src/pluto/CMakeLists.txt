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
    ${VTK_INCLUDE_DIRS}

    ${PROJECT_SOURCE_DIR}/src/common
)

set(SRC_PLUTOVTK
    pluto_vtk_extract_iolib.cpp
    pluto_vtk_extract_iolib.h
)

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
endif()

if(WITH_OPENVDB OR WITH_NANOVDB)
    list(APPEND SRC_PLUTOVTK
        pluto_vtk_convert_vdb.cpp
        pluto_vtk_convert_vdb.h
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

if(WITH_VTK)
    add_definitions(-DWITH_VTK)
endif()

include_directories(${INC})

########################################################
add_library(pluto_vtk_converter ${SRC_PLUTOVTK})
target_link_libraries(
   pluto_vtk_converter 
   $<IF:$<CONFIG:Debug>,${OPENVDB_LIBRARIES_DEBUG},${OPENVDB_LIBRARIES}>
   $<IF:$<CONFIG:Debug>,${TBB_LIBRARIES_DEBUG},${TBB_LIBRARIES}>
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
   ${VTK_LIBRARIES}
)

target_link_libraries(pluto_vtk_converter space_common)
    
install (TARGETS pluto_vtk_converter 
    EXPORT space_converter_targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)
