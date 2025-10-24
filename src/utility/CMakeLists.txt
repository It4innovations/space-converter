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
    # ${BOOST_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIRS}
    ${MPI_CXX_HEADER_DIR}
    ${MPI_CXX_COMPILER_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/src/common
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

include_directories(${INC})
########################################################
if(WITH_NANOVDB AND WITH_OPENVDB)
    set(SRC_NANOVDB
        nano2vdb.cpp
        stb_image_write.h
    )

    add_executable(nano2vdb ${SRC_NANOVDB})
    target_link_libraries(
        nano2vdb 
        ${OPENVDB_LIBRARIES}
        ${BOOST_LIBRARIES}
        ${TBB_LIBRARIES}
        ${MPI_LIBRARIES}
        ${MPI_msmpi_LIBRARY}
    )
        
    install (TARGETS nano2vdb DESTINATION bin)
endif()
########################################################
if(WITH_NANOVDB AND WITH_OPENVDB)
    set(SRC_NANOVDB
        vdb2nano.cpp
        stb_image_write.h
    )

    add_executable(vdb2nano ${SRC_NANOVDB})
    target_link_libraries(
        vdb2nano 
        ${OPENVDB_LIBRARIES}
        ${BOOST_LIBRARIES}
        ${TBB_LIBRARIES}
        ${MPI_LIBRARIES}
        ${MPI_msmpi_LIBRARY}
    )
        
    install (TARGETS vdb2nano DESTINATION bin)
endif()
########################################################
if(WITH_OPENVDB)
    set(SRC_VDB_MERGE
        vdb_merge.cpp
    )

    add_executable(vdb_merge ${SRC_VDB_MERGE})
    target_link_libraries(
        vdb_merge 
        ${OPENVDB_LIBRARIES}
        ${BOOST_LIBRARIES}
        ${TBB_LIBRARIES}
        ${MPI_LIBRARIES}
        ${MPI_msmpi_LIBRARY}
    )
        
    install (TARGETS vdb_merge DESTINATION bin)
endif()
########################################################
if(WITH_CUDAKDTREE)
    set(SRC_CUDAKDTREE
        cudakdtree_tool.cu
        cudakdtree_tool.h
    )    

    add_library(cudakdtree_tool STATIC ${SRC_CUDAKDTREE})
    target_include_directories(cudakdtree_tool PUBLIC ${PROJECT_SOURCE_DIR}/submodules/cudaKDTree ${INC})

    target_link_libraries(
        cudakdtree_tool
        MPI::MPI_CXX
    )
        
    install (TARGETS cudakdtree_tool DESTINATION lib)
endif()
########################################################
if(WITH_NANOFLANN)
    set(SRC_NANOFLANN
        nanoflann_tool.cpp
        nanoflann_tool.h
    )    

    add_library(nanoflann_tool STATIC ${SRC_NANOFLANN})
    target_include_directories(nanoflann_tool PUBLIC ${PROJECT_SOURCE_DIR}/submodules/nanoflann/include/)

    target_link_libraries(
        nanoflann_tool
        MPI::MPI_CXX
    )
        
    install (TARGETS nanoflann_tool DESTINATION lib)
endif()
########################################################
if(WITH_OPENVDB)
    set(SRC_VDB2HISTO
        vdb2histo.cpp
    )

    add_executable(vdb2histo ${SRC_VDB2HISTO})
    target_link_libraries(
        vdb2histo 
        ${OPENVDB_LIBRARIES}
        ${BOOST_LIBRARIES}
        ${TBB_LIBRARIES}
        ${MPI_LIBRARIES}
        ${MPI_msmpi_LIBRARY}
    )
        
    install (TARGETS vdb2histo DESTINATION bin)
endif()
########################################################
if(WITH_OPENVDB)
    set(SRC_VDB2HISTO
        vdb2png.cpp
    )

    add_executable(vdb2png ${SRC_VDB2HISTO})
    target_link_libraries(
        vdb2png 
        ${OPENVDB_LIBRARIES}
        ${BOOST_LIBRARIES}
        ${TBB_LIBRARIES}
        ${MPI_LIBRARIES}
        ${MPI_msmpi_LIBRARY}
    )
        
    install (TARGETS vdb2png DESTINATION bin)
endif()
########################################################
