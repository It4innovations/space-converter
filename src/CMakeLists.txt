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

if(WITH_NANOVDB)

    # Check if NANOVDB_INCLUDE_DIRS is not set or is an empty string
    if(NOT NANOVDB_INCLUDE_DIRS OR NANOVDB_INCLUDE_DIRS STREQUAL "")
        # If not set or empty, then set it to the specified value
        set(NANOVDB_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/submodules/openvdb/nanovdb")
    endif()

endif()

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
    main.cpp

    data_processing.h
    data_processing.cpp

    data_communication.h
    data_communication.cpp

    args_processing.h
    args_processing.cpp
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
endif()

if(WITH_MULTIRES)
    add_definitions(-DWITH_MULTIRES)
endif()

if(WITH_NO_DATA_TEMP)
    add_definitions(-DWITH_NO_DATA_TEMP)
endif()

if(WITH_CUDAKDTREE)
    add_definitions(-DWITH_CUDAKDTREE)
    #find_package(CUDA REQUIRED)

    list(APPEND INC
        ${CUDA_INCLUDE_DIRS}
    )

    # if (WITH_CUDAKDTREE_CPU)
    #     add_definitions(-DWITH_CUDAKDTREE_CPU)
    # endif()

endif()

include_directories(${INC})

########################################################
add_subdirectory(common)
add_subdirectory(utility)
add_subdirectory(gadget)
add_subdirectory(changa)
add_subdirectory(csv)
add_subdirectory(haccbin)

if(WITH_GENERICIO)
    add_subdirectory(genericio)
endif()

if(WITH_HDF5)
    add_subdirectory(hdf5)
endif()
########################################################
add_executable(space_converter ${SRC})
target_link_libraries(
   space_converter 
   ${OPENVDB_LIBRARIES}
   ${BOOST_LIBRARIES}
   ${TBB_LIBRARIES}
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
   ${MPI_CUSTOM_LIBRARIES}
)

target_link_libraries(space_converter 
    space_common
    gadget_converter
    gadget_simple_converter
    changa_tipsy_converter
    changa_nchilada_converter
    csv_converter    
    haccbin_converter    
)

if(NOT WIN32 AND NOT APPLE)
    target_link_libraries(space_converter stdc++fs)
endif()

if(WITH_GENERICIO)
    target_link_libraries(space_converter genericio_converter)
endif()

if(WITH_HDF5)
    target_link_libraries(space_converter hdf5_converter)
endif()

# Add this conditional block to link CUDA libraries
# if(WITH_CUDAKDTREE)
#     target_link_libraries(space_converter 
#         cudakdtree_tool
#         ${CUDA_LIBRARIES}
#         ${CUDA_CUDART_LIBRARY}
#     )
# endif()
if(WITH_CUDAKDTREE)
    set_target_properties(space_converter PROPERTIES
         CUDA_SEPARABLE_COMPILATION ON
         CUDA_STANDARD 14
         POSITION_INDEPENDENT_CODE ON
         #CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES}
    )
endif()    


install (TARGETS space_converter DESTINATION bin)

