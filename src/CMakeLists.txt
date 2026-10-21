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
    ${VTK_INCLUDE_DIRS}
    ${TBB_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIRS}
    ${MPI_CXX_HEADER_DIR}
    ${PROJECT_SOURCE_DIR}/src/common
    ${MERIC_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/submodules/braas-hpc-renderengine/src
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

if(WITH_GPU_CUDA)
    add_definitions(-DWITH_GPU_CUDA)
endif()

if(WITH_GPU_HIP)
    add_definitions(-DWITH_GPU_HIP)
    # Legacy guard macro: the GPU code paths in the sources are guarded by
    # WITH_GPU_CUDA; gpu_device_compat.h remaps the CUDA runtime API to HIP.
    add_definitions(-DWITH_GPU_CUDA)
endif()

if(WITH_CUDAKDTREE)
    add_definitions(-DWITH_CUDAKDTREE)
    if(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES)
        list(APPEND INC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    endif()
endif()

if(WITH_CUDA_AWARE_MPI)
    add_definitions(-DWITH_CUDA_AWARE_MPI)
endif()

if(WITH_CUDA_MALLOCMANAGED)
    add_definitions(-DWITH_CUDA_MALLOCMANAGED)
endif()

if(WITH_HIP_AWARE_MPI)
    add_definitions(-DWITH_HIP_AWARE_MPI)
    add_definitions(-DWITH_CUDA_AWARE_MPI)
endif()

if(WITH_HIP_MALLOCMANAGED)
    add_definitions(-DWITH_HIP_MALLOCMANAGED)
    add_definitions(-DWITH_CUDA_MALLOCMANAGED)
endif()

if(WITH_OPENMP)
    add_definitions(-DWITH_OPENMP)
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

if(WITH_VTK)
    add_definitions(-DWITH_VTK)
endif()

if(WITH_TBB)
    add_definitions(-DWITH_TBB)
endif()

if(WITH_HACC)
    add_definitions(-DWITH_HACC)
endif()

if(WITH_IPIC3D)
    add_definitions(-DWITH_IPIC3D)
endif()

if(WITH_PLUTO)
    add_definitions(-DWITH_PLUTO)
endif()

if(WITH_NANOFLANN)
    add_definitions(-DWITH_NANOFLANN)
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

if(WITH_MERIC)
    add_definitions(-DWITH_MERIC)
endif()

include_directories(${INC})

########################################################
add_subdirectory(common)
add_subdirectory(utility)
add_subdirectory(gadget)
add_subdirectory(changa)

if(WITH_HACC)
    add_subdirectory(hacc)
endif()

if(WITH_IPIC3D)
    add_subdirectory(ipic3d)
endif()

if(WITH_PLUTO)
    add_subdirectory(pluto)
endif()

########################################################
add_executable(space_converter ${SRC})
target_link_libraries(
   space_converter 
   $<IF:$<CONFIG:Debug>,${OPENVDB_LIBRARIES_DEBUG},${OPENVDB_LIBRARIES}>
   $<IF:$<CONFIG:Debug>,${TBB_LIBRARIES_DEBUG},${TBB_LIBRARIES}>
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
   ${MPI_CUSTOM_LIBRARIES}
   ${MERIC_LIBRARIES}
   ${OpenMP_omp_LIBRARY}
)

target_link_libraries(space_converter 
    space_common
    gadget_converter
    gadget_simple_converter
    changa_tipsy_converter
    changa_nchilada_converter        
)

if(WITH_VTK)
    target_link_libraries(space_converter ${VTK_LIBRARIES})
endif()

if(NOT WIN32 AND NOT APPLE)
    target_link_libraries(space_converter stdc++fs)
endif()

if(WITH_GPU_HIP)
    target_link_libraries(space_converter hip::host)
    # Link with the C++ compiler (Cray wrapper): it adds libgomp/MPI/GTL for the
    # g++-compiled objects; the non-RDC HIP objects only need libamdhip64.
    set_target_properties(space_converter PROPERTIES LINKER_LANGUAGE CXX)
endif()

if(WITH_IPIC3D)
    target_link_libraries(space_converter ipic3d_hdf5_converter)
endif()

if(WITH_HACC)
    target_link_libraries(space_converter hacc_genericio_converter haccbin_converter)
endif()

if(WITH_PLUTO)
    target_link_libraries(space_converter pluto_vtk_converter)
endif()

if(WITH_CUDAKDTREE)
    set_target_properties(space_converter PROPERTIES
         CUDA_SEPARABLE_COMPILATION ON
         CUDA_STANDARD 14
         POSITION_INDEPENDENT_CODE ON
         #CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES}
    )
endif()    


install (TARGETS space_converter 
    EXPORT space_converter_targets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

