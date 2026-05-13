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
    ${EMBREE_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/src/common
    ${PROJECT_SOURCE_DIR}/src/utility
    ${PROJECT_SOURCE_DIR}/submodules/nanoflann/include
    ${PROJECT_SOURCE_DIR}/submodules/braas-hpc-renderengine/src
)

set(SRC_COMMON
    convert_vdb.cpp
    convert_vdb.h
    convert_common.cpp
    convert_common.h

    dense_common.cpp
    dense_common.h

    sparse_common.h
    sparse_common.cpp

    raw_common.h

    data_common.h

    data_cache.h
    data_cache.cpp

    convert_vdb_kernel.h
    convert_vdb_kernel.cpp

    ${PROJECT_SOURCE_DIR}/submodules/braas-hpc-renderengine/src/renderengine_tcp.cpp
    ${PROJECT_SOURCE_DIR}/submodules/braas-hpc-renderengine/src/renderengine_tcp.h

    ${PROJECT_SOURCE_DIR}/src/utility/dense_utility.h
)

if(WITH_GPU_CUDA)
    add_definitions(-DWITH_GPU_CUDA)
    list(APPEND SRC_COMMON
        data_cache.cu
        convert_vdb_kernel.cu
        sparse_common.cu
        dense_common.cu
    )
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

if(WITH_NANOFLANN)
    add_definitions(-DWITH_NANOFLANN)     
    list(APPEND SRC_COMMON
        ${PROJECT_SOURCE_DIR}/submodules/nanoflann/include/nanoflann.hpp
    )      
endif()

if(WITH_CUDAKDTREE)
    add_definitions(-DWITH_CUDAKDTREE)  
endif()

if(WITH_NO_DATA_TEMP)
    add_definitions(-DWITH_NO_DATA_TEMP)
endif()

include_directories(${INC})
########################################################
add_library(space_common ${SRC_COMMON})

if(WITH_GPU_CUDA)
    set_target_properties(space_common PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
    )
    
    # Add extended lambda support for NanoVDB CUDA features
    target_compile_options(space_common PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
    )
    
    if(WIN32)
        # Add MSVC standard conforming preprocessor flag for CUDA
        target_compile_options(space_common PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
        )
        # Suppress CUDA warnings about __host__ function calls from __host__ __device__ functions
        # and stack-size-cannot-be-statically-determined (1408) caused by deep SPH dispatch chains
        target_compile_options(space_common PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-diag-suppress=20011,1408>
        )
        # Add debug symbols for CUDA device code in Debug builds
        target_compile_options(space_common PRIVATE
            $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>:-g;-G>
        )
    endif()
endif()

target_link_libraries(
   space_common 
   $<IF:$<CONFIG:Debug>,${OPENVDB_LIBRARIES_DEBUG},${OPENVDB_LIBRARIES}>
   $<IF:$<CONFIG:Debug>,${TBB_LIBRARIES_DEBUG},${TBB_LIBRARIES}>
   ${MPI_LIBRARIES}
   ${MPI_msmpi_LIBRARY}
   ${OpenMP_omp_LIBRARY}
)

if(WITH_CUDAKDTREE)
   target_link_libraries(space_common cudakdtree_tool)
endif()

if(WITH_NANOFLANN)
   target_link_libraries(space_common nanoflann_tool)
endif()

install (TARGETS space_common DESTINATION lib)
