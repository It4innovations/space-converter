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

cmake_minimum_required(VERSION 3.20)

project (space_converter)

########################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(WITH_MPI ON)
set(WITH_OPENMP ON)
set(WITH_NANOVDB ON)

option(WITH_OPENVDB  "Enable OPENVDB" ON)

option(WITH_HDF5  "Enable HDF5" OFF)
option(WITH_HACC  "Enable HACC" OFF)
option(WITH_IPIC3D  "Enable IPIC3D" OFF)
option(WITH_PLUTO  "Enable PLUTO" OFF)
option(WITH_VTK   "Enable VTK" OFF)

option(GADGET_MAX_HSML  "Enable loading of limited values from gadget files" OFF)
option(GADGET_READ_ID  "Enable reading of ID" OFF)
option(GADGET_READ_ID64  "Enable reading of 64 bit ID" OFF)

option(WITH_NANOFLANN  "Enable NANOFLANN" OFF)
option(WITH_CUDAKDTREE "Enable cudaKDTree" OFF)

option(WITH_GPU_CUDA "Enable GPU CUDA" OFF)
option(WITH_CUDA_AWARE_MPI "Enable CUDA-aware MPI for RDMA GPU-direct transfers (requires UCX/RDMA-capable MPI)" OFF)

option(WITH_MERIC  "Enable MERIC" OFF)
########################################################

if(WITH_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    if(WIN32)
        # On Windows, use LLVM OpenMP implementation
        string(APPEND CMAKE_C_FLAGS " /openmp:llvm ")
        string(APPEND CMAKE_CXX_FLAGS " /openmp:llvm /bigobj")
    else()
        string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
        string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
    endif()
    
    string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OpenMP_LINKER_FLAGS}")
    string(APPEND CMAKE_MODULE_LINKER_FLAGS " ${OpenMP_LINKER_FLAGS}")

  else()
    set(WITH_OPENMP OFF)
  endif()
endif()

if(WITH_MERIC)   
    set(MERIC_INCLUDE_DIRS "" CACHE PATH "")
    set(MERIC_LIBRARIES "" CACHE FILEPATH "")
endif()

if(WITH_OPENVDB)   
    set(OPENVDB_INCLUDE_DIRS "" CACHE PATH "")
    set(OPENVDB_LIBRARIES "" CACHE FILEPATH "")
    set(OPENVDB_LIBRARIES_DEBUG "" CACHE FILEPATH "")
    set(WITH_TBB ON)
else()
    set(OPENVDB_INCLUDE_DIRS "")
    set(NANOVDB_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/submodules/openvdb/nanovdb")

    set(WITH_TBB OFF)
endif()

if(WITH_VTK)
    # If you know exactly which modules you need, list them here (recommended).
    find_package(VTK REQUIRED COMPONENTS
        CommonCore
        CommonDataModel
        IOLegacy
        IOXML
        IOParallelXML
    )
endif()

if(WITH_MPI)
    find_package(MPI)
    if(NOT MPI_FOUND)
      set(WITH_MPI OFF)
    endif()

    set(MPI_CUSTOM_LIBRARIES "" CACHE FILEPATH "")
endif()

if(WITH_TBB)
  set(TBB_INCLUDE_DIRS "" CACHE PATH "")
  set(TBB_LIBRARIES "" CACHE FILEPATH "")
  set(TBB_LIBRARIES_DEBUG "" CACHE FILEPATH "")
endif()

if(WITH_NANOVDB)   
    set(NANOVDB_INCLUDE_DIRS "" CACHE PATH "")
endif()

if(WITH_GPU_CUDA OR WITH_CUDAKDTREE)
    # ==================================================================
    # configure and enable cuda
    # ==================================================================
    if (CMAKE_CUDA_ARCHITECTURES)
        # set on the cmdline
        set (CUKD_INIT_ARCHS "${CMAKE_CUDA_ARCHITECTURES}")
    elseif (${CMAKE_VERSION} VERSION_LESS 3.24)
        set (CUKD_INIT_ARCHS "70;80")
    else()
        set (CUKD_INIT_ARCHS "all-major")
    endif()
    set(CUKD_CUDA_ARCHITECTURES "${CUKD_INIT_ARCHS}"
        CACHE STRING "CUDA Arch(s) to build against")
    set(CMAKE_CUDA_ARCHITECTURES ${CUKD_CUDA_ARCHITECTURES})
    
    enable_language(CUDA)
endif()

if(WITH_HACC)
    set(BUILD_SHARED OFF)
    set(BUILD_TESTS OFF)
    set(BUILD_TESTS OFF)
    set(BUILD_BENCHMARKS OFF)
    set(BUILD_FUZZERS OFF)

    add_subdirectory(submodules/c-blosc)

    set(BLOSC_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/submodules/c-blosc/include")
    set(BLOSC_LIBRARIES blosc_static)
    
    # Set zlib and zstd include directories for SZ library
    set(ZLIB_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/submodules/c-blosc/internal-complibs/zlib-1.3.1")
    set(ZSTD_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/submodules/c-blosc/internal-complibs/zstd-1.5.6")
endif()

if(WITH_IPIC3D)
    # IPIC3D requires HDF5
    set(WITH_HDF5 ON)
endif()

if(WITH_PLUTO)
    # PLUTO requires VTK
    set(WITH_VTK ON)
endif()

if(WITH_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS C)
    if(HDF5_FOUND)
        set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
        set(HDF5_LIBRARIES ${HDF5_LIBRARIES})
    endif()
endif()

########################################################
add_subdirectory(src)
########################################################
