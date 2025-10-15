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

set(GADGET3_IO_LIB_MPI ON)

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

# ------------------------------------------------------------------
# MPI
# ------------------------------------------------------------------
if (GADGET3_IO_LIB_MPI)
  find_package(MPI REQUIRED)
endif()

set(GADGET3_IO_LIB_SRC
  IO/io.c
  IO/read_ic.c
  
  IO/gadget3_io_lib.h
  IO/gadget3_io_lib.c

  CodeBase/allvars.c
  CodeBase/allvars.h

  CodeBase/begrun.c
  CodeBase/endrun.c

  CoolingSfr/Sfr_LT/lt_utils.c

  Integrator/timestep.c

  System/mymalloc.c
  System/system.c
  System/allocate.c
)

# ------------------------------------------------------------------
# STATIC IO LIB
# ------------------------------------------------------------------
add_library(gadget3_io_lib_static STATIC
    ${GADGET3_IO_LIB_SRC}
  )

target_compile_definitions(gadget3_io_lib_static PUBLIC -DGADGET3_IO_LIB)

if (GADGET3_IO_LIB_MPI)
  target_link_libraries(gadget3_io_lib_static PUBLIC
    MPI::MPI_CXX
    )
  target_compile_definitions(gadget3_io_lib_static PUBLIC -DGADGET3_IO_LIB_MPI)
endif()

if(WITH_HDF5)
  target_link_libraries(gadget3_io_lib_static PUBLIC
    hdf5::hdf5
    )
  target_compile_definitions(gadget3_io_lib_static PUBLIC -DHAVE_HDF5 -DH5_BUILT_AS_DYNAMIC_LIB)
endif()

install (TARGETS gadget3_io_lib_static DESTINATION lib)
install (FILES IO/gadget3_io_lib.h DESTINATION include)