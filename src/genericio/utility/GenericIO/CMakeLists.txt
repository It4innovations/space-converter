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

# ------------------------------------------------------------------
# MPI
# ------------------------------------------------------------------
find_package(MPI REQUIRED)

set(INC
    .
    ../SZ/include
    ${BLOSC_INCLUDE_DIRS}
)

set(GENERICIO_LIB_SRC
    GenericIO.cxx
    GenericIO.h
    GenericFileIO.h
    CRC64.h
)

include_directories(${INC})

# ------------------------------------------------------------------
# STATIC IO LIB
# ------------------------------------------------------------------
add_library(genericio_lib_static STATIC
    ${GENERICIO_LIB_SRC}
  )

target_compile_definitions(genericio_lib_static PUBLIC -DGENERICIO_NO_MPI)

target_link_libraries(genericio_lib_static PUBLIC
  sz_lib
  MPI::MPI_CXX
  ${BLOSC_LIBRARIES}
  )

install (TARGETS genericio_lib_static DESTINATION lib)

# ------------------------------------------------------------------
# PRINT
# ------------------------------------------------------------------
add_executable(genericio_print
    GenericIOPrint.cxx
  )

target_compile_definitions(genericio_print PUBLIC -DGENERICIO_NO_MPI)

target_link_libraries(genericio_print PUBLIC
  sz_lib
  genericio_lib_static
  MPI::MPI_CXX
  ${BLOSC_LIBRARIES}
  )

install (TARGETS genericio_print DESTINATION bin)
