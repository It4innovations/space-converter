#!/bin/bash

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

ml purge

ml cmake
ml cuda/11.8
ml openmpi

ROOT_DIR=${PWD}

echo "===================RUN"

lib_dir=${ROOT_DIR}/install
install=${ROOT_DIR}/install/space_converter
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${lib_dir}/openvdb/lib64

srun -u ${install}/bin/space_converter --data-type CHANGA_NCHILADA --nc-dir ${cosmo25_nc_path}

