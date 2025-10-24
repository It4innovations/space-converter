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


cd ~/he_space/

ROOT_DIR=${PWD}

#ml Forge 

# ml HDF5/1.14.0-gompi-2022b

# ml c-blosc/1.21.0-GCC-10.2.0
# ml Boost/1.79.0-GCC-11.3.0
# ml tbb/2021.10.0-GCCcore-12.2.0
# ml VirtualGL/3.1-GCC-12.3.0

# ml Qt5/5.15.7-GCCcore-12.2.0
ml CMake/3.24.3-GCCcore-12.2.0
ml GCC/12.2.0
# ml Mesa/22.2.4-GCCcore-12.2.0
# ml CUDA/12.4.0
# ml OpenMPI/4.1.4-GCC-12.2.0

 ml OpenMPI/4.1.6-NVHPC-24.1-CUDA-12.4.0

#################################################

echo "===================RUN2"

lib_dir=${ROOT_DIR}/install

install=${ROOT_DIR}/install/space_converter_kar
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${lib_dir}/openvdb_kar/lib64
export LD_LIBRARY_PATH=$lib_dir/lib-linux_x64/openvdb/lib:$lib_dir/lib-linux_x64/tbb/lib:$LD_LIBRARY_PATH

out=/scratch/project/open-30-28/OpenGadget/out #${ROOT_DIR}/out
out_t=${out}/remote
mkdir -p ${out_t}

#srun -u -n 256 ${install}/bin/space_converter --data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/tbatalha/PICCOLO/HR/C0/snapdir_015/snap_015 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.25 --grid-dim 100 --output-path ${out} --info
#srun -u -n 64 ${install}/bin/space_converter --data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/L500N2048/snapdir_014/snap_014 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.25 --grid-dim 100 --output-path ${out} --info

#srun -u -n 16 ${install}/bin/space_converter --data-type GADGET --gadget-file /scratch/project/open-30-28/OpenGadget/data/notsosmall_example/snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --port 5000

#srun -u -n 128 ${install}/bin/space_converter --data-type CHANGA_TIPSY --tipsy-file /scratch/project/open-30-28/milanjaros/CHANGA/data/tipsy_small.dat --grid-dim 100 --output-path ${out} --port 5000

#srun -u -n 1 ${install}/bin/space_converter --data-type CHANGA_TIPSY --tipsy-file /scratch/project/open-30-28/milanjaros/CHANGA/data/very_tipsy_small.dat --port 5000 --grid-dim 1000 --output-path ${out_t}

#srun -u -n 16 ${install}/bin/space_converter --data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/Dianoga25x/D7/snapdir_092/snap_092 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --port 5005 #--export-data 0 1

#srun -u -n 4 ${install}/bin/space_converter --data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/very_small_example/snap_081 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --port 5000 #--export-data 0 1

#out=/scratch/project/open-30-28/milanjaros/unimem/galaxy3300_dense
out=/scratch/project/open-30-28/milanjaros/temp_space
mkdir -p ${out}

#srun -u -n 256 ${install}/bin/space_converter --data-type CHANGA_NCHILADA --nc-dir /scratch/project/open-30-28/milanjaros/space_test/CHANGA/data/cosmological_box/cosmo25_2304/cosmo25.2304g.nc --grid-dim 4000 --output-path ${out} --export-data 0 2 0
#-rw-rw----+ 1 milanjaros open-30-28 53M May 31 22:07 /scratch/project/open-30-28/milanjaros/unimem/dim500/Gas_Vel.vdb

#--calc-radius-neigh 10 #
# srun -u -n 8 --overlap --gres=gpu:8 /mnt/proj3/open-28-64/blender/he_space/src/space_converter/scripts/space_converter/karolina/gpu.sh ${install}/bin/space_converter --data-type CHANGA_TIPSY --grid-dim 3300 --output-path ${out} --tipsy-file /scratch/project/open-30-28/milanjaros/space_test/CHANGA/tipsy/LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --calc-radius-neigh 100 --export-data 0 2 6 --bbox 448.632935 519.566467 461.891541 495.680908 566.614441 508.939514 #--export-data 0 17 1 --bbox 448.632935 519.566467 461.891541 495.680908 566.614441 508.939514
# rank: 0: find bbox local: 0.035389
# rank: 0: find bbox mpi: 0.038704, box_size: 1695.000000
# rank: 0: grid convert: 515.509036
# rank: 0: bbox-sphere coord: -0.132470, -0.133274, 0.779303, bbox-sphere radius: 39.873158
# rank: 0: bbox coord: 448.632935, 519.566467, 461.891541, 495.680908, 566.614441, 508.939514
# rank: 0: find minmax mpi: 181.208407, particles_count: 39841865
# rank: 0: merged time: 276.559265
# rank: 0: minI: 3.596327e-12, maxI: 3.204487e+03
# rank: 0: final grid time: 627.343697
# finished: /scratch/project/open-30-28/milanjaros/unimem/galaxy3300_dense/Gas_Vel.vdb


#448.632935 519.566467 461.891541 495.680908 566.614441 508.939514
#rank: 0: bbox coord: 456.118256, 526.312317, 468.784668, 489.280914, 559.474915, 501.947327

#srun -u -n 256 -c 128 ${install}/bin/space_converter --data-type CHANGA_TIPSY --tipsy-file /scratch/project/open-30-28/milanjaros/CHANGA2/data/galaxy_merger/merger64Xb0/LOW_64XLARGEMHDVERTDENSWend64FBSB1mergb0.00410  --grid-dim 100 --output-path ${out} --port 5000


#srun -u -n 1024 -c 8 /mnt/proj3/open-28-64/blender/he_space/src/vmtouch/vmtouch -e /scratch/project/open-30-28/milanjaros/CHANGA2/data/cosmological_box/cosmo25_2304/cosmo25.2304g.nc

# for i in 1 2 4 8 16 32 64 128
# do
# tasks=$(( 8192 / $i ))
# cores=$(( $i ))
# temp_log=/scratch/project/open-30-28/milanjaros/CHANGA2/tests/log_${tasks}_${cores}
# echo $temp_log
# srun -u -n $tasks -c $cores ${install}/bin/space_converter --data-type CHANGA_NCHILADA --nc-dir /scratch/project/open-30-28/milanjaros/CHANGA2/data/cosmological_box/cosmo25_2304/cosmo25.2304g.nc --grid-dim 500 --export-data 0 2 --output-path ${out} >> ${temp_log}
# ls -lah /scratch/project/open-30-28/OpenGadget/out/Gas_Vel.vdb >> ${temp_log}
# done

# for i in 32 16 8 4 2 1 
# do
# tasks=$(( $i ))
# cores=128 #$(( $i ))
# temp_log=/scratch/project/open-30-28/milanjaros/CHANGA2/tests/log_${tasks}_${cores}
# echo $temp_log
# srun -u -n $tasks -c $cores ${install}/bin/space_converter --data-type CHANGA_NCHILADA --nc-dir /scratch/project/open-30-28/milanjaros/CHANGA2/data/cosmological_box/cosmo25_2304/cosmo25.2304g.nc --grid-dim 500 --export-data 0 2 --output-path ${out} >> ${temp_log}
# ls -lah /scratch/project/open-30-28/OpenGadget/out/Gas_Vel.vdb >> ${temp_log}
# done

#export example_snap="--data-type CHANGA_TIPSY --grid-dim 1000 --output-path ${out} --tipsy-file /mnt/proj3/open-30-28/CHANGA/tipsy/testgalaxyN32-mpi.smp.gcc.ompi-meluxina-p320000-b24-N32-.000050 --grid-dim 100 --port 5000"
#export example_snap="--data-type CHANGA_TIPSY --grid-dim 1000 --output-path ${out} --tipsy-file /mnt/proj3/open-30-28/CHANGA/tipsy/LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --grid-dim 100 --port 5000"

#export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/very_small_example/snap_081 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --port 5000 --bh-count 101 --raw-particles"
# export example_snap="--data-type GADGET --gadget-file /scratch/project/open-30-28/OpenGadget/data/notsosmall_example/snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 100 --output-path ${out} --port 5000 --bh-count 101 --calc-radius-neigh 10" # --raw-particles --calc-radius-neigh 10

#export example_snap="--data-type GADGET_SIMPLE --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/very_small_example/snap_081 --grid-dim 100 --output-path ${out} --info --calc-radius-neigh 10 --cudakdtree"
#srun --overlap -u -n 64 --gres=gpu:8 /mnt/proj3/open-28-64/blender/he_space/src/space_converter/scripts/space_converter/karolina/gpu.sh ${install}/bin/space_converter ${example_snap}

# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e4_B1e9_z19/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles"
# srun -u -n 1400 ${install}/bin/space_converter ${example_snap}

#export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e5_B1e9_z19/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles"
#srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}

#export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e6_B1e9_z19_DF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1"

#export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e7_B1e9_z19_DF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles"
#srun -u -n 700 ${install}/bin/space_converter ${example_snap}

#srun -u -n 16 ${install}/bin/space_converter --data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/Dianoga25x/D7/snapdir_092/snap_092 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --export-data 0 1 0 --bh-count 101 #--info

#srun --overlap -u -n 1 --gpus 8 ${install}/bin/space_converter --data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/galaxy/snap_SQ_mitCGM_mbar5e4 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --export-data 0 1 0 --bh-count 100 --info
#srun --overlap -u -n 1 --gpus 8 ${install}/bin/space_converter --data-type GADGET --gadget-file /mnt/proj3/open-30-28/OpenGadget/data/galaxy/snap_SQ_mitCGM_mbar5e5 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --export-data 0 1 0 --bh-count 100 --info



#srun --overlap -u -n 8 --gpus 64 --gres=gpu:8 ${install}/bin/space_converter --data-type GENERICIO --genericio-file /mnt/proj3/open-30-28/milanjaros/HACC/MiraTitanU/Grid/M000/L2100/HACC000/analysis/Particles/STEP499/m000.mpicosmo.499 --grid-dim 1000 --output-path ${out} --port 5000 --anim 0 8

#srun --overlap -u -n 8 --gpus 64 --gres=gpu:8 ${install}/bin/space_converter --data-type GENERICIO --genericio-file /mnt/proj3/open-30-28/milanjaros/HACC/QContinuum/M000/L1300/HAC000/analysis/Halos/b0168/FOFp/STEP499/m000-499.fofproperties --grid-dim 1000 --output-path ${out} --port 5000 --anim 48 77 --pos fof_halo_mean_x,fof_halo_mean_y,fof_halo_mean_z --vel fof_halo_mean_vx,fof_halo_mean_vy,fof_halo_mean_vz # --info

#export GENERICIO_USE_MPIIO=1
#srun -u -n 1 ${install}/bin/genericio_print --no-data /mnt/proj3/open-30-28/milanjaros/HACC/QContinuum/M000/L1300/HAC000/analysis/Particles/STEP499/m000.mpicosmo.499

#/mnt/proj3/open-30-28/milanjaros/HACC/QContinuum/M000/L1300/HAC000/analysis/Particles/STEP499/m000.mpicosmo.499
#srun --overlap -u -n 1 
#srun --overlap -u -n 1024 ${install}/bin/space_converter --data-type GENERICIO --genericio-file /scratch/project/open-30-28/milanjaros/Particles/STEP499/m000.mpicosmo.499 --grid-dim 1000 --output-path ${out} --port 5000 #--info

#srun --overlap -u -n 1024 ${install}/bin/space_converter --data-type GENERICIO --genericio-file /scratch/project/open-30-28/milanjaros/HACC/MiraTitanU/Grid/M010/L2100/HACC007/analysis/Halos/b0168/Halopart/STEP300/m010-300.bighaloparticles  --grid-dim 1000 --output-path ${out} --port 5000

#srun --overlap -u -n 1024 ${install}/bin/space_converter --data-type GENERICIO --genericio-file /scratch/project/open-30-28/milanjaros/HACC/QContinuum/M000/L1300/HAC000/analysis/Halos/b0168/FOFp/STEP198/m000-198.fofproperties --grid-dim 1000 --output-path ${out} --port 5000 --pos fof_halo_center_x,fof_halo_center_y,fof_halo_center_z --vel fof_halo_mean_vx,fof_halo_mean_vy,fof_halo_mean_vz

#--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 500 --output-path f:\temp\ --port 5005 --nanovdb

srun --overlap -u -n 4096 ${install}/bin/space_converter --data-type CHANGA_TIPSY --grid-dim 1000 --output-path ${out} --tipsy-file /mnt/proj3/open-30-28/CHANGA/tipsy/LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --grid-dim 100 --port 5000 #--nanovdb


#/mnt/proj3/open-30-28/milanjaros/HACC/SciVIsContest2019DATA/Full.cosmo.0

# srun --overlap -u -n 625 ${install}/bin/space_converter --data-type HACCBIN --haccbin-file /mnt/proj3/open-30-28/milanjaros/HACC/SciVIsContest2019DATA/Full.cosmo.{} --output-path ${out_t} --port 5000 --anim 0 624 --nanovdb

#srun -u -n 32 ${install}/bin/space_converter --data-type GENERICIO --genericio-file /scratch/project/open-30-28/milanjaros/HACC/Farpoint/M000p/L1478/HACC000/analysis/Particles/STEP499/m000p-499.select.mpicosmo --grid-dim 1000 --output-path ${out} --pos-names x y z --vel-names vx vy vz --port 5000 --info #--nanovdb

###############################BH
# 13_1e4_B1e9_z19  13_1e4_B1e9_z19_DF  13_1e5_B1e9_z19  13_1e5_B1e9_z19_DF  13_1e6_B1e9_z19_DF  13_1e6_B1e9_z19_NODF 13_1e7_B1e9_z19_DF
# Centers_1e5_DFZ19.txt
#export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e5_B1e9_z19_DF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles"
#srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}

# # Centers_1e5_Z19_NODF.txt
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e5_B1e9_z19/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles"
# srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}

# # Centers_1e6_DFZ19.txt
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e6_B1e9_z19_DF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles"
# srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}

# # Centers_1e6_Z19_NODF.txt
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e6_B1e9_z19_NODF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles"
# srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}

###############################DM
# 13_1e4_B1e9_z19  13_1e4_B1e9_z19_DF  13_1e5_B1e9_z19  13_1e5_B1e9_z19_DF  13_1e6_B1e9_z19_DF  13_1e6_B1e9_z19_NODF 13_1e7_B1e9_z19_DF
# Centers_1e5_DFZ19.txt
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e5_B1e9_z19_DF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --simple-density --export-data 1 2 0 --raw-particles  --bbox-sphere 0 0 0 50" # --export-data 1 2 0 --raw-particles  --bbox-sphere 0 0 0 50

# export fix_filename="/mnt/proj3/open-30-28/adamiano/Centers_1e5_DFZ19.txt"

# srun --overlap -u -n 1400 ${src}/space-converter/scripts/space_converter/karolina/bh_offset.sh ${install}/bin/space_converter ${example_snap}

# i=10
# if [ "$i" -lt 1000 ]; then
#     printf -v formatted "%03d" "$i"
# else
#     formatted="$i"
# fi
# #echo "$formatted"
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e5_B1e9_z19_DF/snapdir_${formatted}/snap_${formatted} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --bh-count 1 --raw-particles --bbox-sphere 0 0 0 50 --export-data 1 2 0" # --bbox-sphere 0 0 0 10 --anim 0 1400 --anim-merge
# srun --overlap -u -n 1 ${install}/bin/space_converter ${example_snap}

# # Centers_1e5_Z19_NODF.txt
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e5_B1e9_z19/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles --bbox-sphere 0 0 0 50"
# srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}

# # Centers_1e6_DFZ19.txt
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e6_B1e9_z19_DF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles --bbox-sphere 0 0 0 50"
# srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}

# # Centers_1e6_Z19_NODF.txt
# export example_snap="--data-type GADGET --gadget-file /mnt/proj3/open-30-28/adamiano/halos_DM/13_1e6_B1e9_z19_NODF/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 --bh-count 1 --anim-merge --raw-particles --bbox-sphere 0 0 0 50"
# srun --overlap -u -n 1400 ${install}/bin/space_converter ${example_snap}