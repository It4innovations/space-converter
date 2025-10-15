/*
 * Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
 *
 * This program is free software : you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "args_processing.h"
#include <iostream>

namespace space_converter {

	void usage()
	{
		std::cout << "./space_converter --data-type [GADGET, GADGET_SIMPLE, CHANGA_TIPSY, CHANGA_NCHILADA, CSV, GENERICIO, HDF5, HACCBIN] <options> <args>" << std::endl;

		std::cout << "options:" << std::endl;
		std::cout << "\t--grid-dim X" << std::endl;
		std::cout << "\t--output-path X" << std::endl;
		std::cout << "\t--server X" << std::endl;
		std::cout << "\t--port X" << std::endl;
		std::cout << "\t--info" << std::endl;
		std::cout << "\t--nanovdb" << std::endl;
		//std::cout << "\t--dense" << std::endl;
		std::cout << "\t--dense2file" << std::endl;
		std::cout << "\t--anim START END" << std::endl;
		std::cout << "\t--anim-merge" << std::endl;		
		std::cout << "\t--raw-particles" << std::endl;
		std::cout << "\t--rawpart2vdb" << std::endl;
		std::cout << "\t--export-data TYPE DATASET" << std::endl;
		std::cout << "\t--dense-type X" << std::endl;
		std::cout << "\t--bbox-sphere x y z r" << std::endl;
		std::cout << "\t--simple-density" << std::endl;
		std::cout << "\t--offset-position X Y Z" << std::endl;

#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)		
		std::cout << "\t--calc-radius-neigh N" << std::endl;
		std::cout << "\t--calc-radius-neigh-rho-kernel X" << std::endl;
#endif
		std::cout << "\t--calc-radius-neigh-file X" << std::endl;

#ifdef WITH_CUDAKDTREE
		std::cout << "\t--cudakdtree" << std::endl;
		std::cout << "\t--cudakdtree-cpu" << std::endl;
#endif

#ifdef WITH_NANOFLANN
		std::cout << "\t--nanoflann" << std::endl;
#endif
	
		std::cout << "\t--multires" << std::endl;
		std::cout << "\t--bbox x1 y1 z1 x2 y2 z2" << std::endl;
		std::cout << "\t--save-mpi-rank" << std::endl;

		std::cout << "GADGET args:" << std::endl;
		std::cout << "\t--param-file X" << std::endl;
		std::cout << "\t--gadget-file X" << std::endl;
		std::cout << "\t--max-mem-size X" << std::endl;
		std::cout << "\t--buffer-size X" << std::endl;
		std::cout << "\t--part-alloc-factor X" << std::endl;
		std::cout << "\t--bh-count X" << std::endl;

		std::cout << "GADGET_SIMPLE args:" << std::endl;
		std::cout << "\t--gadget-file X" << std::endl;		
		
		std::cout << "CHANGA_TIPSY args:" << std::endl;
		std::cout << "\t--tipsy-file X" << std::endl;

		std::cout << "CHANGA_NCHILADA args:" << std::endl;
		std::cout << "\t--nc-dir X" << std::endl;

		std::cout << "CSV args:" << std::endl;
		std::cout << "\t--csv-file X" << std::endl;

		std::cout << "GENERICIO args:" << std::endl;
		std::cout << "\t--genericio-file X" << std::endl;
		std::cout << "\t--pos-names x y z" << std::endl;
		std::cout << "\t--vel-names vx vy vz" << std::endl;

		std::cout << "HDF5 args:" << std::endl;
		std::cout << "\t--hdf5-file X" << std::endl;

		std::cout << "HACCBIN args:" << std::endl;
		std::cout << "\t--haccbin-file X" << std::endl;

		exit(0);
	}

	void parse_args(FromCL& from_cl, common::SpaceData &space_data, int argc, char** argv)
	{
		if (argc < 3) {
			usage();
		}

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--data-type") {
				from_cl.data_type = argv[++i];
			}
			else if (arg == "--grid-dim") {
				space_data.bbox_dim = std::stoi(argv[++i]);
			}
			else if (arg == "-o" || arg == "--output-path") {
				from_cl.output_path = argv[++i];
			}
			else if (arg == "--server") {
				from_cl.server = argv[++i];
			}
			else if (arg == "--port") {
				from_cl.port = std::stoi(argv[++i]);
			}
			else if (arg == "--info") {
				from_cl.info = true;
				from_cl.remote = false;
			}
#ifdef WITH_OPENVDB
			else if (arg == "--nanovdb") {
				from_cl.use_nanovdb = true;
			}
#endif
			//else if (arg == "--dense") {
			//	//from_cl.use_dense = true;
			//	from_cl.export_extracted_type = 1; // eDense
			//}
			else if (arg == "--dense2file") {
				from_cl.use_dense2file = true;
			}
			else if (arg == "--anim") {
				if(space_data.anim_type == common::SpaceData::AnimType::eNone)
					space_data.anim_type = common::SpaceData::AnimType::eAllPath; // eAllPath

				space_data.anim_start = std::stoi(argv[++i]);
				space_data.anim_end = std::stoi(argv[++i]);
			}
			else if (arg == "--anim-merge") {
				space_data.anim_type = common::SpaceData::AnimType::eAllMerge;
			}
			else if (arg == "--raw-particles") {
				//from_cl.use_raw_particles = true;
				space_data.extracted_type = common::SpaceData::ExtractedType::eParticle; // eRawParticles
			}
			else if (arg == "--save-mpi-rank") {
				from_cl.use_save_mpirank = true;
			}			
			else if (arg == "--rawpart2vdb") {
				from_cl.use_rawpart2vdb = true;
			}			
			else if (arg == "--export-data") {
				space_data.particle_type = std::stoi(argv[++i]);
				space_data.block_name_id = std::stoi(argv[++i]);
				//from_cl.export_dense_type = std::stoi(argv[++i]);				
				from_cl.remote = false;
			}
			else if (arg == "--dense-type") {
				space_data.dense_type = (common::SpaceData::DenseType)std::stoi(argv[++i]);
				if (space_data.dense_type ==  common::SpaceData::DenseType::eNone) //eNone
					space_data.extracted_type = common::SpaceData::ExtractedType::eSparse; // eSparse
				else
					space_data.extracted_type = common::SpaceData::ExtractedType::eDense; // eDense
			}
			else if (arg == "--dense-norm") {
				space_data.dense_norm = (common::SpaceData::DenseNorm)std::stoi(argv[++i]);
			}
			else if (arg == "--bbox-sphere") {
				space_data.use_bbox_sphere = true;
				space_data.bbox_sphere_pos[0] = std::stof(argv[++i]);
				space_data.bbox_sphere_pos[1] = std::stof(argv[++i]);
				space_data.bbox_sphere_pos[2] = std::stof(argv[++i]);
				space_data.bbox_sphere_r = std::stof(argv[++i]);
			}
			else if (arg == "--bbox") {
				//space_data.use_bbox = true;
				space_data.bbox_min[0] = std::stof(argv[++i]);
				space_data.bbox_min[1] = std::stof(argv[++i]);
				space_data.bbox_min[2] = std::stof(argv[++i]);

				space_data.bbox_max[0] = std::stof(argv[++i]);
				space_data.bbox_max[1] = std::stof(argv[++i]);
				space_data.bbox_max[2] = std::stof(argv[++i]);
			}
			else if (arg == "--simple-density") {
				space_data.use_simple_density = true;
				if (space_data.dense_type == common::SpaceData::DenseType::eNone) //eNone
					space_data.dense_type = common::SpaceData::DenseType::eCubic; // eType1

				space_data.extracted_type = common::SpaceData::ExtractedType::eDense; // eDense
			}
			else if (arg == "--offset-position") {
				space_data.offset_position[0] = std::stof(argv[++i]);
				space_data.offset_position[1] = std::stof(argv[++i]);
				space_data.offset_position[2] = std::stof(argv[++i]);
			}
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
			else if (arg == "--calc-radius-neigh") {
				space_data.calc_radius_neigh = std::stoi(argv[++i]);
			}
			else if (arg == "--calc-radius-neigh-rho-kernel") {
				space_data.calc_radius_neigh_rho_kernel = (common::SpaceData::DenseType)std::stoi(argv[++i]);
			}			
#endif			
			else if (arg == "--calc-radius-neigh-file") {
				space_data.calc_radius_neigh_file = argv[++i];
			}
#ifdef WITH_CUDAKDTREE			
			else if (arg == "--cudakdtree") {
				from_cl.use_cudakdtree = true;
			}
			else if (arg == "--cudakdtree-cpu") {
				from_cl.use_cudakdtree_cpu = true;
			}
#endif
#ifdef WITH_NANOFLANN
			else if (arg == "--nanoflann") {
				from_cl.use_nanoflann = true;
			}
#endif
			else if (arg == "--multires") {
				from_cl.use_multires = true;
			}			
			else if (arg == "-h" || arg == "--help") {
				usage();
			}
		}
	}
}

//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\galaxy_merger.dat  --export-data 0 1 0 --port 5000
//--data-type CHANGA_NCHILADA --grid-dim 1000 --output-path e:\temp\changa\tmp\ --nc-dir e:\temp\changa\galaxy_merger.dat.data
//--data-type GADGET --gadget-file f:\temp\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\
//--data-type GADGET --gadget-file f:\temp\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\  --port 5000
//--data-type NANOGRID --grid-dim 1000 --output-path e:\temp\changa\tmp\ --nanogrid-file e:\tmp\space\temp\out_2024-07-31-10310288.nvdb
//--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --export-data 0 1 0
//--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 500 --output-path f:\temp\ --port 5000
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110
//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\galaxy_merger.dat --export-data 0 1 0 --output-path f:\temp\


//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --export-data 0 1 0 --output-path f:\temp\ --info

//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --export-data 0 10 0 --output-path f:\temp\ --info
//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --output-path f:\temp\ --port 5000

//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\accretiondisklowresstd --output-path f:\temp\ --port 5000

//--data-type GADGET --gadget-file e:\temp\black_hole\13_1e4_B1e9_z19\snapdir_{}\snap_{} --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 500 --output-path e:\temp\ --port 5000 --anim 0 1 --bh-count 1
//--data-type GADGET --gadget-file e:\temp\black_hole\13_1e4_B1e9_z19\snapdir_{}\snap_{} --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 500 --output-path e:\temp\ --port 5000 --anim 1000 1001 --bh-count 1
//--data-type GADGET --gadget-file e:\temp\black_hole\13_1e4_B1e9_z19\snapdir_{}\snap_{} --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 500 --output-path e:\temp\ --port 5000 --anim 1000 1001 --bh-count 1 --anim-merge --raw-particles

//--data-type CSV --grid-dim 100 --csv-file e:\temp\csv\data.csv --export-data 0 0 0 --output-path f:\temp\ --info
//--data-type CSV --grid-dim 100 --csv-file e:\temp\csv\data3.csv --output-path f:\temp\  --raw2file

//--data-type GENERICIO --genericio-file e:\data\jar091\down\m000.mpicosmo.499 --grid-dim 1000 --output-path f:\temp\ --info
//--data-type GADGET --gadget-file e:\temp\galaxy\snap_SQ_mitCGM_mbar5e5 --max-mem-size 10000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 100 --output-path f:\temp\ --export-data 0 1 0 --bh-count 4 --info


//--data-type HACCBIN --haccbin-file e:\temp\hacc\Full.cosmo.0 --grid-dim 1000 --output-path f:\temp\ --info
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh 2 --bbox 386.881104 590.247009 613.632996 393.591431 596.957336 620.343323 --export-data 0 1 1
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh-file e:\temp\gadget_rad10.bin

//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh-file e:/temp/gadget_rad10.bin --bbox 386.881104 590.247009 613.632996 393.591431 596.957336 620.343323 --export-data 0 1 1
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh 10 --bbox 386.881104 590.247009 613.632996 393.591431 596.957336 620.343323 --export-data 0 1 1
//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh-file e:\temp\gadget_rad10.bin

//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --port 5000 --bh-count 110 --calc-radius-neigh-file e:\temp\gadget_rad10.bin

//--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --export-data 0 1 1 --bh-count 1 --calc-radius-neigh 10 --multires

//--data-type GADGET --gadget-file e:/temp/galaxy/snap_SQ_mitCGM_mbar5e4 --max-mem-size 100000 --buffer-size 500.0 --part-alloc-factor 1.2 --grid-dim 1000 --output-path f:\temp\ --port 5000 --bh-count 4 --calc-radius-neigh 50

//--data-type GADGET --gadget-file F:\work\blender\SPHtoGrid.jl\snap_sedov --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 10 --output-path f:\temp\ --port 5005

//--data-type GADGET --gadget-file F:\work\blender\SPHtoGrid.jl\snap_sedov --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --port 5005 --bh-count 1 --nanovdb

//--data-type GADGET --gadget-file e:\temp\gadget\very_small_example\snap_081 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path f:\temp\ --port 5005 --bh-count 1 --multires


//--data-type GADGET --gadget-file e:\temp\gadget\notsosmall_example\snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path f:\temp\ --info

//--data-type CHANGA_TIPSY --grid-dim 1000 --output-path e:\temp\changa\tmp\ --tipsy-file e:\temp\changa\tipsy\LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --output-path f:\temp\ --port 5000 --multires
//--data-type GENERICIO --genericio-file e:\temp\hacc\farpoint\m000p-499.select.mpicosmo --grid-dim 1000 --output-path f:\temp\ --port 5000
//--data-type GENERICIO --genericio-file e:\data\jar091\down\m000.mpicosmo.499 --grid-dim 1000 --output-path f:\temp\ --


//--data-type GADGET --gadget-file e:\temp\gadget\alice\snapdir_230\snap_230 --max-mem-size 6000 --buffer-size 150.0 --part-alloc-factor 2.5 --grid-dim 1000 --output-path e:\temp\ --bh-count 1 --export-data 5 1 0
//--data-type GADGET_SIMPLE --gadget-file e:\temp\alice\data\snapdir_1000\snap_1000 --output-path e:\temp\ --port 5000 --calc-radius-neigh 32
//--data-type GADGET_SIMPLE --gadget-file e:\temp\alice\data\snapdir_1000\snap_1000 --output-path e:\temp\ --port 5000
//--data-type GADGET_SIMPLE --gadget-file e:\temp\alice\data\H1e5_DF_p03\snapdir_{}\snap_{} --output-path e:\temp\ --port 5000 --anim 0 3