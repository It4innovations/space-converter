#include <stdio.h>
#include <string.h>

#include "gadget3_io_lib.h"

#include <mpi.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    init_lib();

    char param_file[1024];
    char snap_file[1024];

#if 0
    strcpy(param_file, "e:\\temp\\gadget\\leonardo2\\param.txt");
    strcpy(snap_file, "e:\\temp\\gadget\\leonardo2\\snap_014");
#endif

#if 0
    strcpy(param_file, "e:\\temp\\gadget\\very_small_example\\param.txt");
    strcpy(snap_file, "e:\\temp\\gadget\\very_small_example\\snap_081");
#endif

#if 0
    strcpy(param_file, "/mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/very_small_example/param.txt");
    strcpy(snap_file, "/mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/very_small_example/snap_081");
#endif

#if 0
    strcpy(param_file, "e:\\temp\\gadget\\box_32_L15_gas_extrae_1node_2mpitasks\\param.txt");
    strcpy(snap_file, "e:\\temp\\gadget\\box_32_L15_gas_extrae_1node_2mpitasks\\snapdir_026\\snap_026");
#endif

#if 1
    strcpy(param_file, "/mnt/proj3/open-30-28/OpenGadget/data/L500N2048/param.txt");
    strcpy(snap_file, "/mnt/proj3/open-30-28/OpenGadget/data/L500N2048/snapdir_014/snap_014");
#endif


#if 0
    strcpy(param_file, "/mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/notsosmall_example/param.txt");
    strcpy(snap_file, "/mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/notsosmall_example/snap_091");
#endif

#if 0
    strcpy(param_file, "/scratch/project/open-30-28/OpenGadget/data/L500N2048/snapdir_014/param.txt");
    strcpy(snap_file, "/scratch/project/open-30-28/OpenGadget/data/L500N2048/snapdir_014/snap_014");
#endif

    read_parameter_file(param_file, NULL, NULL, NULL, 0);
    mymalloc_init();
    read_ic(snap_file);

    finish_lib();

    MPI_Finalize();		/* clean up & finalize MPI */

    return 0;
}