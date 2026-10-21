#include <stdio.h>
#include <string.h>

#include "gadget3_io_lib.h"

#include <mpi.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    init_lib();

    if (argc < 3) {
        fprintf(stderr, "usage: %s PARAM_FILE SNAP_FILE\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    char param_file[1024];
    char snap_file[1024];

    snprintf(param_file, sizeof(param_file), "%s", argv[1]);
    snprintf(snap_file, sizeof(snap_file), "%s", argv[2]);

    read_parameter_file(param_file, NULL, NULL, NULL, 0);
    mymalloc_init();
    read_ic(snap_file);

    finish_lib();

    MPI_Finalize();		/* clean up & finalize MPI */

    return 0;
}