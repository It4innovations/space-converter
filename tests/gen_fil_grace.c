/*
 * Write a tiny synthetic GRACE volume output (FIL_GRACE reader; GRACE
 * src/IO/hdf5_output.cpp) for the smoke tests.
 *
 * Build: h5cc -o gen_fil_grace gen_fil_grace.c      Usage: gen_fil_grace OUT_DIR
 *
 * Mesh (cell centred, blocks of NB^3 = 4^3 cells, z >= 0 only as in a run with
 * reflection symmetry in z):
 *   level 0: dx 1, the four blocks [-4, 0] / [0, 4] in x and y, z in [0, 4];
 *            the block x < 0, y < 0 is refined
 *   level 1: dx 0.5, its eight children (blocks of 2^3 in size)
 * Blocks in the file: 3 x level 0 (x > 0 y < 0, x < 0 y > 0, x > 0 y > 0), then
 * 8 x level 1: 11 blocks, 704 cells. Mirror z: 1408; levels 0..0: 192; levels
 * 1..1: 512; region (0.5 0.5 0) .. (4 4 4): the block x > 0, y > 0 only: 64.
 * Datasets:
 *   rho  = 1 + level + 0.01 x      alp = 1 - 0.05 level      Bvec = (0.3, 0.4, 0)
 *   (|Bvec| = 0.5), plus Points, Cells and (first file only) Level, Rank.
 * Files: volume_out_000004.h5 (Iteration 4, Time 1, with /Level and /Rank) and
 * volume_out_000008.h5 (Iteration 8, Time 2, without them: the reader then
 * derives the levels from the cell sizes).
 */

#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NB 4
#define NQ 11
#define NC (NB * NB * NB)
#define NP ((NB + 1) * (NB + 1) * (NB + 1))

typedef struct {
	double lo[3];
	double dx;
	int level;
} Block;

static void write_set(hid_t f, const char* name, hid_t ftype, hid_t mtype, int rank, const hsize_t* dims, const void* data) {
	hid_t sp = H5Screate_simple(rank, dims, NULL);
	hid_t pl = H5Pcreate(H5P_DATASET_CREATE);
	hsize_t chunk[2] = { dims[0] < 256 ? dims[0] : 256, rank > 1 ? dims[1] : 1 };
	H5Pset_chunk(pl, rank, chunk);
	hid_t ds = H5Dcreate2(f, name, ftype, sp, H5P_DEFAULT, pl, H5P_DEFAULT);
	H5Dwrite(ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	H5Dclose(ds);
	H5Pclose(pl);
	H5Sclose(sp);
}

static void write_file(const char* path, const Block* b, unsigned int iteration, double time, int extras) {
	hid_t f = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (f < 0) {
		fprintf(stderr, "cannot create %s\n", path);
		exit(1);
	}
	hid_t sp = H5Screate(H5S_SCALAR);
	hid_t a = H5Acreate2(f, "Time", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, H5T_NATIVE_DOUBLE, &time);
	H5Aclose(a);
	a = H5Acreate2(f, "Iteration", H5T_NATIVE_UINT, sp, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, H5T_NATIVE_UINT, &iteration);
	H5Aclose(a);
	H5Sclose(sp);

	double* points = malloc(sizeof(double) * NQ * NP * 3);
	unsigned long long* cells = malloc(sizeof(unsigned long long) * NQ * NC * 8);
	double* rho = malloc(sizeof(double) * NQ * NC);
	double* alp = malloc(sizeof(double) * NQ * NC);
	double* bvec = malloc(sizeof(double) * NQ * NC * 3);
	unsigned int* level = malloc(sizeof(unsigned int) * NQ * NC);
	unsigned int* rank = malloc(sizeof(unsigned int) * NQ * NC);
	static const int vert[8][3] = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} };

	for (int q = 0; q < NQ; q++) {
		for (int k = 0; k <= NB; k++)
			for (int j = 0; j <= NB; j++)
				for (int i = 0; i <= NB; i++) {
					size_t p = (size_t)i + (NB + 1) * (j + (NB + 1) * ((size_t)k + (NB + 1) * q));
					points[3 * p + 0] = b[q].lo[0] + i * b[q].dx;
					points[3 * p + 1] = b[q].lo[1] + j * b[q].dx;
					points[3 * p + 2] = b[q].lo[2] + k * b[q].dx;
				}
		for (int k = 0; k < NB; k++)
			for (int j = 0; j < NB; j++)
				for (int i = 0; i < NB; i++) {
					size_t c = (size_t)i + NB * (j + NB * ((size_t)k + NB * q));
					for (int v = 0; v < 8; v++)
						cells[8 * c + v] = (size_t)(i + vert[v][0]) + (NB + 1) * ((j + vert[v][1]) + (NB + 1) * ((size_t)(k + vert[v][2]) + (NB + 1) * q));
					double x = b[q].lo[0] + (i + 0.5) * b[q].dx;
					rho[c] = 1.0 + b[q].level + 0.01 * x;
					alp[c] = 1.0 - 0.05 * b[q].level;
					bvec[3 * c + 0] = 0.3;
					bvec[3 * c + 1] = 0.4;
					bvec[3 * c + 2] = 0.0;
					level[c] = (unsigned int)b[q].level;
					rank[c] = 0;
				}
	}

	hsize_t dp[2] = { NQ * NP, 3 }, dc[2] = { NQ * NC, 8 }, ds[1] = { NQ * NC }, dv[2] = { NQ * NC, 3 };
	write_set(f, "/Points", H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, 2, dp, points);
	write_set(f, "/Cells", H5T_STD_U64LE, H5T_NATIVE_ULLONG, 2, dc, cells);
	write_set(f, "/rho", H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, 1, ds, rho);
	write_set(f, "/alp", H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, 1, ds, alp);
	write_set(f, "/Bvec", H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, 2, dv, bvec);
	if (extras) {
		write_set(f, "/Level", H5T_NATIVE_UINT, H5T_NATIVE_UINT, 1, ds, level);
		write_set(f, "/Rank", H5T_NATIVE_UINT, H5T_NATIVE_UINT, 1, ds, rank);
	}
	H5Fclose(f);
	free(points); free(cells); free(rho); free(alp); free(bvec); free(level); free(rank);
}

int main(int argc, char** argv) {
	if (argc < 2) {
		fprintf(stderr, "usage: gen_fil_grace OUT_DIR\n");
		return 1;
	}
	char cmd[4096], path[4096];
	snprintf(cmd, sizeof(cmd), "mkdir -p '%s'", argv[1]);
	if (system(cmd) != 0)
		return 1;

	Block b[NQ];
	int q = 0;
	/* level 0: all but the block x < 0, y < 0 */
	for (int j = 0; j < 2; j++)
		for (int i = 0; i < 2; i++) {
			if (i == 0 && j == 0)
				continue;
			b[q].lo[0] = -4.0 + 4.0 * i;
			b[q].lo[1] = -4.0 + 4.0 * j;
			b[q].lo[2] = 0.0;
			b[q].dx = 1.0;
			b[q].level = 0;
			q++;
		}
	/* level 1: the children of the block x < 0, y < 0 */
	for (int k = 0; k < 2; k++)
		for (int j = 0; j < 2; j++)
			for (int i = 0; i < 2; i++) {
				b[q].lo[0] = -4.0 + 2.0 * i;
				b[q].lo[1] = -4.0 + 2.0 * j;
				b[q].lo[2] = 2.0 * k;
				b[q].dx = 0.5;
				b[q].level = 1;
				q++;
			}

	snprintf(path, sizeof(path), "%s/volume_out_000004.h5", argv[1]);
	write_file(path, b, 4, 1.0, 1);
	snprintf(path, sizeof(path), "%s/volume_out_000008.h5", argv[1]);
	write_file(path, b, 8, 2.0, 0);
	printf("wrote %s/volume_out_00000{4,8}.h5: %d blocks of %d^3 cells\n", argv[1], NQ, NB);
	return 0;
}
