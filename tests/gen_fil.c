/*
 * Write a tiny synthetic Carpet HDF5 data set (FIL / Einstein Toolkit
 * CarpetIOHDF5 3D output) for the smoke tests.
 *
 * Build: h5cc -o gen_fil gen_fil.c      Usage: gen_fil OUT_DIR
 *
 * Grid (vertex centred, NGH = 2 ghost zones):
 *   level 0: delta 1, points -4..4 per axis (9^3), split into two components
 *            along x at the process boundary x = 0.5: c=0 owns x = -4..0,
 *            c=1 owns x = 1..4; each carries NGH = 2 ghost points across that
 *            face (cctk_bbox = 0 there), i.e. c=0 is stored for x = -4..2 and
 *            c=1 for x = -1..4
 *   level 1: delta 0.5, points -2..2 per axis (9^3), one component, all faces
 *            refinement boundaries (cctk_bbox = 1)
 * Finest level wins: a level-0 point is dropped when its cell [x - 0.5, x + 0.5]^3
 * lies inside the level-1 region [-2.25, 2.25]^3, i.e. x, y, z in {-1, 0, 1}:
 *   27 points -> 729 - 27 + 729 = 1431 points; levels 0..0: 729;
 *   ghost zones kept: level 0 has 13 x 81 = 1053 stored points, 54 covered -> 999 + 729 = 1728.
 * Variables, iterations 0 and 4 (time 0 and 1):
 *   HYDROBASE::rho     = (1 + level) * (1 + it) + 0.01 x      (hydrobase-rho.file_0/1.h5: per process)
 *   HYDROBASE::vel[i]  = (0.3, 0.4, 0)                        (hydrobase-vel.h5: one file per group, |vel| = 0.5)
 *   ADMBASE::alp       = 1 - 0.05 level                       (admbase-lapse.h5)
 * plus hydrobase-rho.xy.h5 (2D output with the same dataset names; must be skipped).
 */

#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NGH 2

typedef struct {
	int level, comp;
	int ilo, ihi;      /* stored x index range (level-0 index space), inclusive */
	int bbox[6];
	double delta;
	int lo, hi;        /* stored y/z (and x for level 1) range, inclusive, in units of delta */
} Patch;

static void attr_d(hid_t ds, const char* name, const double* v, hsize_t n) {
	hid_t sp = H5Screate_simple(1, &n, NULL);
	hid_t a = H5Acreate2(ds, name, H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, H5T_NATIVE_DOUBLE, v);
	H5Aclose(a);
	H5Sclose(sp);
}

static void attr_i(hid_t ds, const char* name, const int* v, hsize_t n) {
	hid_t sp = H5Screate_simple(1, &n, NULL);
	hid_t a = H5Acreate2(ds, name, H5T_NATIVE_INT, sp, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(a, H5T_NATIVE_INT, v);
	H5Aclose(a);
	H5Sclose(sp);
}

/* var: 0 rho, 1..3 vel[0..2], 4 alp */
static double value(int var, int level, int it, double x) {
	switch (var) {
	case 0: return (1.0 + level) * (1.0 + it) + 0.01 * x;
	case 1: return 0.3;
	case 2: return 0.4;
	case 3: return 0.0;
	default: return 1.0 - 0.05 * level;
	}
}

static void write_patch(hid_t file, const char* varname, int var, int it, const Patch* p, int rank2d) {
	const int nx = p->ihi - p->ilo + 1, ny = p->hi - p->lo + 1, nz = rank2d ? 1 : ny;
	double* buf = malloc(sizeof(double) * nx * ny * nz);
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++)
				buf[i + nx * (j + ny * k)] = value(var, p->level, it, (p->ilo + i) * p->delta);

	char name[256];
	if (p->comp >= 0)
		snprintf(name, sizeof(name), "%s it=%d tl=0 rl=%d c=%d", varname, it, p->level, p->comp);
	else
		snprintf(name, sizeof(name), "%s it=%d tl=0 rl=%d", varname, it, p->level);
	hsize_t dims3[3] = { (hsize_t)nz, (hsize_t)ny, (hsize_t)nx };
	hsize_t dims2[2] = { (hsize_t)ny, (hsize_t)nx };
	hid_t sp = rank2d ? H5Screate_simple(2, dims2, NULL) : H5Screate_simple(3, dims3, NULL);
	hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

	const double origin[3] = { p->ilo * p->delta, p->lo * p->delta, p->lo * p->delta };
	const double delta[3] = { p->delta, p->delta, p->delta };
	const int iorigin[3] = { p->ilo, p->lo, p->lo };
	const int ngh[3] = { NGH, NGH, NGH };
	const double time = it * 0.25;
	attr_d(ds, "origin", origin, rank2d ? 2 : 3);
	attr_d(ds, "delta", delta, rank2d ? 2 : 3);
	attr_i(ds, "iorigin", iorigin, rank2d ? 2 : 3);
	attr_i(ds, "level", &p->level, 1);
	attr_d(ds, "time", &time, 1);
	attr_i(ds, "cctk_nghostzones", ngh, 3);
	attr_i(ds, "cctk_bbox", p->bbox, 6);

	H5Dclose(ds);
	H5Sclose(sp);
	free(buf);
}

static hid_t create(const char* dir, const char* name) {
	char path[1024];
	snprintf(path, sizeof(path), "%s/%s", dir, name);
	hid_t f = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (f < 0) {
		fprintf(stderr, "cannot create %s\n", path);
		exit(1);
	}
	/* Carpet writes this group into every file; the reader must ignore it */
	hid_t g = H5Gcreate2(f, "Parameters and Global Attributes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(g);
	return f;
}

int main(int argc, char** argv) {
	if (argc < 2) {
		fprintf(stderr, "usage: %s OUT_DIR\n", argv[0]);
		return 1;
	}
	const char* dir = argv[1];
	char cmd[1100];
	snprintf(cmd, sizeof(cmd), "mkdir -p '%s'", dir);
	if (system(cmd) != 0)
		return 1;

	/* level 0: two components with ghost zones across x = 0.5 */
	const Patch c0 = { 0, 0, -4, 0 + NGH, { 1, 0, 1, 1, 1, 1 }, 1.0, -4, 4 };
	const Patch c1 = { 0, 1, 1 - NGH, 4, { 0, 1, 1, 1, 1, 1 }, 1.0, -4, 4 };
	/* level 1: one component (no c= in the dataset name), refinement boundary everywhere */
	const Patch f0 = { 1, -1, -4, 4, { 1, 1, 1, 1, 1, 1 }, 0.5, -4, 4 };
	const int its[2] = { 0, 4 };

	/* rho: per process files (c0 + level 1 on process 0, c1 on process 1) */
	hid_t r0 = create(dir, "hydrobase-rho.file_0.h5");
	hid_t r1 = create(dir, "hydrobase-rho.file_1.h5");
	hid_t vel = create(dir, "hydrobase-vel.h5");
	hid_t alp = create(dir, "admbase-lapse.h5");
	hid_t xy = create(dir, "hydrobase-rho.xy.h5");
	const char* vel_names[3] = { "HYDROBASE::vel[0]", "HYDROBASE::vel[1]", "HYDROBASE::vel[2]" };
	for (int n = 0; n < 2; n++) {
		const int it = its[n];
		write_patch(r0, "HYDROBASE::rho", 0, it, &c0, 0);
		write_patch(r0, "HYDROBASE::rho", 0, it, &f0, 0);
		write_patch(r1, "HYDROBASE::rho", 0, it, &c1, 0);
		for (int c = 0; c < 3; c++) {
			write_patch(vel, vel_names[c], 1 + c, it, &c0, 0);
			write_patch(vel, vel_names[c], 1 + c, it, &c1, 0);
			write_patch(vel, vel_names[c], 1 + c, it, &f0, 0);
		}
		write_patch(alp, "ADMBASE::alp", 4, it, &c0, 0);
		write_patch(alp, "ADMBASE::alp", 4, it, &c1, 0);
		write_patch(alp, "ADMBASE::alp", 4, it, &f0, 0);
		/* 2D slice output: same dataset names, rank 2 */
		write_patch(xy, "HYDROBASE::rho", 0, it, &c0, 1);
		write_patch(xy, "HYDROBASE::rho", 0, it, &f0, 1);
	}
	H5Fclose(r0);
	H5Fclose(r1);
	H5Fclose(vel);
	H5Fclose(alp);
	H5Fclose(xy);
	printf("wrote synthetic Carpet HDF5 data to %s\n", dir);
	return 0;
}
