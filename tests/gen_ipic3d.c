/*
 * Generate a small synthetic IPIC3D_HDF5 dataset for smoke tests.
 *
 * Usage:   gen_ipic3d <output_dir>
 * Build:   h5cc -o gen_ipic3d gen_ipic3d.c
 *
 * Writes (layout expected by src/ipic3d/ipic3d_hdf5_extract_iolib.cpp):
 *   <dir>/restart0.hdf
 *     /particles/species_0/{x,y,z,u,v,w,q}/cycle_0   1D double [400]  (electrons, q < 0)
 *     /particles/species_0/ID/cycle_0                1D uint64 [400]
 *     /particles/species_1/{x,y,z,u,v,w,q}/cycle_0   1D double [200]  (ions, q > 0)
 *     /particles/species_1/ID/cycle_0                1D uint64 [200]
 *     /fields/{Ex,Ey,Ez,Bx,By,Bz}/cycle_0            3D double [10][10][10]
 *     /moments/species_{0,1}/{rho,pXX,pXY,pXZ,pYY,pYZ,pZZ}/cycle_0
 *                                                    3D double [10][10][10]
 *     /topology/cartesian_coord                      1D int [3] = {0,0,0}
 *   <dir>/settings.hdf
 *     /collective/{Dx,Dy,Dz,Lx,Ly,Lz}                1-element double
 *     /collective/{Nxc,Nyc,Nzc}, /topology/{XLEN,YLEN,ZLEN}  1-element int
 *
 * The 10^3 field/moment tiles follow iPIC3D's node-centered layout with a
 * 1-cell ghost layer on each side (8 interior cells per axis = Nxc/Nyc/Nzc
 * in settings.hdf). The reader turns every grid point into one synthetic
 * "particle" per species, so with both species present it sees
 * 400 + 1000 (species_0) and 200 + 1000 (species_1) particles.
 *
 * All values come from simple deterministic formulas (no rand()): re-running
 * produces identical file contents.
 */

#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#define N_ELECTRONS 400          /* species_0 */
#define N_IONS      200          /* species_1 */
#define Q_ELECTRON  (-0.001)     /* per-particle charge, species_0 */
#define Q_ION       (0.0025)     /* per-particle charge, species_1 (|q| unequal on purpose) */

#define NXC 8                    /* interior cells per axis (settings.hdf Nxc/Nyc/Nzc) */
#define NG  (NXC + 2)            /* stored tile size incl. 1-cell ghost layer = 10 */
#define DX  1.0                  /* grid spacings - unequal so the box is non-cubic */
#define DY  1.5
#define DZ  2.0
#define LX  (NXC * DX)           /* domain size: 8 x 12 x 16 */
#define LY  (NXC * DY)
#define LZ  (NXC * DZ)

/* Write one dataset, creating intermediate groups automatically. */
static void write_dataset(hid_t file, const char *path, hid_t mem_type, hid_t file_type,
                          int rank, const hsize_t *dims, const void *data)
{
    hid_t space = H5Screate_simple(rank, dims, NULL);
    hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl, 1);
    hid_t dset = H5Dcreate2(file, path, file_type, space, lcpl, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        fprintf(stderr, "gen_ipic3d: failed to create dataset %s\n", path);
        exit(1);
    }
    if (H5Dwrite(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
        fprintf(stderr, "gen_ipic3d: failed to write dataset %s\n", path);
        exit(1);
    }
    H5Dclose(dset);
    H5Pclose(lcpl);
    H5Sclose(space);
}

static void write_doubles(hid_t file, const char *path, int rank, const hsize_t *dims,
                          const double *data)
{
    write_dataset(file, path, H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE, rank, dims, data);
}

static void write_uint64s(hid_t file, const char *path, int rank, const hsize_t *dims,
                          const uint64_t *data)
{
    write_dataset(file, path, H5T_NATIVE_UINT64, H5T_STD_U64LE, rank, dims, data);
}

static void write_ints(hid_t file, const char *path, int rank, const hsize_t *dims,
                       const int *data)
{
    write_dataset(file, path, H5T_NATIVE_INT, H5T_STD_I32LE, rank, dims, data);
}

/* 1-element convenience writers for settings.hdf */
static void write_scalar_double(hid_t file, const char *path, double v)
{
    hsize_t one = 1;
    write_doubles(file, path, 1, &one, &v);
}

static void write_scalar_int(hid_t file, const char *path, int v)
{
    hsize_t one = 1;
    write_ints(file, path, 1, &one, &v);
}

/* frac(i * c): cheap deterministic pseudo-uniform sequence in [0, 1).
   (v is always positive and small here, so truncation == floor; this keeps
   the program libm-free and h5cc needs no extra -lm.) */
static double scramble(size_t i, double c)
{
    double v = (double)(i + 1) * c;
    return v - (double)(long long)v;
}

/* Write /particles/species_<sp>/{x,y,z,u,v,w,q,ID}/cycle_0 for n particles. */
static void write_species_particles(hid_t file, int sp, size_t n, double q)
{
    double *buf = (double *)malloc(n * sizeof(double));
    uint64_t *ids = (uint64_t *)malloc(n * sizeof(uint64_t));
    hsize_t dims = (hsize_t)n;
    char path[128];
    size_t i;

    /* positions: golden-ratio-ish scrambles fill the box quasi-uniformly */
    snprintf(path, sizeof(path), "/particles/species_%d/x/cycle_0", sp);
    for (i = 0; i < n; i++) buf[i] = LX * scramble(i, 0.6180339887 + 0.1 * sp);
    write_doubles(file, path, 1, &dims, buf);

    snprintf(path, sizeof(path), "/particles/species_%d/y/cycle_0", sp);
    for (i = 0; i < n; i++) buf[i] = LY * scramble(i, 0.7548776662 + 0.1 * sp);
    write_doubles(file, path, 1, &dims, buf);

    snprintf(path, sizeof(path), "/particles/species_%d/z/cycle_0", sp);
    for (i = 0; i < n; i++) buf[i] = LZ * scramble(i, 0.8191725134 + 0.1 * sp);
    write_doubles(file, path, 1, &dims, buf);

    /* velocities: small, varying, sign-alternating */
    snprintf(path, sizeof(path), "/particles/species_%d/u/cycle_0", sp);
    for (i = 0; i < n; i++) buf[i] = 0.1 * (scramble(i, 0.331) - 0.5);
    write_doubles(file, path, 1, &dims, buf);

    snprintf(path, sizeof(path), "/particles/species_%d/v/cycle_0", sp);
    for (i = 0; i < n; i++) buf[i] = 0.1 * (scramble(i, 0.557) - 0.5);
    write_doubles(file, path, 1, &dims, buf);

    snprintf(path, sizeof(path), "/particles/species_%d/w/cycle_0", sp);
    for (i = 0; i < n; i++) buf[i] = 0.1 * (scramble(i, 0.773) - 0.5);
    write_doubles(file, path, 1, &dims, buf);

    /* constant per-species charge; the reader uses |q| as the mass proxy */
    snprintf(path, sizeof(path), "/particles/species_%d/q/cycle_0", sp);
    for (i = 0; i < n; i++) buf[i] = q;
    write_doubles(file, path, 1, &dims, buf);

    snprintf(path, sizeof(path), "/particles/species_%d/ID/cycle_0", sp);
    for (i = 0; i < n; i++) ids[i] = (uint64_t)(sp * 1000000 + (int)i + 1);
    write_uint64s(file, path, 1, &dims, ids);

    free(ids);
    free(buf);
}

/* Fill one NG^3 node grid from a per-node formula and write it as a 3D dataset. */
typedef double (*grid_formula)(int ix, int iy, int iz, int sp);

static void write_grid(hid_t file, const char *path, grid_formula f, int sp, double *buf)
{
    hsize_t dims[3] = { NG, NG, NG };
    size_t idx = 0;
    int ix, iy, iz;
    for (ix = 0; ix < NG; ix++)
        for (iy = 0; iy < NG; iy++)
            for (iz = 0; iz < NG; iz++)
                buf[idx++] = f(ix, iy, iz, sp);
    write_doubles(file, path, 3, dims, buf);
}

/* Field/moment formulas: smooth, distinct per component, all deterministic. */
static double f_Ex(int ix, int iy, int iz, int sp) { (void)sp; return 0.01 * ix + 0.001 * iy + 0.0001 * iz; }
static double f_Ey(int ix, int iy, int iz, int sp) { (void)sp; return 0.01 * iy + 0.001 * iz + 0.0001 * ix; }
static double f_Ez(int ix, int iy, int iz, int sp) { (void)sp; return 0.01 * iz + 0.001 * ix + 0.0001 * iy; }
static double f_Bx(int ix, int iy, int iz, int sp) { (void)sp; (void)ix; (void)iz; return 0.1 + 0.002 * iy; }
static double f_By(int ix, int iy, int iz, int sp) { (void)sp; (void)ix; (void)iy; return 0.2 + 0.002 * iz; }
static double f_Bz(int ix, int iy, int iz, int sp) { (void)sp; (void)iy; (void)iz; return 0.3 + 0.002 * ix; }
/* charge density: negative for electrons (sp 0), positive for ions (sp 1) */
static double f_rho(int ix, int iy, int iz, int sp)
{
    double base = 0.4 + 0.001 * (ix + iy + iz);
    return (sp == 0) ? -base : base;
}
static double f_pXX(int ix, int iy, int iz, int sp) { (void)iy; (void)iz; return 0.01 + 0.0001 * (ix + sp); }
static double f_pYY(int ix, int iy, int iz, int sp) { (void)ix; (void)iz; return 0.02 + 0.0001 * (iy + sp); }
static double f_pZZ(int ix, int iy, int iz, int sp) { (void)ix; (void)iy; return 0.03 + 0.0001 * (iz + sp); }
static double f_pXY(int ix, int iy, int iz, int sp) { (void)iz; return 0.0001 * (ix - iy) + 0.001 * sp; }
static double f_pXZ(int ix, int iy, int iz, int sp) { (void)iy; return 0.0001 * (ix - iz) + 0.001 * sp; }
static double f_pYZ(int ix, int iy, int iz, int sp) { (void)ix; return 0.0001 * (iy - iz) + 0.001 * sp; }

static void write_species_moments(hid_t file, int sp, double *buf)
{
    char path[128];
#define MOMENT(name, formula) \
    do { \
        snprintf(path, sizeof(path), "/moments/species_%d/" name "/cycle_0", sp); \
        write_grid(file, path, formula, sp, buf); \
    } while (0)
    MOMENT("rho", f_rho);
    MOMENT("pXX", f_pXX);
    MOMENT("pXY", f_pXY);
    MOMENT("pXZ", f_pXZ);
    MOMENT("pYY", f_pYY);
    MOMENT("pYZ", f_pYZ);
    MOMENT("pZZ", f_pZZ);
#undef MOMENT
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s <output_dir>\n", argv[0]);
        return 1;
    }
    const char *out_dir = argv[1];
    if (mkdir(out_dir, 0775) != 0 && errno != EEXIST) {
        fprintf(stderr, "gen_ipic3d: cannot create directory %s\n", out_dir);
        return 1;
    }

    char file_path[1024];
    double *grid_buf = (double *)malloc((size_t)NG * NG * NG * sizeof(double));

    /* ------------------------------------------------ restart0.hdf ------- */
    snprintf(file_path, sizeof(file_path), "%s/restart0.hdf", out_dir);
    hid_t restart = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (restart < 0) {
        fprintf(stderr, "gen_ipic3d: cannot create %s\n", file_path);
        return 1;
    }

    write_species_particles(restart, 0, N_ELECTRONS, Q_ELECTRON);
    write_species_particles(restart, 1, N_IONS, Q_ION);

    write_grid(restart, "/fields/Ex/cycle_0", f_Ex, 0, grid_buf);
    write_grid(restart, "/fields/Ey/cycle_0", f_Ey, 0, grid_buf);
    write_grid(restart, "/fields/Ez/cycle_0", f_Ez, 0, grid_buf);
    write_grid(restart, "/fields/Bx/cycle_0", f_Bx, 0, grid_buf);
    write_grid(restart, "/fields/By/cycle_0", f_By, 0, grid_buf);
    write_grid(restart, "/fields/Bz/cycle_0", f_Bz, 0, grid_buf);

    write_species_moments(restart, 0, grid_buf);
    write_species_moments(restart, 1, grid_buf);

    /* single-file dataset: this tile sits at topology coordinate (0,0,0) */
    {
        int coord[3] = { 0, 0, 0 };
        hsize_t three = 3;
        write_ints(restart, "/topology/cartesian_coord", 1, &three, coord);
    }
    H5Fclose(restart);
    printf("wrote %s: species_0 %d particles (q=%g), species_1 %d particles (q=%g), "
           "grids %dx%dx%d\n", file_path, N_ELECTRONS, Q_ELECTRON, N_IONS, Q_ION, NG, NG, NG);

    /* ------------------------------------------------ settings.hdf ------- */
    snprintf(file_path, sizeof(file_path), "%s/settings.hdf", out_dir);
    hid_t settings = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (settings < 0) {
        fprintf(stderr, "gen_ipic3d: cannot create %s\n", file_path);
        return 1;
    }
    write_scalar_double(settings, "/collective/Dx", DX);
    write_scalar_double(settings, "/collective/Dy", DY);
    write_scalar_double(settings, "/collective/Dz", DZ);
    write_scalar_double(settings, "/collective/Lx", LX);
    write_scalar_double(settings, "/collective/Ly", LY);
    write_scalar_double(settings, "/collective/Lz", LZ);
    write_scalar_int(settings, "/collective/Nxc", NXC);
    write_scalar_int(settings, "/collective/Nyc", NXC);
    write_scalar_int(settings, "/collective/Nzc", NXC);
    write_scalar_int(settings, "/topology/XLEN", 1);
    write_scalar_int(settings, "/topology/YLEN", 1);
    write_scalar_int(settings, "/topology/ZLEN", 1);
    H5Fclose(settings);
    printf("wrote %s: box %gx%gx%g, cells %dx%dx%d, spacing (%g,%g,%g)\n",
           file_path, LX, LY, LZ, NXC, NXC, NXC, DX, DY, DZ);

    free(grid_buf);
    return 0;
}
