#ifndef GADGET3_IO_LIB_H
#define GADGET3_IO_LIB_H


#ifdef _WIN32

#ifdef GADGET3_IO_LIB
#	define GADGET3_IO_LIB_API __declspec(dllexport)
#else
#	define GADGET3_IO_LIB_API __declspec(dllimport)
#endif

#else
#	define GADGET3_IO_LIB_API
#endif

#ifdef __cplusplus
extern "C"
{
#endif

	GADGET3_IO_LIB_API void init_lib(int world_rank, int world_size);
	GADGET3_IO_LIB_API void finish_lib();
	GADGET3_IO_LIB_API void set_parameter(int ICFormat, int SnapFormat, int NumFilesWrittenInParallel, int MaxMemSize, double BufferSize, double PartAllocFactor, int TotBHs);
	GADGET3_IO_LIB_API void read_parameter_file(char* fname, char* tag[], void** addr, int* id, int nt);
	GADGET3_IO_LIB_API void read_ic(char* fname);

	GADGET3_IO_LIB_API void mymalloc_init();

#ifdef __cplusplus
}
#endif
#endif