set(SZ_Sources
    src/ByteToolkit.c
    src/dataCompression.c
    src/DynamicIntArray.c
    src/iniparser.c
    src/CompressElement.c
    src/DynamicByteArray.c
    src/rw.c
    src/TightDataPointStorageI.c
    src/TightDataPointStorageD.c
    src/TightDataPointStorageF.c
    src/conf.c
    src/DynamicDoubleArray.c
    src/TypeManager.c
    src/dictionary.c
    src/DynamicFloatArray.c
    src/VarSet.c
    src/callZlib.c
    src/Huffman.c
    src/sz_float.c
    src/sz_double.c
    src/sz_int8.c
    src/sz_int16.c
    src/sz_int32.c
    src/sz_int64.c
    src/sz_uint8.c
    src/sz_uint16.c
    src/sz_uint32.c
    src/sz_uint64.c
    src/szd_uint8.c
    src/szd_uint16.c
    src/szd_uint32.c
    src/szd_uint64.c
    src/szd_float.c
    src/szd_double.c
    src/szd_int8.c
    src/szd_int16.c
    src/szd_int32.c
    src/szd_int64.c
    src/utility.c
    src/sz.c
    src/sz_float_pwr.c
    src/sz_double_pwr.c
    src/szd_float_pwr.c
    src/szd_double_pwr.c
    src/sz_double_ts.c
    src/sz_float_ts.c
    src/szd_double_ts.c
    src/szd_float_ts.c

    include/ByteToolkit.h
    include/dataCompression.h
    include/DynamicIntArray.h
    include/iniparser.h
    include/CompressElement.h
    include/DynamicByteArray.h
    include/rw.h
    include/TightDataPointStorageI.h
    include/TightDataPointStorageD.h
    include/TightDataPointStorageF.h
    include/conf.h
    include/DynamicDoubleArray.h
    include/TypeManager.h
    include/dictionary.h
    include/DynamicFloatArray.h
    include/VarSet.h
    include/callZlib.h
    include/Huffman.h
    include/sz_float.h
    include/sz_double.h
    include/sz_int8.h
    include/sz_int16.h
    include/sz_int32.h
    include/sz_int64.h
    include/sz_uint8.h
    include/sz_uint16.h
    include/sz_uint32.h
    include/sz_uint64.h
    include/szd_uint8.h
    include/szd_uint16.h
    include/szd_uint32.h
    include/szd_uint64.h
    include/szd_float.h
    include/szd_double.h
    include/szd_int8.h
    include/szd_int16.h
    include/szd_int32.h
    include/szd_int64.h
    include/utility.h
    include/sz.h
    include/sz_float_pwr.h
    include/sz_double_pwr.h
    include/szd_float_pwr.h
    include/szd_double_pwr.h
    include/sz_double_ts.h
    include/sz_float_ts.h
    include/szd_double_ts.h
    include/szd_float_ts.h
)
set(SZ_include_dirs
    .
    include
    ${ZLIB_INCLUDE_DIRS}
    ${ZSTD_INCLUDE_DIRS}
)

include_directories(${SZ_include_dirs})
add_library(sz_lib STATIC ${SZ_Sources})

# target_link_libraries(
#    sz_lib 
#    ${ZLIB_LIBRARIES}
#    ${ZSTD_LIBRARIES}
# )

#set_property(TARGET sz_lib PROPERTY C_STANDARD 99)
#target_include_directories(sz_lib PRIVATE ${SZ_include_dirs})
#target_include_directories(sz_lib INTERFACE 
#    blosc
#    include
#    )
#target_compile_definitions(blosc PRIVATE HAVE_LZ4 HAVE_SNAPPY HAVE_ZLIB HAVE_ZSTD HAVE_UNISTD_H)
#if(OpenMP_FOUND)
#    target_link_libraries(blosc PRIVATE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
#endif()
