#pragma once

#define HAVE_XDR_HYPER
//#undef HAVE_XDR_HYPER

#ifdef _WIN32

#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

typedef unsigned long long u_int64_t;

#endif