#ifndef PRECISION_H
#define PRECISION_H

#ifndef GADGET3_IO_LIB
#include "../gadgetconfig.h"
#endif
#include "../CodeBase/switches.h"


typedef int sort_type;
#define sizesort sizeof(sort_type)

#ifndef LONGIDS
typedef unsigned int MyIDType;
#else
typedef unsigned long long MyIDType;
#endif


#ifndef DOUBLEPRECISION		/* default is single-precision */
typedef float MyFloat;
typedef float MyDouble;
typedef float MyLongDouble;
typedef double MyAtLeastDouble;
#elif (DOUBLEPRECISION+0) == 1	/* everything double-precision, more memory but faster */
typedef double MyFloat;
typedef double MyDouble;
typedef double MyLongDouble;
typedef double MyAtLeastDouble;
#elif (DOUBLEPRECISION+0) == 2	/* everything long double-precision */
typedef long double MyFloat;
typedef long double MyDouble;
typedef long double MyLongDouble;
typedef long double MyAtLeastDouble;
#elif (DOUBLEPRECISION+0) == 3	/* mixed precision, very low memory, moderate precission */
typedef float MyFloat;
typedef float MyDouble;
typedef double MyLongDouble;
typedef double MyAtLeastDouble;
#elif (DOUBLEPRECISION+0) == 4	/* low memory, some extra high precissin for special particle values (like Pos) */
typedef float MyFloat;
typedef double MyDouble;
typedef long double MyLongDouble;
typedef double MyAtLeastDouble;
#elif (DOUBLEPRECISION+0) == 5	/* low memory, extra high precission for local calculations and special particle values */
typedef float MyFloat;
typedef double MyDouble;
typedef long double MyLongDouble;
typedef long double MyAtLeastDouble;
#else /* Not allowd any longer! */
#error DOUBLEPRECISSION has to be a value between 1 and 5 if set !
#endif

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else
typedef float MyOutputFloat;
#endif

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else
typedef float MyInputFloat;
#endif

#endif
