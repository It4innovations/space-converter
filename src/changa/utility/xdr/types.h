/* @(#)types.h	2.3 88/08/15 4.0 RPCSRC */
/*
 * Sun RPC is a product of Sun Microsystems, Inc. and is provided for
 * unrestricted use provided that this legend is included on all tape
 * media and as a part of the software program in whole or part.  Users
 * may copy or modify Sun RPC without charge, but are not authorized
 * to license or distribute it to anyone else except as part of a product or
 * program developed by the user.
 *
 * SUN RPC IS PROVIDED AS IS WITH NO WARRANTIES OF ANY KIND INCLUDING THE
 * WARRANTIES OF DESIGN, MERCHANTIBILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE, OR ARISING FROM A COURSE OF DEALING, USAGE OR TRADE PRACTICE.
 *
 * Sun RPC is provided with no support and without any obligation on the
 * part of Sun Microsystems, Inc. to assist in its use, correction,
 * modification or enhancement.
 *
 * SUN MICROSYSTEMS, INC. SHALL HAVE NO LIABILITY WITH RESPECT TO THE
 * INFRINGEMENT OF COPYRIGHTS, TRADE SECRETS OR ANY PATENTS BY SUN RPC
 * OR ANY PART THEREOF.
 *
 * In no event will Sun Microsystems, Inc. be liable for any lost revenue
 * or profits or other special, indirect and consequential damages, even if
 * Sun has been advised of the possibility of such damages.
 *
 * Sun Microsystems, Inc.
 * 2550 Garcia Avenue
 * Mountain View, California  94043
 */
/*      @(#)types.h 1.18 87/07/24 SMI      */

/* fixincludes should not add extern "C" to this file */
/*
 * Rpc additions to <sys/types.h>
 */
#ifndef _RPC_TYPES_H
#define _RPC_TYPES_H 1

typedef int bool_t;
typedef int enum_t;

#define        __dontcare__    -1

#ifndef FALSE
#      define  FALSE   (0)
#endif

#ifndef TRUE
#      define  TRUE    (1)
#endif

#ifndef NULL
#      define  NULL 0
#endif

#include <stdlib.h>		/* For malloc decl.  */
#define mem_alloc(bsize)	malloc(bsize)
#define mem_free(ptr, bsize)	free(ptr)

#if defined __APPLE_CC__ || defined __FreeBSD__
# define __u_char_defined
# define __daddr_t_defined
#endif

#ifdef _WIN32
#ifndef __u_char_defined
typedef unsigned char u_char;
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;
typedef long long quad_t;
typedef unsigned long long u_quad_t;
typedef struct {
    int val[2];
} fsid_t;
# define __u_char_defined
#endif

#ifndef __daddr_t_defined
typedef long daddr_t;
typedef char* caddr_t;
# define __daddr_t_defined
#endif

#endif

#ifndef __u_char_defined
typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;
# define __u_char_defined
#endif
#ifndef __daddr_t_defined
typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;
# define __daddr_t_defined
#endif

#include <sys/types.h>

#ifndef _WIN32
#include <sys/time.h>
#include <sys/param.h>
#endif

/*#include <netinet/in.h>*/

#ifndef INADDR_LOOPBACK
#define       INADDR_LOOPBACK         (u_long)0x7F000001
#endif


#ifndef MAXHOSTNAMELEN
#define        MAXHOSTNAMELEN  64
#endif

#endif /* rpc/types.h */
