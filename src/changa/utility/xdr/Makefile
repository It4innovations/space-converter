#
#  Makefile for custom-compile XDR
#
#CC = xlc
CC_FLAGS=-O

LIBTARGETS = xdr.o xdr_stdio.o xdr_float.o xdr_mem.o htonl.o

libxdr.a: $(LIBTARGETS)
	-rm libxdr.a
	$(AR) q libxdr.a $(LIBTARGETS)
	ranlib libxdr.a

clean:
	rm -f $(LIBTARGETS)

xdr.o: xdr.c xdr.h types.h
	$(CC) $(CC_FLAGS) -c xdr.c

xdr_stdio.o: xdr_stdio.c xdr.h types.h
	$(CC) $(CC_FLAGS) -c xdr_stdio.c

xdr_float.o: xdr_float.c xdr.h types.h
	$(CC) $(CC_FLAGS) -c xdr_float.c

xdr_mem.o: xdr_mem.c
	$(CC) $(CC_FLAGS) -c xdr_mem.c

htonl.o: htonl.c 
	$(CC) $(CC_FLAGS) -c htonl.c

