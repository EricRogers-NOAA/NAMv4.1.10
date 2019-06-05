SHELL=	/bin/bash
FC	=ifort
FFLAGS	= -convert big_endian -assume byterecl
OBJS	= nam_create_trak_fnl.o
CMD	= ${OUTDIR}/nam_create_trak_fnl.exe

$(CMD):$(OBJS)
	$(FC) -o  $(@) $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FFLAGS) $< 
clean:
	rm -rf *.o *.exe
