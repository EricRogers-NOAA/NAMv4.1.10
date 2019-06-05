SHELL=	/bin/bash
LIBS=
FC=	ifort
LDFLAGS=
FFLAGS=	-convert big_endian -assume byterecl -g -traceback
OBJS=	nam_merge_nest_b.o interp_coef.o
CMD=	${OUTDIR}/nam_merge_nest_b.exe

all:$(CMD)

$(CMD):$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FFLAGS) $<

clean:
	-rm -rf *.o *.exe
