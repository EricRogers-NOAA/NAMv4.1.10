SHELL=	/bin/bash
LIBS=
FC=	ifort
LDFLAGS=
FFLAGS=	-convert big_endian -assume byterecl
OBJS=	nam_inter_4to6.o interp_coef.o
CMD=	${OUTDIR}/nam_inter_4to6.exe

all:$(CMD)

$(CMD):$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FFLAGS) $<

clean:
	-rm -rf *.o *.lst *.exe

