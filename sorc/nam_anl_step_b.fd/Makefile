SHELL=	/bin/bash
LIBS=	
FC	=ifort
FFLAGS=	-convert big_endian -assume byterecl -g -traceback
OBJS	=nam_anl_step_b.o Split_Sym_Asy.o correct_mat_1.o interp_coef.o
CMD	= ${OUTDIR}/nam_anl_step_b.exe

$(CMD):$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LIBS)

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FFLAGS) $< 
clean:
	-rm -rf *.o *.exe

