SHELL   =/bin/bash
EXEC    = $(OUTDIR)/nam_nems_bin_io_degree.exe
FC      =  ifort
FOPTS  = -convert big_endian -assume byterecl -mkl -traceback -I${NEMSIO_INC}
LIBS   = ${NEMSIO_LIB} ${BACIO_LIB4} ${W3EMC_LIB4} ${W3NCO_LIB4}
OBJS = interp_coef_bgrid.o nems_bin_io_degree.o
SRCS = interp_coef_bgrid.f90 nems_bin_io_degree.f90
# *************************************************************************
$(EXEC): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) $(EXEC) $(OBJS) *.mod

#$(OBJS): $(SRCS)
#	$(FC)  $(FOPTS) $(SRCS) -o $(EXEC) $(LIBS)

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FOPTS) $<
