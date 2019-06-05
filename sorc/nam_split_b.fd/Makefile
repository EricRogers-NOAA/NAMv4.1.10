SHELL=	/bin/bash
LIBS=	${W3NCO_LIBd} ${W3EMC_LIBd} ${SP_LIBd}
FC=	ifort
LDFLAGS=
FFLAGS	=-O3 -convert big_endian -DLINUX -fp-model precise -assume byterecl
FFLAGS_DOUBLE	 = -integer-size 32 -real-size 64
ESSLLIB        =  -mkl=sequential

SRCS =	modules.f wrf_split_new2.f BOUND.f fft99.f sig_p_convt1.f\
	SEPAR.f WNLIT.f FDUMP.f H12.f I1MACH.f J4SAVE.f XGETUA.f\
	WNLSM.f WNNLS.f XERABT.f XERCTL.f XERPRT.f XERROR.f XERRWV.f \
	XERSAV.f srotm.f srotmg.f amatrix.f rodist.f landcmsk.f\
        create_rel_domain.f

OBJS =	modules.o wrf_split_new2.o BOUND.o fft99.o sig_p_convt1.o\
	SEPAR.o WNLIT.o FDUMP.o H12.o I1MACH.o J4SAVE.o XGETUA.o\
	WNLSM.o WNNLS.o XERABT.o XERCTL.o XERPRT.o XERROR.o XERRWV.o \
	XERSAV.o srotm.o srotmg.o amatrix.o rodist.o landcmsk.o\
        create_rel_domain.o

CMD=	${OUTDIR}/nam_split_b.exe

all:$(CMD)

$(CMD):$(OBJS)
	$(FC) $(FFLAGS) $(FFLAGS_DOUBLE) -o $(@) $(OBJS) $(LIBS) $(ESSLLIB)

clean:
	-rm -rf *.o *.lst *.exe
