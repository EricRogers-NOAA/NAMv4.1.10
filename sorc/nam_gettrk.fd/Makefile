SHELL=	/bin/bash
LIBS=	${W3NCO_LIBd} ${W3EMC_LIBd} ${G2_LIBd} ${BACIO_LIB4} ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}
ISIZE = 4
RSIZE = 8
FCOMP=  ifort -nofree
CCOMP=  icc
INC = ${G2_INCd}
LDFLAGS= 
FFLAGS= -O2 -fpe0 -I $(INC)   -i$(ISIZE) -r$(RSIZE)
CFLAGS= -O2

gettrk:      gettrk_main.f gettrk_modules.o module_waitfor.o cwaitfor.o
	@echo " "
	@echo "  Compiling the main tracking program and subroutines....."
	$(FCOMP) $(FFLAGS) $(LDFLAGS) gettrk_modules.o module_waitfor.o cwaitfor.o gettrk_main.f $(LIBS) -o ${OUTDIR}/nam_gettrk
	@echo " "

cwaitfor.o: cwaitfor.c
	@echo " "
	@echo "  Compiling the waitfor C routine...."
	$(CCOMP) $(CFLAGS) -c cwaitfor.c -o cwaitfor.o

module_waitfor.o: module_waitfor.f
	@echo " "
	@echo "  Compiling the waitfor fortran module...."
	$(FCOMP) $(FFLAGS) -c module_waitfor.f -o module_waitfor.o

gettrk_modules.o:    gettrk_modules.f
	@echo " "
	@echo "  Compiling the regular tracker fortran modules....."
	$(FCOMP) $(FFLAGS) -c gettrk_modules.f -o gettrk_modules.o
	@echo " "

clean:
	-rm -f  *.o  *.mod
