SHELL=  /bin/sh
ISIZE = 4
RSIZE = 8
COMP=	ifort
#LIBS=   -L/nwprod/lib -lw3emc_d -lw3nco_d
LIBS=	$(W3NCO_LIBd) $(W3EMC_LIBd)
LDFLAGS= 
FFLAGS= -O0 -g -check bounds -check format -i$(ISIZE) -r$(RSIZE)

supvit:     supvit_main.f supvit_modules.o
	@echo " "
	@echo "  Compiling program that sorts and updates vitals records...."
	$(COMP) $(FFLAGS) $(LDFLAGS) supvit_modules.o supvit_main.f $(LIBS) -o $(OUTDIR)/nam_supvit
	@echo " "

supvit_modules.o:   supvit_modules.f
	@echo " "
	@echo "  Compiling the modules....."
	$(COMP) -c supvit_modules.f -o supvit_modules.o
	@echo " "

clean:
	-rm -f  *.o  *.mod 

