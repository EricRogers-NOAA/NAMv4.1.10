all:
	ifort -c -O -convert big_endian -list -auto latlon_to_ij_h_bgrid_tc.F90
	ifort latlon_to_ij_h_bgrid_tc.o \
	       -o ${OUTDIR}/nam_latlon_to_ij_h_bgrid_tc

clean:
	-rm -rf *.o *.lst

