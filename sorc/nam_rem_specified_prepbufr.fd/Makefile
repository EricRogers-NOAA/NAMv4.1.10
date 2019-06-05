all:
	ifort -c -O2 -convert big_endian -list -auto nam_rem_specified_prepbufr.f
	ifort nam_rem_specified_prepbufr.o \
	       -o ${OUTDIR}/nam_rem_specified_prepbufr ${W3NCO_LIB4} ${W3EMC_LIB4} ${BUFR_LIB4}

clean:
	-rm -rf *.o *.lst

