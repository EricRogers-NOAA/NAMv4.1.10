# add modules just before calling the jjob;
# a few *ecf files don't call jjobs by "$NWROOT..." and must be edite manually;
# find them by:
#   for file in `find . -name "*ecf" ` ; do grep -q ^\$NWROOT $file ; [ $? -ne 0 ] && echo $file ; done 
#
for file in `find . -name "*ecf" ` ; do 
cp $file $file.orig
sed s/^\$NWROOT/\
\\n\
module\ purge\\n\
module\ load\ prod_util\\n\
module\ load\ grib_util\\/\${grib_util_ver:?}\\n\
module\ load\ ics\\/\${ics_ver:?}\\n\
module\ load\ NetCDF\\/4.2\\/serial\\n\
module\ load\ ibmpe\\n\
module\ load\ lsf\\/\${lsf_ver:?}\\n\
module\ load\ gempak\\/ncep\\n\
module\ load\ hpss\\/4.0.1.2\\n\
\\n\
\$NWROOT/g   $file.orig  > $file
rm $file.orig
done
exit
