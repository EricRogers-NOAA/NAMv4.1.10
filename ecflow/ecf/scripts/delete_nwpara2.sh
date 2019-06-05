for file in `find . -name "*ecf" ` ; do
cp $file $file.orig
sed s%module\ use\ \-a\ /nwpara2/modulefiles/%%g $file.orig  > $file
done

