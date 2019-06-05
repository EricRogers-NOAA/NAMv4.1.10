for file in `find . -name "*ecf" ` ; do
cp $file $file.bak
sed  s%module\ load\ gempak/ncep%%g  $file.bak  > $file
done

