for file in `find . -name "*ecf" ` ; do
cp $file $file.bak
sed  s/module\ load\ prod_util//g  $file.bak | sed s%module\ load\ hpss/4.0.1.2%%g > $file
done

