for file in `find . -name "*ecf" ` ; do 
cp $file $file.orig
awk '!/export\ MP/ && !/export\ OMP/ && !/export\ KMP/ {print}'  $file.orig > $file
rm $file.orig
done
exit

