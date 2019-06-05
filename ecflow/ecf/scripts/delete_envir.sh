for file in `find . -name "*ecf" ` ; do
cp $file $file.orig
sed s/export\ envir=%ENVIR%//g $file.orig  > $file
done

