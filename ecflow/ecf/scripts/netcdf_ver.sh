for file in `find . -name "*ecf" ` ; do
        cp $file $file.bak
        sed  s%4.2\/serial%\${netcdf_ver:?}\/serial%g $file.bak  > $file
done

