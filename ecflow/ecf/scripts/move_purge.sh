for file in `find . -name "*ecf" ` ; do
    if grep -q "module purge" $file ; then
        cp $file $file.bak
        sed  s/module\ purge//g $file.bak | sed  s/%include\ \<head.h\>/module\ purge\\n%include\ \<head.h\>/g  > $file
    fi
done

