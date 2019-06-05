#!/bin/ksh
#****************************************************************
#check if there is any hurricane vitals issued by NHC. If true,
#pick out number of storms, storm number, lats and lons
#usage: sh testTCvital.sh 2012080400
#Feb. 2014  - Guang Ping Lou
#****************************************************************
set -x

SDATE=$1

set -x
synfile2=${COMINnhc}/storm-data
ystart=`echo $SDATE | cut -c1-4`
mstart=`echo $SDATE | cut -c5-6`
dstart=`echo $SDATE | cut -c7-8`
hstart=`echo $SDATE | cut -c9-10`

EXEC=$EXECnam/nam_latlon_to_ij_h_bgrid_tc

rm -f tmpvit numb_tmpvit tmpvitall tmpvit1 tmpvit3
  grep "$ystart$mstart$dstart ${hstart}00" ${synfile} | \
      grep -v TEST | grep -v INVEST | awk 'substr($0,6,1) !~ /[8-9]/ {print $0}' | \
        grep NHC | grep -v ONE | grep -v TWO | grep -v THREE | grep -v FOUR | \
         grep -v FIVE | grep -v SIX | grep -v SEVEN | grep -v EIGHT | grep -v NINE | \
         grep -v TEN | grep -v ELEVEN | grep -v TWELVE | grep -v TEEN | \
         grep -v TWENTY | grep -v THIRTY >tmpvit3
    sort -k1,3 -u tmpvit3 > tmpvit1
lines=` cat tmpvit1 |wc -l`
if [ $lines -eq 0 ]; then
 echo No tropical storm ...  exit here 
  echo 0 > numb_tmpvit
   exit
 else
  >tmpvit2
  while read atcfrec1
   do
     stmnum=` echo $atcfrec1 | awk '{ print $2}' `
     grep -i $stmnum  tmpvit3 | tail -1 >> tmpvit2
   done < tmpvit1
fi

lines=` cat tmpvit2 |wc -l`
if [ $lines -eq 0 ]; then
 echo No tropical storm ...  exit here 
  echo 0 > numb_tmpvit
   exit
 else
  numbstm=1
  while read atcfrec
   do
     stmnum=` echo $atcfrec | awk '{ print $2}' `
     stmnam=` echo $atcfrec | awk '{ print $3}' `
     stmdev=` grep -h -i $stmnam ${synfile2}/*${ystart}/*${ystart}.com | \
               grep CARQ | grep "34, NEQ," | \
               awk 'substr($0,5,1) !~ /[8-9]/ {print $0}' | \
               awk '{ print $11}' | sed 's/,//' `
     if [ ${stmdev} = 'LO' -o ${stmdev} = 'DB' -o ${stmdev} = 'EX' -o ${stmdev} = 'WV' -o ${stmdev} = 'XX' ]; then
       echo "not a tropical storm..."
     else
     cenlat=` echo $atcfrec | awk '{ print $6}' | sed 's/[A-Z]//g' `
     cenlon=` echo $atcfrec | awk '{ print $7}' | sed 's/[A-Z]//g' `
     radius=` echo $atcfrec | awk '{ print $12}' `
     maxwnd=` echo $atcfrec | awk '{ print $13}' `
     cenlat=$((10#$cenlat))
     cenlat=$(echo "scale=1; $cenlat/10" |bc)
     cenlon=$((10#$cenlon))
     cenlon=$(echo "scale=1; $cenlon/10" |bc)
cat > itag <<EOF
$cenlat
-$cenlon
EOF
rm fort.51
ln -s inside fort.51
$EXEC < itag
INSIDE=`cut -c 1-1 inside`
if [ $INSIDE -eq 1 -a $maxwnd -ge 18 ] ; then
     echo $numbstm $stmnum $cenlat $cenlon $radius
     echo $numbstm > numb_tmpvit
     echo $stmnum $cenlat "-"$cenlon $radius > $GESDIR/${RUN}.t${cyc}z.latlons.${numbstm}.${tmmark}
     echo "$atcfrec" > $GESDIR/${RUN}.t${cyc}z.tmpvit.${numbstm}.${tmmark}
     cat $GESDIR/${RUN}.t${cyc}z.tmpvit.${numbstm}.${tmmark} >> tmpvitall
     let numbstm=$numbstm+1
fi
fi
done < tmpvit2

fi
  let numbstm=$numbstm-1
  echo "number of storms= " $numbstm

if [ $numbstm -gt 0 ] ; then
  cp tmpvitall $GESDIR/${RUN}.t${cyc}z.tmpvit.${tmmark}
  cp tmpvitall tmpvit
fi
