#!/bin/bash

#BSUB -P FV3GFS-T2O
#BSUB -o namgsi.o%J
#BSUB -e namgsi.o%J
#BSUB -J namgsi
#BSUB -q dev
#BSUB -W 00:30
#BSUB -R span[ptile=7]
#BSUB -R affinity[core(4)]
#BSUB -x 
#BSUB -n 490

set -ax
date

export NODES=70
ntasks=490
ptile=7
threads=4

npe_node_max=28

module purge
module load EnvVars/1.0.2
module load lsf/10.1
module load ips/18.0.1.163
module load impi/18.0.1
module load prod_util/1.1.0
module load prod_envir/1.0.2
module load crtm/2.2.5
module load CFP/2.0.1
module load grib_util/1.1.0

module list


export OMP_NUM_THREADS=$threads
export KMP_STACKSIZE=2048M
export KMP_AFFINITY=scatter
export FORT_BUFFERED=true


# Set experiment name and analysis date
adate=2019011512
exp=test.$adate
TM=00

inputdir=/gpfs/dell2/emc/modeling/noscrub/emc.glopara/CASES/input_forgsitest 

# Set YES to use background from EMC paralel
# NO = take from /com/gfs/para
use_emc_para=NO


# Set YES to use ensemble, NO=standard 3dvar
DOHYBVAR=YES

# Set YES to run 4D-EnsVar.  NO=3D-EnsVar or 3DVAR
DO4DENSVAR=YES

# Set YES to use smoothed enkf forecasts
SMOOTH_ENKF=YES

# Set new radiance bias correction flag
export UNCOMPRESS=gunzip

# Generate diagnostic files
USE_RADSTAT=YES
GENDIAG=YES
DIAG_SUFFIX=""
CDATE=$adate
DIAG_COMPRESS=YES
COMPRESS=gzip
DIAG_TARBALL=YES
USE_CFP=YES





# Set path/file for gsi executable
gsiexec=/gpfs/dell2/emc/modeling/noscrub/emc.glopara/nam.v4.1.8/exec/nam_gsi
EXECnam=/gpfs/dell2/emc/modeling/noscrub/emc.glopara/nam.v4.1.8/exec

# Set the JCAP resolution which you want.
# All resolutions use LEVS=64
export JCAP=766
export JCAP_B=1534
export LEVS=64

# Set runtime and save directories
PTMP=/gpfs/dell2/ptmp
DATA=$PTMP/$LOGNAME/tmp${JCAP}/${exp}
SAVDIR=$PTMP/$LOGNAME/out${JCAP}/${exp}



# Specify GSI fixed field
fixgsi=/gpfs/dell2/emc/modeling/noscrub/emc.glopara/git/gsi/master/fix
FIXnam=/gpfs/dell2/emc/modeling/noscrub/emc.glopara/nam.v4.1.8/fix

FIXCRTM=$CRTM_FIX   # CRTM_FIX defined by crtm module


if [ $TM = 06 ] ; then
    v1=0.6
    fs=.false.
    berror=$FIXnam/nam_nmmstat_na_glberror.gcv
    anavinfo=$FIXnam/anavinfo_nems_nmmb_glb
else
    v1=1.0
    fs=.true.
    berror=$FIXnam/nam_nmmstat_na.gcv
    anavinfo=$FIXnam/anavinfo_nems_nmmb
fi

# Set GSD cloud analysis flag
i_gsdcldanal_type=2


# Set variables used in script
#   CLEAN up $DATA when finished (YES=remove, NO=leave alone)
#   ndate is a date manipulation utility
#   ncp is cp replacement, currently keep as /bin/cp

CLEAN=NO
NDATE=${NDATE:-/nwprod/util/exec/ndate}
export wc=${wc:-/usr/bin/wc}
ncpc=/bin/cp
ncpl="ln -fs"


# Set up $DATA
rm -rf $DATA
mkdir -p $DATA
cd $DATA
rm -rf core*

rm -f gsiparm.anl

## only use in master, not DA-FV3-IMPL
## binary_diag=$binary_diag,netcdf_diag=$netcdf_diag,
nhr_assimilation=03

cat <<EOF > gsiparm.anl
 &SETUP
   miter=2,niter(1)=50,niter(2)=50,niter_no_qc(1)=20,
   write_diag(1)=.true.,write_diag(2)=.false.,write_diag(3)=.true.,
   gencode=78,qoption=2,
   factqmin=0.0,factqmax=0.0,
   iguess=-1,use_gfs_ozone=.true.,
   oneobtest=.false.,retrieval=.false.,
   nhr_assimilation=$nhr_assimilation,l_foto=.false.,
   use_pbl=.false.,gpstop=30.,
   use_gfs_nemsio=.true.,
   print_diag_pcg=.true.,
   newpc4pred=.true., adp_anglebc=.true., angord=4,
   passive_bc=.true., use_edges=.false., emiss_bc=.true.,
   diag_precon=.true., step_start=1.e-3,
 /
 &GRIDOPTS
  wrf_nmm_regional=.false.,wrf_mass_regional=.false.,nems_nmmb_regional=.true.,diagnostic_reg=.false.,
  nmmb_reference_grid='H',grid_ratio_nmmb=3.0,
  filled_grid=.false.,half_grid=.false.,netcdf=.false.,nvege_type=20,
 /
 &BKGERR
   hzscl=0.373,0.746,1.50,
   vs=0.6,bw=0.,fstat=.true.,
 /
 &ANBKGERR
   anisotropic=.false.,
 /
 &JCOPTS
 /
 &STRONGOPTS
   nstrong=0,nvmodes_keep=20,period_max=3.,
    baldiag_full=.true.,baldiag_inc=.true.,  
 /
 &OBSQC
   dfact=0.75,dfact1=3.0,noiqc=.false.,c_varqc=0.02,
   vadfile='prepbufr',njqc=.false.,vqc=.true.,
 /
 &OBS_INPUT
   dmesh(1)=120.0,time_window_max=1.5,ext_sonde=.true.,
 /
OBS_INPUT::
!  dfile          dtype       dplat       dsis                  dval    dthin  dsfcalc
   prepbufr       ps          null        ps                  0.0     0     0
   prepbufr       t           null        t                   0.0     0     0
   prepbufr_profl t           null        t                   0.0     0     0
   prepbufr       q           null        q                   0.0     0     0
   prepbufr_profl q           null        q                   0.0     0     0
   prepbufr       pw          null        pw                  0.0     0     0
   prepbufr       uv          null        uv                  0.0     0     0
   prepbufr_profl uv          null        uv                  0.0     0     0
   satwndbufr     uv          null        uv                  0.0     0     0
   prepbufr       spd         null        spd                 0.0     0     0
   prepbufr       dw          null        dw                  0.0     0     0
   l2rwbufr       rw          null        rw                  0.0     0     0
   prepbufr       sst         null        sst                 0.0     0     0
   nsstbufr       sst         nsst        sst                 0.0     0     0
   gpsrobufr      gps_bnd     null        gps                 0.0     0     0
   hirs3bufr      hirs3       n17         hirs3_n17           0.0     1     0
   hirs4bufr      hirs4       metop-a     hirs4_metop-a       0.0     1     1
   gimgrbufr      goes_img    g11         imgr_g11            0.0     1     0
   gimgrbufr      goes_img    g12         imgr_g12            0.0     1     0
   airsbufr       airs        aqua        airs281SUBSET_aqua  0.0     1     1
   amsuabufr      amsua       n15         amsua_n15           0.0     1     1
   amsuabufr      amsua       n18         amsua_n18           0.0     1     1
   amsuabufr      amsua       metop-a     amsua_metop-a       0.0     1     1
   airsbufr       amsua       aqua        amsua_aqua          0.0     1     1
   amsubbufr      amsub       n17         amsub_n17           0.0     1     1
   mhsbufr        mhs         n18         mhs_n18             0.0     1     1
   mhsbufr        mhs         metop-a     mhs_metop-a         0.0     1     1
   ssmitbufr      ssmi        f14         ssmi_f14            0.0     1     0
   ssmitbufr      ssmi        f15         ssmi_f15            0.0     1     0
   amsrebufr      amsre_low   aqua        amsre_aqua          0.0     1     0
   amsrebufr      amsre_mid   aqua        amsre_aqua          0.0     1     0
   amsrebufr      amsre_hig   aqua        amsre_aqua          0.0     1     0
   ssmisbufr      ssmis       f16         ssmis_f16           0.0     1     0
   ssmisbufr      ssmis       f17         ssmis_f17           0.0     1     0
   ssmisbufr      ssmis       f18         ssmis_f18           0.0     1     0
   ssmisbufr      ssmis       f19         ssmis_f19           0.0     1     0
   gsnd1bufr      sndrd1      g12         sndrD1_g12          0.0     1     0
   gsnd1bufr      sndrd2      g12         sndrD2_g12          0.0     1     0
   gsnd1bufr      sndrd3      g12         sndrD3_g12          0.0     1     0
   gsnd1bufr      sndrd4      g12         sndrD4_g12          0.0     1     0
   gsnd1bufr      sndrd1      g11         sndrD1_g11          0.0     1     0
   gsnd1bufr      sndrd2      g11         sndrD2_g11          0.0     1     0
   gsnd1bufr      sndrd3      g11         sndrD3_g11          0.0     1     0
   gsnd1bufr      sndrd4      g11         sndrD4_g11          0.0     1     0
   gsnd1bufr      sndrd1      g13         sndrD1_g13          0.0     1     0
   gsnd1bufr      sndrd2      g13         sndrD2_g13          0.0     1     0
   gsnd1bufr      sndrd3      g13         sndrD3_g13          0.0     1     0
   gsnd1bufr      sndrd4      g13         sndrD4_g13          0.0     1     0
   iasibufr       iasi        metop-a     iasi616_metop-a     0.0     1     1
   omibufr        omi         aura        omi_aura            0.0     2     0
   hirs4bufr      hirs4       n19         hirs4_n19           0.0     1     1
   amsuabufr      amsua       n19         amsua_n19           0.0     1     1
   mhsbufr        mhs         n19         mhs_n19             0.0     1     1
   tcvitl         tcp         null        tcp                 0.0     0     0
   seviribufr     seviri      m08         seviri_m08          0.0     1     0
   seviribufr     seviri      m09         seviri_m09          0.0     1     0
   seviribufr     seviri      m10         seviri_m10          0.0     1     0
   hirs4bufr      hirs4       metop-b     hirs4_metop-b       0.0     1     1
   amsuabufr      amsua       metop-b     amsua_metop-b       0.0     1     1
   mhsbufr        mhs         metop-b     mhs_metop-b         0.0     1     1
   iasibufr       iasi        metop-b     iasi616_metop-b     0.0     1     1
   atmsbufr       atms        npp         atms_npp            0.0     1     0
   crisbufr       cris        npp         cris_npp            0.0     1     0
   gsnd1bufr      sndrd1      g14         sndrD1_g14          0.0     1     0
   gsnd1bufr      sndrd2      g14         sndrD2_g14          0.0     1     0
   gsnd1bufr      sndrd3      g14         sndrD3_g14          0.0     1     0
   gsnd1bufr      sndrd4      g14         sndrD4_g14          0.0     1     0
   gsnd1bufr      sndrd1      g15         sndrD1_g15          0.0     1     0
   gsnd1bufr      sndrd2      g15         sndrD2_g15          0.0     1     0
   gsnd1bufr      sndrd3      g15         sndrD3_g15          0.0     1     0
   gsnd1bufr      sndrd4      g15         sndrD4_g15          0.0     1     0
   oscatbufr      uv          null        uv                  0.0     0     0
   mlsbufr        mls30       aura        mls30_aura          0.0     0     0
   avhambufr      avhrr       metop-a     avhrr3_metop-a      0.0     1     0
   avhpmbufr      avhrr       n18         avhrr3_n18          0.0     1     0
   prepbufr       mta_cld     null        mta_cld             1.0     0     0     
   prepbufr       gos_ctp     null        gos_ctp             1.0     0     0     
   lgycldbufr     larccld     null        larccld             1.0     0     0
   lghtnbufr      lghtn       null        lghtn               1.0     0     0
::
 &SUPEROB_RADAR
   del_azimuth=5.,del_elev=.25,del_range=5000.,del_time=.5,elev_angle_max=5.,minnum=50,range_max=100000.,
   l2superob_only=.false.,
 /
 &LAG_DATA
 /
 &HYBRID_ENSEMBLE
   l_hyb_ens=.true.,
   n_ens=81,
   uv_hyb_ens=.true.,
   beta1_inv=0.25,
   s_ens_h=300,
   s_ens_v=5,
   generate_ens=.false.,
   regional_ensemble_option=1,
   aniso_a_en=.false.,
   nlon_ens=1536,
   nlat_ens=770,
   jcap_ens=766,
   l_ens_in_diff_time=.true.,
   jcap_ens_test=0,coef_bw=0.5,
   full_ensemble=.true.,betaflg=.true.,pwgtflg=.true.,
   ensemble_path="",
 /
 &RAPIDREFRESH_CLDSURF
   i_gsdcldanal_type=$i_gsdcldanal_type,
   dfi_radar_latent_heat_time_period=20.0,
   l_use_hydroretrieval_all=.false.,
   metar_impact_radius=10.0,
   metar_impact_radius_lowCloud=4.0,
   l_gsd_terrain_match_surfTobs=.false.,
   l_sfcobserror_ramp_t=.false.,
   l_sfcobserror_ramp_q=.false.,
   l_PBL_pseudo_SurfobsT=.false.,
   l_PBL_pseudo_SurfobsQ=.false.,
   l_PBL_pseudo_SurfobsUV=.false.,
   pblH_ration=0.75,
   pps_press_incr=20.0,
   l_gsd_limit_ocean_q=.false.,
   l_pw_hgt_adjust=.false.,
   l_limit_pw_innov=.false.,
   max_innov_pct=0.1,
   l_cleanSnow_WarmTs=.false.,
   r_cleanSnow_WarmTs_threshold=5.0,
   l_conserve_thetaV=.false.,
   i_conserve_thetaV_iternum=3,
   l_cld_bld=.false.,
   cld_bld_hgt=1200.0,
   build_cloud_frac_p=0.50,
   clear_cloud_frac_p=0.1,
   iclean_hydro_withRef=1,
   iclean_hydro_withRef_allcol=0,
 /
 &CHEM
 /
 &SINGLEOB_TEST
   maginnov=0.1,magoberr=0.1,oneob_type='t',
   oblat=45.,oblon=270.,obpres=850.,obdattim=${adate},
   obhourset=0.,
 /
EOF


# Set fixed files
#   berror   = forecast model background error statistics
#   specoef  = CRTM spectral coefficients
#   trncoef  = CRTM transmittance coefficients
#   emiscoef = CRTM coefficients for IR sea surface emissivity model
#   aerocoef = CRTM coefficients for aerosol effects
#   cldcoef  = CRTM coefficients for cloud effects
#   satinfo  = text file with information about assimilation of brightness temperatures
#   satangl  = angle dependent bias correction file (fixed in time)
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (optional)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)
#   aeroinfo = text file with information about assimilation of aerosol data

emiscoef_IRwater=$FIXCRTM/Nalli.IRwater.EmisCoeff.bin
emiscoef_IRice=$FIXCRTM/NPOESS.IRice.EmisCoeff.bin
emiscoef_IRland=$FIXCRTM/NPOESS.IRland.EmisCoeff.bin
emiscoef_IRsnow=$FIXCRTM/NPOESS.IRsnow.EmisCoeff.bin
emiscoef_VISice=$FIXCRTM/NPOESS.VISice.EmisCoeff.bin
emiscoef_VISland=$FIXCRTM/NPOESS.VISland.EmisCoeff.bin
emiscoef_VISsnow=$FIXCRTM/NPOESS.VISsnow.EmisCoeff.bin
emiscoef_VISwater=$FIXCRTM/NPOESS.VISwater.EmisCoeff.bin
emiscoef_MWwater=$FIXCRTM/FASTEM6.MWwater.EmisCoeff.bin
aercoef=$FIXCRTM/AerosolCoeff.bin
cldcoef=$FIXCRTM/CloudCoeff.bin
satinfo=$FIXnam/nam_regional_satinfo.txt
scaninfo=$FIXnam/nam_regional_scaninfo.txt
pcpinfo=$FIXnam/nam_global_pcpinfo.txt
ozinfo=$FIXnam/nam_global_ozinfo.txt
errtable=$FIXnam/nam_errtable.r3dv
convinfo=$FIXnam/nam_regional_convinfo.txt
mesonetuselist=$FIXnam/nam_mesonet_uselist.txt
stnuselist=$FIXnam/nam_mesonet_stnuselist.txt
qdaylist=$FIXnam/q_day_rejectlist
qnightlist=$FIXnam/q_night_rejectlist
tdaylist=$FIXnam/t_day_rejectlist
unightlist=$FIXnam/t_night_rejectlist
wbinuselist=$FIXnam/wbinuselist
atmsfilter=$fixgsi/atms_beamwidth.txt


# Copy observational data to $DATA
$ncpc $inputdir/* .


# Copy executable and fixed files to $DATA
$ncpc $gsiexec ./gsi.x
$ncpc $anavinfo ./anavinfo
$ncpc $berror   ./berror_stats
$ncpc $errtable ./errtable
$ncpc $emiscoef_IRwater ./Nalli.IRwater.EmisCoeff.bin
$ncpc $emiscoef_IRice ./NPOESS.IRice.EmisCoeff.bin
$ncpc $emiscoef_IRsnow ./NPOESS.IRsnow.EmisCoeff.bin
$ncpc $emiscoef_IRland ./NPOESS.IRland.EmisCoeff.bin
$ncpc $emiscoef_VISice ./NPOESS.VISice.EmisCoeff.bin
$ncpc $emiscoef_VISland ./NPOESS.VISland.EmisCoeff.bin
$ncpc $emiscoef_VISsnow ./NPOESS.VISsnow.EmisCoeff.bin
$ncpc $emiscoef_VISwater ./NPOESS.VISwater.EmisCoeff.bin
$ncpc $emiscoef_MWwater ./FASTEM6.MWwater.EmisCoeff.bin
$ncpc $aercoef  ./AerosolCoeff.bin
$ncpc $cldcoef  ./CloudCoeff.bin
$ncpc $satinfo  ./satinfo
$ncpc $scaninfo ./scaninfo
$ncpc $pcpinfo  ./pcpinfo
$ncpc $ozinfo   ./ozinfo
$ncpc $convinfo ./convinfo
$ncpc $mesonetuselist ./mesonetuselist
$ncpc $stnuselist ./mesonet_stnuselist
$ncpc $qdaylist ./q_day_rejectlist
$ncpc $qnightlist ./q_night_rejectlist
$ncpc $tdaylist ./t_day_rejectlist
$ncpc $tnightlist ./t_night_rejectlist
$ncpc $wbinuselist ./wbinuselist

mydbzinfile=refd3d.t12z.grb2f00

# Copy CRTM coefficient files based on entries in satinfo file

set +x
for file in `awk '{if($1!~"!"){print $1}}' satinfo | sort | uniq` ;do
   $ncpc $FIXCRTM/${file}.SpcCoeff.bin ./
   $ncpc $FIXCRTM/${file}.TauCoeff.bin ./
done
set -x

# Copy observational data to $DATA
##$ncpc $inputdir/* .


# Copy bias correction, atmospheric and surface files


if [ $i_gsdcldanal_type -eq 2 ] ; then
  echo "RUNNING CLOUD ANALYSIS PREP `date`"
  rm -f anavinfo
  $ncpc $FIXnam/anavinfo_nems_nmmb_cld  ./anavinfo
  outfile=ref3d
  $WGRIB2 -s $mydbzinfile | grep ":REFD:" | $WGRIB2 -i $mydbzinfile -text $outfile
  $EXECnam/nam_ref2nemsio wrf_inout${nhr_assimilation} > ref_log >> ref_stderr
  rc==$?
  echo "CLOUD ANALYSIS PREP DONE `date`"
fi


# Run gsi under Parallel Operating Environment (poe) on NCEP IBM

printenv > stdout.env

APRUN="mpirun -n $ntasks"
date
$APRUN $DATA/gsi.x < gsiparm.anl 1> stdout 2> stderr
rc=$?
date

exit

##cat fort.2* > gdas.t${cyca}z.gsistat

if [[ "$GENDIAG" = "NO" ]] ; then
  date
  exit
fi

# Save output
mkdir -p $SAVDIR
cat stdout fort.2* > $SAVDIR/stdout.anl.$adate
cat fort.2*        > $SAVDIR/${prefix_obs}.gsistat
cat fort.2*        > $SAVDIR/gsistat.$dumpobs.$adate
$ncpc siganl          $SAVDIR/gfnanl.$dumpobs.$adate
$ncpc satbias_out     $SAVDIR/biascr.$dumpobs.$adate
$ncpc satbias_pc.out  $SAVDIR/biascr_pc.$dumpobs.$adate
$ncpc satbias_out.int $SAVDIR/biascr.int.$dumpobs.$adate

CNVSTAT=$SAVDIR/cnvstat.gdas.$adate
PCPSTAT=$SAVDIR/pcpstat.gdas.$adate
OZNSTAT=$SAVDIR/oznstat.gdas.$adate
RADSTAT=$SAVDIR/radstat.gdas.$adate

rm -f $CNVSTAT
rm -f $PCPSTAT
rm -f $OZNSTAT
rm -f $RADSTAT


cd $DATA    # we should already be in $DATA, but extra cd to be sure.
rm -rf diag_*

echo "before GENDIAG= $GENDIAG at `date`"

# If requested, generate diagnostic files
if [ $GENDIAG = "YES" ] ; then

   # Set up lists and variables for various types of diagnostic files.
   ntype=3

   diagtype[0]="conv conv_gps conv_ps conv_q conv_sst conv_t conv_uv"
   diagtype[1]="pcp_ssmi_dmsp pcp_tmi_trmm"
   diagtype[2]="sbuv2_n16 sbuv2_n17 sbuv2_n18 sbuv2_n19 gome_metop-a gome_metop-b omi_aura mls30_aura ompsnp_npp ompstc8_npp"
   diagtype[3]="hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep sndr_g12_prep sndr_g13_prep sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13 sndrd3_g13 sndrd4_g13 sndrd1_g14 sndrd2_g14 sndrd3_g14 sndrd4_g14 sndrd1_g15 sndrd2_g15 sndrd3_g15 sndrd4_g15 hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua amsua_aqua imgr_g08 imgr_g11 imgr_g12 imgr_g14 imgr_g15 ssmi_f13 ssmi_f15 hirs4_n18 hirs4_metop-a amsua_n18 amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_f16 ssmis_f17 ssmis_f18 ssmis_f19 ssmis_f20 iasi_metop-a hirs4_n19 amsua_n19 mhs_n19 seviri_m08 seviri_m09 seviri_m10 seviri_m11 cris_npp cris-fsr_npp cris-fsr_n20 atms_npp atms_n20 hirs4_metop-b amsua_metop-b mhs_metop-b iasi_metop-b avhrr_n18 avhrr_metop-a amsr2_gcom-w1 gmi_gpm saphir_meghat ahi_himawari8"

   diaglist[0]=listcnv
   diaglist[1]=listpcp
   diaglist[2]=listozn
   diaglist[3]=listrad

   diagfile[0]=$CNVSTAT
   diagfile[1]=$PCPSTAT
   diagfile[2]=$OZNSTAT
   diagfile[3]=$RADSTAT

   numfile[0]=0
   numfile[1]=0
   numfile[2]=0
   numfile[3]=0

   # Set diagnostic file prefix based on lrun_subdirs variable
   if [ $lrun_subdirs = ".true." ]; then
      prefix=" dir.*/"
   else
      prefix="pe*"
   fi

   if [ $USE_CFP = "YES" ]; then
      rm $DATA/diag.sh $DATA/mp_diag.sh
      cat > $DATA/diag.sh << EOFdiag
#!/bin/sh
lrun_subdirs=\$1
binary_diag=\$2
type=\$3
loop=\$4
string=\$5
CDATE=\$6
DIAG_COMPRESS=\$7
DIAG_SUFFIX=\$8
if [ \$lrun_subdirs = ".true." ]; then
   prefix=" dir.*/"
else
   prefix="pe*"
fi
file=diag_\${type}_\${string}.\${CDATE}\${DIAG_SUFFIX}
if [ \$binary_diag = ".true." ]; then
   cat \${prefix}\${type}_\${loop}* > \$file
else
   $catexec -o \$file \${prefix}\${type}_\${loop}*
fi
if [ \$DIAG_COMPRESS = "YES" ]; then
   $COMPRESS \$file
fi
EOFdiag
      chmod 755 $DATA/diag.sh
   fi

   # Collect diagnostic files as a function of loop and type.
   # Loop over first and last outer loops to generate innovation
   # diagnostic files for indicated observation types (groups)
   #
   # NOTE:  Since we set miter=2 in GSI namelist SETUP, outer
   #        loop 03 will contain innovations with respect to
   #        the analysis.  Creation of o-a innovation files
   #        is triggered by write_diag(3)=.true.  The setting
   #        write_diag(1)=.true. turns on creation of o-g
   #        innovation files.

   loops="01 03"
   for loop in $loops; do
      case $loop in
         01) string=ges;;
         03) string=anl;;
          *) string=$loop;;
      esac
      echo $(date) START loop $string >&2
      n=-1
      while [ $((n+=1)) -le $ntype ] ;do
         for type in $(echo ${diagtype[n]}); do
            count=$(ls ${prefix}${type}_${loop}* | wc -l)
            if [ $count -gt 0 ]; then
               if [ $USE_CFP = "YES" ]; then
                  echo "$DATA/diag.sh $lrun_subdirs $binary_diag $type $loop $string $CDATE $DIAG_COMPRESS $DIAG_SUFFIX" | tee -a $DATA/mp_diag.sh
               else
                  cat ${prefix}${type}_${loop}* > diag_${type}_${string}.${CDATE}${DIAG_SUFFIX}
               fi
               echo "diag_${type}_${string}.${CDATE}*" >> ${diaglist[n]}
               numfile[n]=$(expr ${numfile[n]} + 1)
            fi
         done
      done
      echo $(date) END loop $string >&2
   done

   # We should already be in $DATA, but extra cd to be sure.
   cd $DATA

   # If requested, compress diagnostic files
   if [ $DIAG_COMPRESS = "YES" -a $USE_CFP = "NO" ]; then
      echo $(date) START $COMPRESS diagnostic files >&2
      for file in $(ls diag_*${CDATE}${DIAG_SUFFIX}); do
         $COMPRESS $file
      done
      echo $(date) END $COMPRESS diagnostic files >&2
   fi

   if [ $USE_CFP = "YES" ] ; then
      chmod 755 $DATA/mp_diag.sh
      ncmd=$(cat $DATA/mp_diag.sh | wc -l)
      if [ $ncmd -gt 0 ]; then
         ncmd_max=$((ncmd < npe_node_max ? ncmd : npe_node_max))
         APRUNCFP_DIAG=$(eval echo $APRUNCFP)
         $APRUNCFP_DIAG $DATA/mp_diag.sh
      fi
   fi

   # If requested, create diagnostic file tarballs
   if [ $DIAG_TARBALL = "YES" ]; then
      echo $(date) START tar diagnostic files >&2
      n=-1
      while [ $((n+=1)) -le $ntype ] ;do
         TAROPTS="-uvf"
         if [ ! -s ${diagfile[n]} ]; then
            TAROPTS="-cvf"
         fi
         if [ ${numfile[n]} -gt 0 ]; then
            tar $TAROPTS ${diagfile[n]} $(cat ${diaglist[n]})
         fi
      done

      # Restrict CNVSTAT
      chmod 750 $CNVSTAT
      ${CHGRP_CMD} $CNVSTAT

      # Restrict RADSTAT
      chmod 750 $RADSTAT
      ${CHGRP_CMD} $RADSTAT

      echo $(date) END tar diagnostic files >&2
   fi

fi # End diagnostic file generation block - if [ $GENDIAG = "YES" ]

echo "after GENDIAG= $GENDIAG at `date`"

# If requested, clean up $DATA
if [[ "$CLEAN" = "YES" ]];then
   if [[ $rc -eq 0 ]];then
      rm -rf $DATA
      cd $DATA
      cd ../
      rmdir $DATA
   fi
fi
date

exit
