#!/bin/bash
year=2020
mnth=03
itdeb=20
itjum=1  # To jump 1 day (1: every 3 hours)
nit=167 
itfin=$(( itdeb + ( nit -1 ) * itjum ))
istri=1   # Jumps in latitude and longitude (stride)
version=QBOi
NSIDE=64
NPIX=49152
IIM=$(( 360 / istri ))
JJM=$(( 180 / istri + 1 )) 
# use the location where ICON data are stored
ICON=...
echo 'nit=' $nit 'itfin=' $itfin
echo 'nside='$NSIDE
echo 'iim='$IIM


declare -a nday=([1]=31 28 31 30 31 30 31 31 30 31 30 31)
let "reste = ${year} % 4"
if (($reste < 1))
then
   nday[2]=29
fi
echo $nday[2]

for mth in 1;do

itfin=$(( itdeb + ( nit -1 ) * itjum ))

\rm ../netsto/gwd_healpix_${version}_${year}${mnth}_${nit}it.* *.h *.f90 *.F  *.nc *.o *.mod input_data laun_gwd_healpix 

cat << EOF0 > dimensions.h

       INTEGER npix,llm, iim, jjm,llm_era5, ndm,nqmx
       INTEGER nside
       PARAMETER (iim=${IIM},jjm=${JJM},llm_era5=69,ndm=1,nqmx=2)
       PARAMETER (nside=${NSIDE},npix=${NPIX},llm=90)
EOF0
cat << EOF1 > input_data
$itdeb
$itfin
$itjum
$istri
EOF1
cat <<EOF2 > header_ctl
DSET       ../netsto/gwd_healpix_${version}_${year}${mnth}_${nit}it.nc
DTYPE      NETCDF                                              
TITLE      gwd_healpix_${version}_${year}${mnth}_${nit}it
EOF2


#Inputs:  #spring data

ln -s $ICON/nextgems_3dvariables_n64.nc VITUV_PRESS_TEMP.nc
ln -s $ICON/nextgems_precipflux_n64.nc PREC.nc
ln -s $ICON/nextgems_vvor_n64.nc VORT.nc
ln -s $ICON/products_winds_temp.nc products_winds.nc


ln -s ../prog/*.h .
ln -s ../prog/acama_gwd_rando_${version}.f90 acama_gwd_rando.f90
ln -s ../prog/flott_gwd_rando_${version}.f90 flott_gwd_rando.f90
ln -s ../prog/orografi_strato_${version}.f90 orografi_strato.f90
ln -s ../prog/laun_gwd_healpix.f90 laun_gwd_healpix.f90

make 

./laun_gwd_healpix < input_data


cp gwd_healpix_temp.nc ../netsto/gwd_healpix_${version}_${year}${mnth}_${nit}it.nc
cp gwd_healpix_temp_lonlat.nc ../netsto/gwd_healpix_lonlat_${version}_${year}${mnth}_${nit}it.nc



echo 'File:' gwd_healpix_${version}_${year}${mnth}
echo 'Iter:' ${nit}




done
