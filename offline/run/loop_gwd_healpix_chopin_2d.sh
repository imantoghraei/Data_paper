#!/bin/bash
year=2020
mnth=12
itdeb=1
itjum=1  # To jump 1 day (1: every 3 hours)
nit=167
itfin=$(( itdeb + ( nit -1 ) * itjum ))
istri=1   # Jumps in latitude and longitude (stride)
version=QBOi
NSIDE=64
NPIX=49152
IIM=$(( 360 / istri ))
JJM=$(( 180 / istri + 1 )) 
#vesri=/Data/dsk1/flott/VESRI
#vesri=/sauvegarde/flott/dsk1/VESRI
vesri=/dsk1/flott/VESRI
echo 'nit=' $nit 'itfin=' $itfin
echo 'nside='$NSIDE
echo 'iim='$IIM

#declare -a nday=([1]=31 28 31 30 31 30 31 31 30 31 30 31)
declare -a nday=([1]=31 28 31 30 31 30 31 31 30 31 30 31)
let "reste = ${year} % 4"
if (($reste < 1))
then
   nday[2]=29
fi
echo $nday[2]

for mth in 1;do
# for mth in 1 2 3 4 5 6 7 8 9 10 11 12;do
#printf -v mnth %02d $mth
#let "nit = ${nday[$mth]}"
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
#ln -s $vesri/era5/data data   
#ln -s $vesri/laura/coarse_grained/nextgems_3dvariables_n64.nc VITUV_PRESS_TEMP.nc
#ln -s $vesri/laura/coarse_grained/nextgems_w_level_full_n64.nc VITW.nc
#ln -s $vesri/laura/coarse_june24/nextgems_precipflux_n64_spring.nc PREC.nc
#ln -s $vesri/laura/coarse_grained/vvor/nextgems_vvor_n64_spring.nc VORT.nc
#ln -s $vesri/laura/coarse_june24/products_winds_temp_20200320_20200409.nc products_winds.nc

#Inputs:  #summer data
#ln -s $vesri/laura/coarse_grained/nextgems_3dvariables_n64_summer.nc VITUV_PRESS_TEMP.nc
#ln -s $vesri/laura/coarse_grained/nextgems_precipflux_n64_summer.nc PREC.nc
#ln -s $vesri/laura/coarse_grained/vvor/nextgems_vvor_n64_summer.nc VORT.nc
#ln -s $vesri/laura/coarse_grained/products_winds_temp_20200620_20200710.nc products_winds.nc

#Inputs:  #autumn data
#ln -s $vesri/laura/coarse_grained/nextgems_3dvariables_n64_autumn.nc VITUV_PRESS_TEMP.nc
#ln -s $vesri/laura/coarse_june24/nextgems_precipflux_n64_autumn.nc PREC.nc
#ln -s $vesri/laura/coarse_grained/vvor/nextgems_vvor_n64_autumn.nc VORT.nc
#ln -s $vesri/laura/coarse_june24/products_winds_temp_20200920_20201010.nc products_winds.nc

#Inputs:  #winter data
ln -s $vesri/laura/coarse_grained/nextgems_3dvariables_n64_winter.nc VITUV_PRESS_TEMP.nc
ln -s $vesri/laura/coarse_june24/nextgems_precipflux_n64_winter.nc PREC.nc
ln -s $vesri/laura/coarse_grained/vvor/nextgems_vvor_n64_winter.nc VORT.nc
ln -s $vesri/laura/coarse_june24/products_winds_temp_20201220_20210109.nc products_winds.nc

ln -s ../prog/*.h .
ln -s ../prog/acama_gwd_rando_${version}.f90 acama_gwd_rando.f90
ln -s ../prog/flott_gwd_rando_${version}.f90 flott_gwd_rando.f90
ln -s ../prog/orografi_strato_${version}.f90 orografi_strato.f90
ln -s ../prog/laun_gwd_healpix.f90 laun_gwd_healpix.f90

make 

./laun_gwd_healpix < input_data


#echo sed -i 's/'${IIM}'/1/g' gwd_era5_temp.ctl 
#sed -i 's/'${IIM}'/1/g' gwd_era5_temp.ctl 
#cat header_ctl gwd_era5_temp.ctl > ../netsto/gwd_era5_${version}_${year}${mnth}_${nit}it.ctl
#ncwa -a longitude ../netsto/gwd_era5_temp.nc ../netsto/gwd_era5_${version}_${year}${mnth}_${nit}it.nc

cp gwd_healpix_temp.nc ../netsto/gwd_healpix_${version}_${year}${mnth}_${nit}it.nc
cp gwd_healpix_temp_lonlat.nc ../netsto/gwd_healpix_lonlat_${version}_${year}${mnth}_${nit}it.nc
cp gwd_zonal_lonlat.nc ../netsto/gwd_zonal_lonlat_${version}_${year}${mnth}_${nit}it.nc
cp gwd_zonal_lonlat_p.nc ../netsto/gwd_zonal_lonlat_p_${version}_${year}${mnth}_${nit}it.nc
cp gwd_zonal_lonlat_n.nc ../netsto/gwd_zonal_lonlat_n_${version}_${year}${mnth}_${nit}it.nc


#cp gwd_EPF_lonlat.nc ../netsto/gwd_EPF_lonlat${version}_${year}${mnth}_${nit}it.nc
#cdo zonmean gwd_healpix_temp_lonlat.nc ../netsto/gwd_healpix_lonlat_zonmean_${version}_${year}${mnth}_${nit}it.nc
#\rm ../netsto/gwd_era5_temp.nc 
echo 'File:' gwd_healpix_${version}_${year}${mnth}
echo 'Iter:' ${nit}


#\rm  *.h *.f90  *.nc *.o *.mod input_data laun_gwd_era5 data *ctl ../netsto/gwd_era5_temp.nc 

done
