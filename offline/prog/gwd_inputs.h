!!! lonlat fields
      REAL              :: zstd_cr
      REAL	        :: rupwp(iim,jjm,klev),ustror(iim,jjm,klev)
      REAL	        :: ustrfr(iim,jjm,klev),ustrpr(iim,jjm,klev)

      REAL              :: rupwp_or(iim,jjm,klev)  
      REAL              :: ustror_or(iim,jjm,klev)
      REAL              :: ustrpr_or(iim,jjm,klev)
      REAL              :: ustrfr_or(iim,jjm,klev)

      REAL              :: rupwp_nor_land(iim,jjm,klev)  
      REAL              :: ustror_nor_land(iim,jjm,klev)
      REAL              :: ustrpr_nor_land(iim,jjm,klev)
      REAL              :: ustrfr_nor_land(iim,jjm,klev)

      REAL              :: rupwp_nor_sea(iim,jjm,klev)  
      REAL              :: ustror_nor_sea(iim,jjm,klev)
      REAL              :: ustrpr_nor_sea(iim,jjm,klev)
      REAL              :: ustrfr_nor_sea(iim,jjm,klev)

      REAL              :: rupwp_bar_z(klev)  
      REAL              :: ustror_bar_z(klev)
      REAL              :: ustrpr_bar_z(klev)
      REAL              :: ustrfr_bar_z(klev)

      REAL              :: rupwp_bar_z_w(klev)  
      REAL              :: ustror_bar_z_w(klev)
      REAL              :: ustrpr_bar_z_w(klev)
      REAL              :: ustrfr_bar_z_w(klev)

      REAL              :: rupwp_zonal(jjm,klev)  
      REAL              :: ustror_zonal(jjm,klev)
      REAL              :: ustrpr_zonal(jjm,klev)
      REAL              :: ustrfr_zonal(jjm,klev)

      REAL              :: rupwp_or_bar_z(klev)  
      REAL              :: ustror_or_bar_z(klev)
      REAL              :: ustrpr_or_bar_z(klev)
      REAL              :: ustrfr_or_bar_z(klev)

      REAL              :: rupwp_or_bar_z_w(klev)  
      REAL              :: ustror_or_bar_z_w(klev)
      REAL              :: ustrpr_or_bar_z_w(klev)
      REAL              :: ustrfr_or_bar_z_w(klev)

      REAL              :: rupwp_or_zonal(jjm,klev)  
      REAL              :: ustror_or_zonal(jjm,klev)
      REAL              :: ustrpr_or_zonal(jjm,klev)
      REAL              :: ustrfr_or_zonal(jjm,klev)

      REAL              :: rupwp_nor_land_bar_z(klev)  
      REAL              :: ustror_nor_land_bar_z(klev)
      REAL              :: ustrpr_nor_land_bar_z(klev)
      REAL              :: ustrfr_nor_land_bar_z(klev)

      REAL              :: rupwp_nor_land_bar_z_w(klev)  
      REAL              :: ustror_nor_land_bar_z_w(klev)
      REAL              :: ustrpr_nor_land_bar_z_w(klev)
      REAL              :: ustrfr_nor_land_bar_z_w(klev)

      REAL              :: rupwp_nor_land_zonal(jjm,klev)  
      REAL              :: ustror_nor_land_zonal(jjm,klev)
      REAL              :: ustrpr_nor_land_zonal(jjm,klev)
      REAL              :: ustrfr_nor_land_zonal(jjm,klev)

      REAL              :: rupwp_nor_sea_bar_z(klev)  
      REAL              :: ustror_nor_sea_bar_z(klev)
      REAL              :: ustrpr_nor_sea_bar_z(klev)
      REAL              :: ustrfr_nor_sea_bar_z(klev)

      REAL              :: rupwp_nor_sea_bar_z_w(klev)  
      REAL              :: ustror_nor_sea_bar_z_w(klev)
      REAL              :: ustrpr_nor_sea_bar_z_w(klev)
      REAL              :: ustrfr_nor_sea_bar_z_w(klev)

      REAL              :: rupwp_nor_sea_zonal(jjm,klev)  
      REAL              :: ustror_nor_sea_zonal(jjm,klev)
      REAL              :: ustrpr_nor_sea_zonal(jjm,klev)
      REAL              :: ustrfr_nor_sea_zonal(jjm,klev)

      REAL              :: rupwp_mean(jjm,klev)  
      REAL              :: ustror_mean(jjm,klev)
      REAL              :: ustrpr_mean(jjm,klev)
      REAL              :: ustrfr_mean(jjm,klev)



      REAL              :: rupwp_f2_1(iim,jjm,klev)  
      REAL              :: rupwp_f2_2(iim,jjm,klev)
      REAL              :: rupwp_f2_3(iim,jjm,klev)
      REAL              :: rupwp_f2_4(iim,jjm,klev)
      REAL              :: rupwp_f2_5(iim,jjm,klev)
      REAL              :: rupwp_f2_6(iim,jjm,klev)

