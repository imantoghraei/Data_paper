
!  NETCDF STORAGE:      


!     Here storage
      status = nf90_inq_varid(ncio,"level",rvarid)
      status = nf90_put_var(ncio,rvarid,zlev(:),start=(/1/), &
                            count=(/llm/))

     status = nf90_inq_varid(ncio,"rupwp_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustror_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))
     status = nf90_inq_varid(ncio,"rupwp_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustror_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))


!  2-D Variable:

      status = nf90_inq_varid(ncio,"rupwp_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

     status = nf90_inq_varid(ncio,"ustror_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrpr_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrfr_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))





     status = nf90_inq_varid(ncio,"rupwp_or_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_or_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustror_or_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_or_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_or_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_or_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_or_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_or_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))
     status = nf90_inq_varid(ncio,"rupwp_or_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_or_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustror_or_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_or_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_or_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_or_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_or_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_or_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))


!  2-D Variable:

      status = nf90_inq_varid(ncio,"rupwp_or_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_or_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

 
      status = nf90_inq_varid(ncio,"ustror_or_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_or_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrpr_or_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_or_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrfr_or_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_or_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

     

     status = nf90_inq_varid(ncio,"rupwp_nor_land_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_nor_land_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

 
    status = nf90_inq_varid(ncio,"ustror_nor_land_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_nor_land_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_nor_land_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_nor_land_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_nor_land_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_nor_land_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

     status = nf90_inq_varid(ncio,"rupwp_nor_land_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_nor_land_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

 
    status = nf90_inq_varid(ncio,"ustror_nor_land_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_nor_land_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_nor_land_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_nor_land_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_nor_land_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_nor_land_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))



!  2-D Variable:

      status = nf90_inq_varid(ncio,"rupwp_nor_land_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_nor_land_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))
 
      status = nf90_inq_varid(ncio,"ustror_nor_land_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_nor_land_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrpr_nor_land_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_nor_land_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrfr_nor_land_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_nor_land_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      



     status = nf90_inq_varid(ncio,"rupwp_nor_sea_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_nor_sea_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))


    status = nf90_inq_varid(ncio,"ustror_nor_sea_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_nor_sea_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_nor_sea_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_nor_sea_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_nor_sea_bar_z",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_nor_sea_bar_z,start=(/1,iter/), &
                            count=(/llm,1/))
     status = nf90_inq_varid(ncio,"rupwp_nor_sea_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_nor_sea_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))


    status = nf90_inq_varid(ncio,"ustror_nor_sea_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_nor_sea_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrpr_nor_sea_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_nor_sea_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))

    status = nf90_inq_varid(ncio,"ustrfr_nor_sea_bar_z_w",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_nor_sea_bar_z_w,start=(/1,iter/), &
                            count=(/llm,1/))


!  2-D Variable:

      status = nf90_inq_varid(ncio,"rupwp_nor_sea_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_nor_sea_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

  
      status = nf90_inq_varid(ncio,"ustror_nor_sea_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_nor_sea_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrpr_nor_sea_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_nor_sea_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrfr_nor_sea_zonal",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_nor_sea_zonal,start=(/1,1,iter/), &
                            count=(/jjm,llm,1/))

!  3-D Variables

      status = nf90_inq_varid(ncio,"rupwp",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

  
      status = nf90_inq_varid(ncio,"ustror",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrpr",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrfr",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp-ustr",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp-ustror-ustrpr-ustrfr,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp_or",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_or,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

  
      status = nf90_inq_varid(ncio,"ustror_or",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_or,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrpr_or",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_or,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrfr_or",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_or,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp_nor",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_nor_land+rupwp_nor_sea,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

  
      status = nf90_inq_varid(ncio,"ustror_nor",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustror_nor_land+ustror_nor_sea,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrpr_nor",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrpr_nor_land+ustrpr_nor_sea,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"ustrfr_nor",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             ustrfr_nor_land+ustrfr_nor_sea,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))


      status = nf90_inq_varid(ncio,"rupwp_f2_1",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_f2_1,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp_f2_2",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_f2_2,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp_f2_3",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_f2_3,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp_f2_4",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_f2_4,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp_f2_5",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_f2_5,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp_f2_6",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             rupwp_f2_6,start=(/1,1,1,iter/), &
                            count=(/iim,jjm,llm,1/))

