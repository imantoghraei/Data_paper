   ! Dimensions ID
       status = nf90_def_dim(ncio, "longitude",iim, londimid)
       status = nf90_def_dim(ncio, "latitude",jjm, latdimid)
       status = nf90_def_dim(ncio, "level",llm, levdimid)
       status = nf90_def_dim(ncio, "time",nf90_unlimited, timdimid)
!  Variables 1D with attributes
       status = nf90_def_var(ncio, "longitude",nf90_float, &
           (/londimid/),lonvarid)
      status = nf90_put_att(ncio,lonvarid,"units", &
           "degrees east")
       status = nf90_def_var(ncio, "latitude",nf90_float, &
           (/latdimid/),latvarid)
      status = nf90_put_att(ncio,latvarid,"units", &
           "degrees north")
       status = nf90_def_var(ncio, "level",nf90_float, &
           (/levdimid/),levvarid)
      status = nf90_put_att(ncio,levvarid,"units", &
           "km")
      status = nf90_put_att(ncio,levvarid,"Description", &
           "ICON geometric height at full level center")


       status = nf90_def_var(ncio, "rupwp_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustror_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "rupwp_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustror_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

!  Variables 2D with attributes:
 
    status = nf90_def_var(ncio,"rupwp_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
  
      status = nf90_def_var(ncio,"ustror_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrpr_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrfr_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1



  
      status = nf90_def_var(ncio, "rupwp_or_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1
  
       status = nf90_def_var(ncio, "ustror_or_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_or_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_or_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1
     status = nf90_def_var(ncio, "rupwp_or_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1
  
       status = nf90_def_var(ncio, "ustror_or_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_or_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_or_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1
 
     status = nf90_def_var(ncio,"rupwp_or_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
 
  
      status = nf90_def_var(ncio,"ustror_or_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrpr_or_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrfr_or_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1



    


      status = nf90_def_var(ncio, "rupwp_nor_land_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustror_nor_land_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_nor_land_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_nor_land_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1
     status = nf90_def_var(ncio, "rupwp_nor_land_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustror_nor_land_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_nor_land_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_nor_land_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

     status = nf90_def_var(ncio,"rupwp_nor_land_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustror_nor_land_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrpr_nor_land_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrfr_nor_land_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1




     status = nf90_def_var(ncio, "rupwp_nor_sea_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustror_nor_sea_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_nor_sea_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_nor_sea_bar_z",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1
     status = nf90_def_var(ncio, "rupwp_nor_sea_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustror_nor_sea_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrpr_nor_sea_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1

       status = nf90_def_var(ncio, "ustrfr_nor_sea_bar_z_w",nf90_float, &
           (/levdimid,timdimid/),uvarid)
      nvar=nvar+1


!  Variables 2D with attributes:
 

     status = nf90_def_var(ncio,"rupwp_nor_sea_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustror_nor_sea_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrpr_nor_sea_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrfr_nor_sea_zonal",nf90_float, &
            (/latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1



!  Variables 3D with attributes:
 

     status = nf90_def_var(ncio,"rupwp",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustror",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrpr",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrfr",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"rupwp-ustr",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

     status = nf90_def_var(ncio,"rupwp_or",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustror_or",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrpr_or",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrfr_or",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

     status = nf90_def_var(ncio,"rupwp_nor",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustror_nor",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrpr_nor",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_def_var(ncio,"ustrfr_nor",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1


    status = nf90_def_var(ncio,"rupwp_f2_1",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

    status = nf90_def_var(ncio,"rupwp_f2_2",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

    status = nf90_def_var(ncio,"rupwp_f2_3",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

    status = nf90_def_var(ncio,"rupwp_f2_4",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

    status = nf90_def_var(ncio,"rupwp_f2_5",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

    status = nf90_def_var(ncio,"rupwp_f2_6",nf90_float, &
            (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1

      status = nf90_enddef(ncio)

!  End of data definition, entering data mode

      status = nf90_inq_varid(ncio,"longitude",rvarid)
      status = nf90_put_var(ncio,rvarid,xlon(:),start=(/1/), &
                            count=(/iim/))

      status = nf90_inq_varid(ncio,"latitude",rvarid)
      status = nf90_put_var(ncio,rvarid,ylat(:),start=(/1/), &
                            count=(/jjm/))
      


!  Closing data
      status = nf90_close(ncio)
      print *,'erreur close?',status


       




       




