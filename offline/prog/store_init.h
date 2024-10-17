! Dimensions ID
       status = nf90_def_dim(ncio, "Healpix_grids",npix, healpixdimid)
       status = nf90_def_dim(ncio, "level",llm, levdimid)
       status = nf90_def_dim(ncio, "time",nf90_unlimited, timdimid)
!  Variables 1D with attributes
    
       status = nf90_def_var(ncio, "level",nf90_float, &
           (/levdimid/),levvarid)
      status = nf90_put_att(ncio,levvarid,"units", &
           "km?")
      status = nf90_put_att(ncio,levvarid,"Description", &
           "Log pressure altitude: -7*ln(bp+ap/101325)")
      status = nf90_def_var(ncio, "time",nf90_float, &
           (/timdimid/),timvarid)
      status = nf90_put_att(ncio,timvarid,"units", &
           "seconds since 1970-01-01 00:00:0.0")

!  End of data definition, entering data mode
      
      status = nf90_inq_varid(ncio,"level",rvarid)
      status = nf90_put_var(ncio,rvarid,zlev(:),start=(/1/), &
                            count=(/llm/))

!  Closing data
      status = nf90_close(ncio)
      print *,'erreur close?',status


       




       




