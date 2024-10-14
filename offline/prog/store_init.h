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
!  Variables 2D with attributes:
      status = nf90_def_var(ncio,"prec",nf90_float, &
            (/healpixdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "mm/s")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Total precipitation")
      status = nf90_def_var(ncio,"zustrhi",nf90_float, &
            (/healpixdimid,timdimid/),healpixvarid)
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "GW surf. Zon. Stress (Pa)")
      nvar=nvar+1
      status = nf90_def_var(ncio,"zvstrhi",nf90_float, &
            (/healpixdimid,timdimid/),healpixvarid)
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "GW surf. Mer. Stress (Pa)")
      nvar=nvar+1
      status = nf90_def_var(ncio,"zmea",nf90_float, &
            (/healpixdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "zmean in healpix grids")

!  Variables 3D with attributes:
      status = nf90_def_var(ncio,"vitu",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "m/s")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Zonal wind")
      status = nf90_def_var(ncio,"vitv",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "m/s")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Meridional wind")
      status = nf90_def_var(ncio,"temp",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "Kelvin")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Temperature")
      status = nf90_def_var(ncio,"rot",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "s-1")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Rotational")
      status = nf90_def_var(ncio,"rho",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "rho(:,z)")
      status = nf90_def_var(ncio,"rupwp",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "uprime_wprime")
      status = nf90_def_var(ncio,"d_u_mount",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "GW and blocked flow drag due to mountains")
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "m/s/day")
      status = nf90_def_var(ncio,"d_v_mount",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "GW and blocked flow drag due to mountains")
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "m/s/day")
      status = nf90_def_var(ncio,"d_u_preci",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "m/s/day")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "GW drag due to precips")
      status = nf90_def_var(ncio,"d_u_front",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "m/s/day")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "GW drag due to fronts")
      status = nf90_def_var(ncio,"east_gwstress",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "East GW stress")
      status = nf90_def_var(ncio,"west_gwstress",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "West GW stress")
      status = nf90_def_var(ncio,"u_orstress",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Orographic stress as a function of z")
      status = nf90_def_var(ncio,"u_prstress",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Precipitation stress as a function of z")
      status = nf90_def_var(ncio,"u_frstress",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Front stress as a function of z")
      status = nf90_def_var(ncio,"u_stress_orpr",nf90_float, &
           (/healpixdimid,levdimid,timdimid/),healpixvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,healpixvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,healpixvarid,"Description", &
           "Total stress as a function of z")
      status = nf90_enddef(ncio)
!  End of data definition, entering data mode
      
      status = nf90_inq_varid(ncio,"level",rvarid)
      status = nf90_put_var(ncio,rvarid,zlev(:),start=(/1/), &
                            count=(/llm/))

!  Closing data
      status = nf90_close(ncio)
      print *,'erreur close?',status


       




       




