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
           "ICON level altitude")
      status = nf90_def_var(ncio, "time",nf90_float, &
           (/timdimid/),timvarid)
      status = nf90_put_att(ncio,timvarid,"units", &
           "seconds since 1970-01-01 00:00:0.0")
!  Variables 2D with attributes:
      status = nf90_def_var(ncio,"prec",nf90_float, &
            (/londimid,latdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/hr")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Total precipitation")
      status = nf90_def_var(ncio,"rlatera5",nf90_float, &
            (/londimid,latdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "radian")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "latitude")
      status = nf90_def_var(ncio,"rlat",nf90_float, &
            (/londimid,latdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "radian")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "latitude")
      status = nf90_def_var(ncio,"zustrhi",nf90_float, &
            (/londimid,latdimid,timdimid/),uvarid)
      status = nf90_put_att(ncio,uvarid,"Description", &
           "GW surf. Zon. Stress (Pa)")
      nvar=nvar+1
      status = nf90_def_var(ncio,"zvstrhi",nf90_float, &
            (/londimid,latdimid,timdimid/),uvarid)
      status = nf90_put_att(ncio,uvarid,"Description", &
           "GW surf. Mer. Stress (Pa)")
      nvar=nvar+1
!  Variables 3D with attributes:
      status = nf90_def_var(ncio,"vitu",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/s")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Zonal wind")
      status = nf90_def_var(ncio,"vitv",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/s")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Meridional wind")
      status = nf90_def_var(ncio,"vitw",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/s")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "w")
      status = nf90_def_var(ncio,"temp",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Kelvin")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Temperature")
      status = nf90_def_var(ncio,"rot",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "s-1")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Rotational")
      status = nf90_def_var(ncio,"pfull",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Air Pressure ")
      status = nf90_def_var(ncio,"rho",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "rho(:,z)")
      status = nf90_def_var(ncio,"rupwp",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "uprime_wprime")
      status = nf90_def_var(ncio,"rvptp",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "vprime_tprime")
      status = nf90_def_var(ncio,"upwp",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "uprime_wprime")
      status = nf90_def_var(ncio,"vptp",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "vprime_tprime")
      status = nf90_def_var(ncio,"EPF",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Eliassen-palm flux")
      status = nf90_def_var(ncio,"d_u_mount",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"Description", &
           "GW and blocked flow drag due to mountains")
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/s/day")
     status = nf90_def_var(ncio,"d_v_mount",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"Description", &
           "GW and blocked flow drag due to mountains")
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/s/day")
      status = nf90_def_var(ncio,"d_u_preci",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/s/day")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "GW drag due to precips")
      status = nf90_def_var(ncio,"d_u_front",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "m/s/day")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "GW drag due to fronts")
      status = nf90_def_var(ncio,"east_gwstress",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "East GW stress")
      status = nf90_def_var(ncio,"west_gwstress",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "West GW stress")
      status = nf90_def_var(ncio,"u_orstress",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Orographic stress as a function of z")
      status = nf90_def_var(ncio,"u_prstress",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Precipitation stress as a function of z")
      status = nf90_def_var(ncio,"u_frstress",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Front stress as a function of z")
      status = nf90_def_var(ncio,"u_stress_orpr",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Total stress(without front) as a function of z")
      status = nf90_def_var(ncio,"u_stress_orprfr",nf90_float, &
           (/londimid,latdimid,levdimid,timdimid/),uvarid)
      nvar=nvar+1
      status = nf90_put_att(ncio,uvarid,"units", &
           "Pa")
      status = nf90_put_att(ncio,uvarid,"Description", &
           "Total stress as a function of z")
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


       




       




