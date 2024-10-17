!  NETCDF STORAGE:      


!     Here storage

      status = nf90_inq_varid(ncio,"time",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             var1d(1),start=(/iter/))



