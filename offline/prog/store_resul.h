!  NETCDF STORAGE:      


!     Here storage

      status = nf90_inq_varid(ncio,"time",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             var1d(1),start=(/iter/))
!  2-D Variable:
      status = nf90_inq_varid(ncio,"prec",rvarid)
      var1d_hor_healp=prec
      status = nf90_put_var(ncio,rvarid,var1d_hor_healp,start=(/1,iter/), &
                            count=(/npix,1/))
      status = nf90_inq_varid(ncio,"zustrhi",rvarid)
      var1d_hor_healp=zustrhi
      status = nf90_put_var(ncio,rvarid,var1d_hor_healp,start=(/1,iter/), &
                            count=(/npix,1/))
      status = nf90_inq_varid(ncio,"zvstrhi",rvarid)
      var1d_hor_healp=zvstrhi
      status = nf90_put_var(ncio,rvarid,var1d_hor_healp,start=(/1,iter/), &
                            count=(/npix,1/))
      status = nf90_inq_varid(ncio,"zmea",rvarid)
      var1d_hor_healp=zmea
      status = nf90_put_var(ncio,rvarid,var1d_hor_healp,start=(/1/), &
                            count=(/npix/))
!  3-D Variables
      status = nf90_inq_varid(ncio,"vitu",rvarid)
      var2d_healp(:,:)=u(:,:)
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"vitv",rvarid)
      var2d_healp(:,:)=v(:,:)
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"temp",rvarid)
      var2d_healp(:,:)=t(:,:)
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                           count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"rot",rvarid)
      var2d_healp(:,:)=rot(:,:)
     status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"rho",rvarid)
      var2d_healp(:,:)=rho_z(:,:)
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))
 
      status = nf90_inq_varid(ncio,"rupwp",rvarid)
      var2d_healp(:,:)=rupwp(:,:)
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"d_u_mount",rvarid)
      var2d_healp(:,:)=d_u_oro(:,:) 
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"d_v_mount",rvarid)
      var2d_healp(:,:)=d_v_oro(:,:) 
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"d_u_preci",rvarid)
      var2d_healp(:,:)=d_u_lot(:,:)
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"d_u_front",rvarid)
      var2d_healp(:,:)=d_u_hin(:,:)
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"east_gwstress",rvarid)
       var2d_healp=east_acama+east_flott
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"west_gwstress",rvarid)
      var2d_healp=west_acama+west_flott
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"u_orstress",rvarid)
      var2d_healp=zustror_z
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"u_prstress",rvarid)
      var2d_healp=zustrpr_z
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"u_frstress",rvarid)
      var2d_healp=zustrfr_z
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))

      status = nf90_inq_varid(ncio,"u_stress_orpr",rvarid)
      var2d_healp=zustr_z
      status = nf90_put_var(ncio,rvarid,var2d_healp,start=(/1,1,iter/),&
                            count=(/npix,llm,1/))


