!  NETCDF STORAGE:      


!     Here storage
      status = nf90_inq_varid(ncio,"level",rvarid)
      status = nf90_put_var(ncio,rvarid,zlev(:),start=(/1/), &
                            count=(/llm/))

      status = nf90_inq_varid(ncio,"time",rvarid)
      status = nf90_put_var(ncio,rvarid, &
                             var1d(1),start=(/iter/))
!  2-D Variable:
      status = nf90_inq_varid(ncio,"prec",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var2d(ii,jj)=prec_lonlat(ij)
      enddo;enddo;
      status = nf90_put_var(ncio,rvarid,var2d,start=(/1,1,iter/), &
                            count=(/iim,jjm,1/))
      status = nf90_inq_varid(ncio,"rlatera5",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var2d(ii,jj)=rlat_era5(ij)
      enddo;enddo;
      status = nf90_put_var(ncio,rvarid,var2d,start=(/1,1,iter/), &
                            count=(/iim,jjm,1/))
      status = nf90_inq_varid(ncio,"rlat",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var2d(ii,jj)=rlat_lonlat(ij)
      enddo;enddo;
      status = nf90_put_var(ncio,rvarid,var2d,start=(/1,1,iter/), &
                            count=(/iim,jjm,1/))
      status = nf90_inq_varid(ncio,"zustrhi",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var2d(ii,jj)=zustrhi_lonlat(ij)
      enddo;enddo;
      status = nf90_put_var(ncio,rvarid,var2d,start=(/1,1,iter/),&
                            count=(/iim,jjm,1/))
      status = nf90_inq_varid(ncio,"zvstrhi",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var2d(ii,jj)=zvstrhi_lonlat(ij)
      enddo;enddo;
      status = nf90_put_var(ncio,rvarid,var2d,start=(/1,1,iter/),&
                            count=(/iim,jjm,1/))
!  3-D Variables
      status = nf90_inq_varid(ncio,"vitu",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=u_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"vitv",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=v_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"vitw",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=w_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"temp",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=t_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"pfull",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=pplay_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rho",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=rho_z_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rot",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=rot_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rupwp",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=rupwp_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"rvptp",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=rvptp_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"upwp",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=upwp_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"vptp",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=vptp_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"EPF",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=EPF_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"EPF_phi",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;var3d(ii,jj,:)=EPF_phi_lonlat(ij,:)
      enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"d_u_mount",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=d_u_oro_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"d_v_mount",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=d_v_oro_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"d_u_preci",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=d_u_lot_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"d_u_front",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=d_u_hin_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"east_gwstress",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=east_acama_lonlat(ij,:)+east_flott_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"west_gwstress",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=west_acama_lonlat(ij,:)+west_flott_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))


      status = nf90_inq_varid(ncio,"u_orstress",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=ustror_z_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"u_prstress",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=ustrpr_z_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"u_frstress",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=ustrfr_z_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"du_orstressdz",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=dustrordz_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"du_prstressdz",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=dustrprdz_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"du_frstressdz",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=dustrfrdz_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))


      status = nf90_inq_varid(ncio,"dEPFdz",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=dEPFdz_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"drupwpdz",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=drupwpdz_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))

      status = nf90_inq_varid(ncio,"dEPFdphi",rvarid)
      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim
      var3d(ii,jj,:)=dEPFdphi_lonlat(ij,:);enddo;enddo
      status = nf90_put_var(ncio,rvarid,var3d,start=(/1,1,1,iter/),&
                            count=(/iim,jjm,llm,1/))
