      PROGRAM LAUN_SSO 

      use netcdf
      IMPLICIT none
!NETCDF STUFF
      integer ncid,status,varid,dimid
      integer start(4),count(4),stride(4),ilmdz
!NETCDF STUFF FOR OUTPUT FILE:
      integer ncio
      integer londimid,latdimid,levdimid,timdimid
      integer lonvarid,latvarid,levvarid,timvarid
      integer uvarid,vvarid,tvarid,rvarid
!FORMATTED I/O STUFF
      character*1 entete
!LOOPS ON TIMES
      integer it,nit,nvar,istri



include "dimensions.h" ! IT: New dimension.h for the healpix grids
include "dimphy.h"   ! IT: New dimphy.h for the healpix grids
include "YOEGWD.h"
      
! Those thinks are on the dynamical grid?

       integer ii,jj,ij
      real, allocatable :: var1d(:),var2d(:,:),var3d(:,:,:) 
      integer i,j,l,ll,index




!  IMPUT DYNAMICAL FIELDS THAT SHOULD BE PROVIDED OBVIOUSLY
!  PLACED ON THE PHYSICAL GRID:

      integer zlev(klev)
    


! TIME STEP OF THE PHYSICS, time of the data
!
!      REAL dtime
!      INTEGER itime
!

!==============================IT======================================
!!! lonlat fields

      
      REAL, allocatable :: upwp_lonlat(:,:,:,:)  ! IT: uprime*wprime from Laura's data
      REAL, allocatable :: ustror_z_lonlat(:,:,:,:)
      REAL, allocatable :: ustrpr_z_lonlat(:,:,:,:)
      REAL, allocatable :: ustrfr_z_lonlat(:,:,:,:)
      REAL, allocatable :: ustrorpr_z_lonlat(:,:,:,:)
      REAL, allocatable :: ustrorprfr_z_lonlat(:,:,:,:)

      REAL, allocatable :: upwp_bar_zt(:,:)  ! IT: uprime*wprime from Laura's data
      REAL, allocatable :: ustror_bar_zt(:,:)
      REAL, allocatable :: ustrpr_bar_zt(:,:)
      REAL, allocatable :: ustrfr_bar_zt(:,:)
      REAL, allocatable :: ustrorpr_bar_zt(:,:)
      REAL, allocatable :: ustrorprfr_bar_zt(:,:)

      REAL, allocatable :: upwp_bar_z(:)  ! IT: uprime*wprime from Laura's data
      REAL, allocatable :: ustror_bar_z(:)
      REAL, allocatable :: ustrpr_bar_z(:)
      REAL, allocatable :: ustrfr_bar_z(:)
      REAL, allocatable :: ustrorpr_bar_z(:)
      REAL, allocatable :: ustrorprfr_bar_z(:)

      REAL, allocatable :: upwp_zonal_t(:,:,:)  ! IT: uprime*wprime from Laura's data
      REAL, allocatable :: ustror_zonal_t(:,:,:)
      REAL, allocatable :: ustrpr_zonal_t(:,:,:)
      REAL, allocatable :: ustrfr_zonal_t(:,:,:)
      REAL, allocatable :: ustrorpr_zonal_t(:,:,:)
      REAL, allocatable :: ustrorfr_zonal_t(:,:,:)

      REAL, allocatable :: upwp_zonal(:,:)  ! IT: uprime*wprime from Laura's data
      REAL, allocatable :: ustror_zonal(:,:)
      REAL, allocatable :: ustrpr_zonal(:,:)
      REAL, allocatable :: ustrfr_zonal(:,:)
      REAL, allocatable :: ustrorpr_zonal(:,:)
      REAL, allocatable :: ustrorprfr_zonal(:,:)


!!!! For calculating the log pressure vert. coordinate
    REAL ZH(klon,klev+1)
    REAL PSEC ! Security to avoid division by 0 pressure
    REAL H0 ! Characteristic Height of the atmosphere
    REAL PR , TR ! Reference Pressure and Temperature

!======================================================================

! OFF-DEB
!   THIS IS THE TIME-STEP THE PHYSICS SEE, IT IS
!  PUT TO 30 MINUTES HERE RATHER ABITRARILY

      DTIME=30.*60.
      DTIME=3.*3600.        !!! Healpix data are every 3 hours

!  Load the number of dates you want to proceed:

      read(*,*)nit
      read(*,*)istri



!  READING THE FIELDS FROM THE STORED DATA:

!  LOADING STRESSES
         
        status = NF90_OPEN( &
      'gwd_healpix_temp_lonlat.nc',NF90_NOWRITE,ncid)
      status = NF90_INQ_VARID(ncid,'longitude',varid)
      status = NF90_INQ_DIMID(ncid,'longitude',dimid)
      status = NF90_INQUIRE_DIMENSION(ncid,dimid,len=ii)
      if(ii.ne.iim) then;print *,'x-Dim wrong?';
      status = &
        NF90_GET_VAR(ncid, varid, xlon,start=(/1/),count=(/iim/), &
                     stride=(/istri/))
      status = NF90_INQ_VARID(ncid,'latitude',varid)
      status = NF90_INQ_DIMID(ncid,'latitude',dimid)
      status = NF90_INQUIRE_DIMENSION(ncid,dimid,len=jj)
      if(jj.ne.jjm) then;print *,'y-Dim wrong?'
      status = &
        NF90_GET_VAR(ncid, varid, ylat,start=(/1/),count=(/jjm/), &
                     stride=(/istri/))

      status = NF90_INQ_VARID(ncid,'level',varid)
      status = &
        NF90_GET_VAR(ncid, varid, zlev,start=(/1/),count=(/llm/), &
                     stride=(/istri/))

      status = NF90_INQ_VARID(ncid,'time',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var1d,start=(/1/),count=(/nit/), &
                     stride=(/istri/))

      
      allocate(var2d(iim,jjm)) 
      allocate(var3d(iim,jjm,llm)) 
      allocate(upwp_lonlat(iim,jjm,llm,nit)) 
      allocate(ustror_z_lonlat(iim,jjm,llm,nit)) 
      allocate(ustrpr_z_lonlat(iim,jjm,llm,nit)) 
      allocate(ustrfr_z_lonlat(iim,jjm,llm,nit)) 
      allocate(ustrorpr_z_lonlat(iim,jjm,llm,nit)) 
      allocate(ustrorprfr_z_lonlat(iim,jjm,llm,nit)) 

      allocate(upwp_bar_zt(llm,nit)) 
      allocate(ustror_bar_zt(llm,nit)) 
      allocate(ustrpr_bar_zt(llm,nit)) 
      allocate(ustrfr_bar_zt(llm,nit)) 
      allocate(ustrorpr_bar_zt(llm,nit))
      allocate(ustrorprfr_bar_zt(llm,nit))

      allocate(upwp_bar_z(llm)) 
      allocate(ustror_bar_z(llm)) 
      allocate(ustrpr_bar_z(llm)) 
      allocate(ustrfr_bar_z(llm)) 
      allocate(ustrorpr_bar_z(llm)) 
      allocate(ustrorprfr_bar_z(llm)) 

      allocate(upwp_zonal_t(iim,llm,nit)) 
      allocate(ustror_zonal_t(iim,llm,nit)) 
      allocate(ustrpr_zonal_t(iim,llm,nit)) 
      allocate(ustrfr_zonal_t(iim,llm,nit)) 
      allocate(ustrorpr_zonal_t(iim,llm,nit))
      allocate(ustrorprfr_zonal_t(iim,llm,nit))

      allocate(upwp_zonal(iim,llm)) 
      allocate(ustror_zonal(iim,llm)) 
      allocate(ustrpr_zonal(iim,llm)) 
      allocate(ustrfr_zonal(iim,llm)) 
      allocate(ustrorpr_zonal(iim,llm)) 
      allocate(ustrorprfr_zonal(iim,llm)) 


   !LOOP ON THE DAYS


      DO it=1,nit



      status = NF90_INQ_VARID(ncid,'upwp',varid)
      status = &
              NF90_GET_VAR(ncid, varid, var3d, &
              start=(/1,1,1,it/),count=(/iim,jjm,klev,1/), &
              stride=(/istri,istri,1,1/))
             upwp_lonlat(:,:,:,it) = var3d(:,:,:)

       status = NF90_INQ_VARID(ncid,'u_orstress',varid)
       status = &
              NF90_GET_VAR(ncid, varid, var3d, &
              start=(/1,1,1,it/),count=(/iim,jjm,klev,1/), &
              stride=(/istri,istri,1,1/))
             zustror_z_lonlat(:,:,:,it) = var3d(:,:,:)

       status = NF90_INQ_VARID(ncid,'u_prstress',varid)
      status = &
              NF90_GET_VAR(ncid, varid, var3d, &
              start=(/1,1,1,it/),count=(/iim,jjm,klev,1/), &
              stride=(/istri,istri,1,1/))
             ustrpr_z_lonlat(:,:,:,it) = var3d(:,:,:)

       status = NF90_INQ_VARID(ncid,'u_frstress',varid)
      status = &
              NF90_GET_VAR(ncid, varid, var3d, &
              start=(/1,1,1,it/),count=(/iim,jjm,klev,1/), &
              stride=(/istri,istri,1,1/))
             ustrfr_z_lonlat(:,:,:,it) = var3d(:,:,:)

       status = NF90_INQ_VARID(ncid,'u_stress_orprcd',varid)
      status = &
              NF90_GET_VAR(ncid, varid, var3d, &
              start=(/1,1,1,it/),count=(/iim,jjm,klev,1/), &
              stride=(/istri,istri,1,1/))
             ustrorpr_z_lonlat(:,:,:,it) = var3d(:,:,:)

       status = NF90_INQ_VARID(ncid,'u_stress_orprfr',varid)
       status = &
              NF90_GET_VAR(ncid, varid, var3d, &
              start=(/1,1,1,it/),count=(/iim,jjm,klev,1/), &
              stride=(/istri,istri,1,1/))
             ustrorprfr_z_lonlat(:,:,:,it) = var3d(:,:,:)

    
!  END OF LOOP ON DAYS

      ENDDO
                 status = nf90_close(ncio)



!  Prepare outputs at grads and netcdf format:
! Creating Netcdf file:
       status = nf90_create("gwd_zonal_lonlat.nc", &
                             nf90_clobber,ncio)
       print *,'Erreur create?',status

include "store_init_zonal.h"

 


   Do it=1,nit
       ! Horizontal averages
        Do ll=1,klev     
           upwp_bar_zt(ll,it) = SUM(upwp_lonlat(:,:,ll,it))/klon_lonlat
           ustror_bar_zt(ll,it) = SUM(ustror_z_lonlat(:,:,ll,it))/klon_lonlat
           ustrpr_bar_zt(ll,it) = SUM(ustrpr_z_lonlat(:,:,ll,it))/klon_lonlat
           ustrfr_bar_zt(ll,it) = SUM(ustrfr_z_lonlat(:,:,ll,it))/klon_lonlat
           ustrorpr_bar_zt(ll,it) = SUM(ustrorpr_z_lonlat(:,:,ll,it))/klon_lonlat
           ustrorprfr_bar_zt(ll,it) = SUM(ustrorprfr_z_lonlat(:,:,ll,it))/klon_lonlat
	Enddo
       ! Zonal averages 
       Do ll=1,klev
         Do ii=1,xlon
             upwp_zonal_t(ii,ll,it) = SUM(upwp_lonlat(ii,:,ll,it))/jjm             
             ustror_zonal_t(ii,ll,it) = SUM(ustror_z_lonlat(ii,:,ll,it))/jjm             
             ustrpr_zonal_t(ii,ll,it) = SUM(ustrpr_z_lonlat(ii,:,ll,it))/jjm             
             ustrfr_zonal_t(ii,ll,it) = SUM(ustrfr_z_lonlat(ii,:,ll,it))/jjm             
             ustrorpr_zonal_t(ii,ll,it) = SUM(ustrorpr_z_lonlat(ii,:,ll,it))/jjm             
             ustrorprfr_zonal_t(ii,ll,it) = SUM(ustrorprfr_z_lonlat(ii,:,ll,it))/jjm 
         Enddo
       Enddo            
   Enddo
   
! time averages

  Do ll=1,klev     
           upwp_bar_z(ll) = SUM(upwp_bar_zt(ll,:))/nit
           ustror_bar_z(ll) = SUM(ustror_bar_zt(ll,:))/nit
           ustrpr_bar_z(ll) = SUM(ustrpr_bar_zt(ll,:))/nit
           ustrfr_bar_z(ll) = SUM(ustrfr_bar_zt(ll,:))/nit
           ustrorpr_bar_z(ll) = SUM(ustrorpr_bar_zt(ll,:))/nit
           ustrorprfr_bar_z(ll) = SUM(ustrorprfr_bar_zt(ll,:))/nit
         Do ii=1,xlon
             upwp_zonal(ii,ll) = SUM(upwp_zonal_t(ii,ll,:))/nit             
             ustror_zonal(ii,ll) = SUM(ustror_zonal_t(ii,ll,:))/nit             
             ustrpr_zonal(ii,ll) = SUM(ustrpr_zonal_t(ii,ll,:))/nit             
             ustrfr_zonal(ii,ll) = SUM(ustrfr_zonal_t(ii,ll,:))/nit             
             ustrorpr_zonal(ii,ll) = SUM(ustrorpr_zonal_t(ii,ll,:))/nit  
             ustrorprfr_zonal(ii,ll) = SUM(ustrorprfr_zonal_t(ii,ll,:))/nit  
         Enddo                      
  Enddo




include  "store_resul_zonal.h"     

      status = nf90_close(ncio)
!==============================IT======================================





  
!       STOP 
        END

