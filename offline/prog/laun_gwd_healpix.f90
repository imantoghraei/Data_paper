      PROGRAM LAUN_SSO 

      use netcdf
      use flott_gwd_rando_m, only: flott_gwd_rando
      use acama_gwd_rando_m, only: acama_gwd_rando
!==============================IT======================================
      use healpix_modules   ! IT: Modules necessary to work with the the healpix grids
!======================================================================
      IMPLICIT none
!NETCDF STUFF
      integer ncid,status,varid,dimid
      integer start(4),count(4),stride(4),ilmdz
!NETCDF STUFF FOR OUTPUT FILE:
      integer ncio
      integer londimid,latdimid,levdimid,timdimid
      integer lonvarid,latvarid,levvarid,timvarid
      integer uvarid,vvarid,tvarid,rvarid
      integer healpixdimid, healpixvarid ! IT: Healpix grids
!FORMATTED I/O STUFF
      character*1 entete
      real amoins,aplus,bmoins,bplus
!LOOPS ON TIMES
      integer it,nit,iter,nvar,itdeb,itfin,itjum,istri

!======================================================================
!
! Author(s) F. LOTT            date: 19990624
!
! AN OFFLINE CALL OF THE SSO DRAG SCHEME by LOTT&MILLER(1997),
! THAT SHOULD MAKE EASY TO IMPLEMENT THE SCHEME IN ANY GCM:
! IT IS A SIMPLE INTERFACE TOWARD ROUTINES THAT ARE CALLED
! AS THEY ARE PROVIDED HERE: IN THE LMD-GCM.
! PARTS THAT ARE TYPICAL OFF THE OFFLINE CALL ARE SURROUNDED
! BY COFF-BEG -- COFF-END
!
! IN MOST GCMS YOU CAN FIGURE THAT AT THE STAGE THIS SEQUENCES
! OF ROUTINES ARE CALLED, YOU ARE IN THE ROUTINE THAT CALL
! ALL THE PHYSICAL PARAMETRIZATION ONE AFTER THE OTHER.
!======================================================================
!
!  Variables (INPUT OF THE PHYSICS   IN ON LINE CALLS)
!            (READ ON EXTERNAL FILES IN OFF-LINE CALL)
!
! klon_era5----input-I-Total number of horizontal points that get into physics in era5
! klev_era5----input-I-Number of vertical levels in era5
! paprs---input-R-Pressure in semi-layers   (Pa)
! pplay---input-R-Pressure inside layers    (Pa)
! pphis---input-R-geopotential at the ground
! u-------input-R-Zonal velocity
! v-------input-R-Meridional velocity
! t-------input-R-Temperature (K)
! dtime---input-R-Time step of the physics
! zmea----input-R-Mean Orography (m)
! zstd_era5-input-R-SSO standard deviation (m) in era5
! zsig_era5-input-R-SSO slope in era5
! zgam_era5-input-R-SSO Anisotropy in era5
! zthe_era5-input-R-SSO Angle in era5
! zpic_era5-input-R-SSO Peacks elevation (m) in era5
! zval_era5-input-R-SSO Valleys elevation (m) in era5
! prec----input-R-Total precipitation kg/s/m^2
! zgpcp---input-R-GPCP Precipitations in mm/day, missing value -99999

!==============================IT======================================
! klon-input-I-Total number of horizontal points in the healpix that get into physics 
! klev-input-I-Number of vertical levels in the healpix grids
! zmea-input-R-Mean Orography (m) in the healpix grids
! zstd-input-R-SSO standard deviation (m) in the healpix grids
! zsig-input-R-SSO slope in the healpix grids
! zgam-input-R-SSO Anisotropy in the healpix grids
! zthe-input-R-SSO Angle in the healpix grids
! zpic-input-R-SSO Peacks elevation (m) in the healpix grids
! zval-input-R-SSO Valleys elevation (m) in the healpix grids
! prec-input-R-Total precipitation kg/s/m^2 in the healpix grids
! zgpcp-input-R-GPCP Precipitations in mm/day in the healpix grids
!======================================================================

!------------Modified quantities---------------------------
!
! d_u_oro-output-R-"u" increment due to SSO DRAG     (m/s)
! d_v_oro-output-R-"v" increment due to SSO DRAG     (m/s)
! d_t_oro-output-R-"t" increment due to SSO DRAG     (K)
! d_u_lif-output-R-"u" increment due to MOUNTAIN LIFT(m/s)
! d_v_lif-output-R-"v" increment due to MOUNTAIN LIFT(m/s)
! d_t_lif-output-R-"t" increment due to MOUNTAIN LIFT(K)
! zulow,zvlow -output-R: Low-level wind
! zustrdr,zvstrdr      : Surface stress due to SSO drag      (Pa)
! zustrli,zvstrli      : Surface stress due to MOUNTAIN LIFT (Pa)
!
! igwd--local-I: Total nb of points where the orography schemes are active
! itest-local-I: Flags to indicate active points
! idx---local-I: Locate the physical location of an active point.
!
!------------Issued from commons---------------------------
!
! iim--common-I: Number of longitude intervals
! jjm--common-I: Number of latitude intervals
! klon-common-I: Number of points seen by the physics
!                (iim+1)*(jjm+1) for instance in the healpix grids
! klev-common-I: Number of vertical layers in the healpix grids
!
! npix--common-I: Total number of horizontal healpix grids




include "dimensions.h" ! IT: New dimension.h for the healpix grids
include "dimphy.h"   ! IT: New dimphy.h for the healpix grids
include "YOEGWD.h"
      


       integer ii,jj,ij
       real xlon(iim),ylat(jjm),zpre(llm)
       real xx(iim),yy(jjm)
       real ap(llm),bp(llm),zlev(llm)
       real, allocatable :: var1d(:),var2d(:,:),var3d(:,:,:) 
       !!! IT: for healpix grids :

       real, allocatable :: var1d_vert_healp(:), var1d_hor_healp(:)
       real, allocatable :: var2d_healp(:,:), var3d_healp(:,:,:)

      integer i,j,l,ll,index

!  This is for precips data retrieval

! TYPICALLY USED FOR THE OFFLINE VERSION ONLY
      integer iplay(llm),zplay(llm)
      character*12 input,output
      character*25 filename

!  IMPUT DYNAMICAL FIELDS THAT SHOULD BE PROVIDED OBVIOUSLY
!  PLACED ON THE PHYSICAL GRID:

      integer level(klev)
      REAL paprs(klon,klev+1),pplay(klon,klev),psol(klon)    
      REAL u(klon,klev),v(klon,klev), w(klon,klev)
      REAL ubar(klev),vbar(klev), wbar(klev) ! IT: bar quantities in each layer
      REAL up(klon,klev),vp(klon,klev), wp(klon,klev) ! IT: prime quantities in each layer
      REAL t(klon,klev),rot(klon,klev)
      REAL rupwp(klon,klev), rvptp(klon,klev) ! IT: From ICON data
      REAL EPF(klon,klev), EPF_phi(klon,klev)
      REAL rho_z(klon,klev)  
      REAL upwp(klon,klev),upvp(klon,klev), vptp(klon,klev) ! IT: From ICON data
      REAL zg(klon,klev) , zghalf(klon,klev+1)  ! IT: geometric height from ICON data
      REAL east_acama(klon,klev),west_acama(klon,klev)
      REAL east_flott(klon,klev),west_flott(klon,klev)
      REAL rlat_era5(klon_era5)
      REAL rlat_lonlat(klon_era5)
!  OUTPUT TENDENCIES (NOT USED IN THE OFF-LINE CALL)

      REAL d_u(klon,klev),d_v(klon,klev),d_t(klon,klev)
!
! TIME STEP OF THE PHYSICS, time of the data
!
      REAL dtime
      INTEGER itime
!
! SUBGRID-SCALE OROGRAPHY PARAMETERS
!
      real zstd_era5(klon_era5),zsig_era5(klon_era5),zmea_era5(klon_era5)
      real zgam_era5(klon_era5),zthe_era5(klon_era5)
      real zpic_era5(klon_era5),zval_era5(klon_era5)

!==============================IT======================================
! IT: SUBGRID-SCALE OROGRAPHY PARAMETERS in healpix grids
!
      real zstd(klon),zsig(klon),zmea(klon)
      real zgam(klon),zthe(klon)
      real zpic(klon),zval(klon)
      real rlat(klon)
!======================================================================
!
!  OUTPUT DIAGNOSTICS:
!     zulow(:),zvlow(:):       LOW-LEVEL WIND
!     zustrdr(:), zvstrdr(:):  LOW LEVEL STRESS DUE TO THE DRAG (Pa)
!     zustrli(:), zvstrli(:):  LOW LEVEL STRESS DUE TO THE LIFT (Pa)

      REAL zulow(klon),zvlow(klon),zustrdr(klon)
      REAL zvstrdr(klon), zustrli(klon), zvstrli(klon)
!==============================IT======================================
      REAL ustror_z(klon,klev+1)       ! IT: DRAG STRESS (Pa) as a function of z
      REAL ustrpr_z(klon,klev+1)       ! IT: Precipitation STRESS (Pa) as a function of z
      REAL ustrfr_z(klon,klev+1)       ! IT: Front STRESS (Pa) as a function of z
      REAL dustrordz(klon,klev+1)      ! IT: doro/dz
      REAL dustrprdz(klon,klev+1)      ! IT: dpr/dz
      REAL dustrfrdz(klon,klev+1)      ! IT: dfr/dz
      REAL dEPFdz(klon,klev+1)      ! IT: dEPF/dz
      REAL drupwpdz(klon,klev+1)      ! IT: d(rupwp)/dz
      REAL rad_e
!==============================IT======================================
      INTEGER igwd,igwdim
      INTEGER idx(klon),itest(klon)

! U,V,T TENDENCIES DUE TO SSO DRAG and LIFT

      real d_u_oro(klon,klev),d_v_oro(klon,klev)
      real d_t_oro(klon,klev)
      real d_u_lif(klon,klev), d_v_lif(klon,klev)
      real d_t_lif(klon,klev)
      real d_u_hin(klon,klev), d_v_hin(klon,klev)
      real d_t_hin(klon,klev)
      real d_u_lot(klon,klev), d_v_lot(klon,klev)
      real d_t_lot(klon,klev)
!  Surface stress due to hines!
      real zustrhi(klon),zvstrhi(klon)
!  Surface stress due to Lott!
      real zustrlo(klon),zvstrlo(klon)
      real bvlow(klon)
! Precip 2D fields needed to compute GWs amplitudes
      real prec(klon)

!==============================IT======================================
!!! lonlat fields
      real prec_lonlat(klon_era5)
      real zustrhi_lonlat(klon_era5),zvstrhi_lonlat(klon_era5)
      REAL u_lonlat(klon_era5,klev),v_lonlat(klon_era5,klev), w_lonlat(klon_era5,klev)
      REAL t_lonlat(klon_era5,klev),rot_lonlat(klon_era5,klev)
      REAL rupwp_lonlat(klon_era5,klev), rvptp_lonlat(klon_era5,klev)
      REAL EPF_lonlat(klon_era5,klev),EPF_phi_lonlat(klon_era5,klev)  ! IT: ICON data
      REAL upwp_lonlat(klon_era5,klev), vptp_lonlat(klon_era5,klev)! IT: ICON data
      REAL rho_z_lonlat(klon_era5,klev), pplay_lonlat(klon_era5,klev) 
      real d_u_oro_lonlat(klon_era5,klev),d_v_oro_lonlat(klon_era5,klev)
      real d_u_lif_lonlat(klon_era5,klev), d_v_lif_lonlat(klon_era5,klev)
      real d_u_hin_lonlat(klon_era5,klev), d_v_hin_lonlat(klon_era5,klev)
      real d_u_lot_lonlat(klon_era5,klev), d_v_lot_lonlat(klon_era5,klev)
      REAL east_acama_lonlat(klon_era5,klev),west_acama_lonlat(klon_era5,klev)
      REAL east_flott_lonlat(klon_era5,klev),west_flott_lonlat(klon_era5,klev)
      REAL ustror_z_lonlat(klon_era5,klev)
      REAL ustrpr_z_lonlat(klon_era5,klev)
      REAL ustrfr_z_lonlat(klon_era5,klev)
      REAL dustrordz_lonlat(klon_era5,klev)
      REAL dustrprdz_lonlat(klon_era5,klev)
      REAL dustrfrdz_lonlat(klon_era5,klev)
      REAL dEPFdz_lonlat(klon_era5,klev)      
      REAL dEPFdphi_lonlat(klon_era5,klev)      
      REAL drupwpdz_lonlat(klon_era5,klev)      

!!!! For calculating the vorticity
    ROMEGA = 0.00007292115

!======================================================================

! OFF-DEB
!   THIS IS THE TIME-STEP THE PHYSICS SEE, IT IS
!  PUT TO 30 MINUTES HERE RATHER ABITRARILY

      DTIME=30.*60.
      DTIME=3.*3600.        !!! Healpix data are every 3 hours

!  Load the number of dates you want to proceed:

      read(*,*)itdeb
      read(*,*)itfin
      read(*,*)itjum
      read(*,*)istri

!==============================IT======================================
! Allocating the dimension of the healpix grids to the variables 
      allocate(var1d_hor_healp(npix))

!======================================================================


!  READ THE DYNAMICAL FIELDS OF THE GCM:
!  IN THE LMD-GCM, THE DATA ARE STORED
!  FROM TOP TO BOTTOM
!
!  WHEN GETTING INTO THE PHYSICS PACKAGE ALL THE DATA THAT
!  ARE READ NEXT SHOULD BE KNOWN, IN AN ON-LINE CONTEXT

!  LOAD THE SSO-PARAMETERS
         
        status = NF90_OPEN( &
      'data/sso_era5_hres.nc',NF90_NOWRITE,ncid)
      status = NF90_INQ_VARID(ncid,'longitude',varid)
      status = NF90_INQ_DIMID(ncid,'longitude',dimid)
      status = NF90_INQUIRE_DIMENSION(ncid,dimid,len=ii)
      if(ii.ne.iim) then;print *,'x-Dim SSOs wrong?';
      print *,ii,iim;endif
      status = &
        NF90_GET_VAR(ncid, varid, xlon,start=(/1/),count=(/iim/), &
                     stride=(/istri/))
      print *,'xlon(',iim,'):',xlon(1),xlon(iim/2),xlon(iim)
      xx(:)=xlon(:)*acos(-1.)/180.
      status = NF90_INQ_VARID(ncid,'latitude',varid)
      status = NF90_INQ_DIMID(ncid,'latitude',dimid)
      status = NF90_INQUIRE_DIMENSION(ncid,dimid,len=jj)
      if(jj.ne.jjm) then;print *,'y-Dim SSOs wrong?'
       print *,jj,jjm;endif
      status = &
        NF90_GET_VAR(ncid, varid, ylat,start=(/1/),count=(/jjm/), &
                     stride=(/istri/))
      print *,'ylat(',jjm,'):',ylat(1),ylat(jjm/2),ylat(jjm)
         if(ylat(2).lt.ylat(1))then; print *,'reorganize';stop;endif
      yy(:)=ylat(:)*acos(-1.)/180.
      index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
      rlat_era5(index)=ylat(jj); enddo;enddo          


      allocate(var2d(iim,jjm)) 
      allocate(var3d(iim,jjm,llm)) 
      status = NF90_INQ_VARID(ncid,'zmea',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var2d,start=(/1,1/),count=(/iim,jjm/), &
                     stride=(/istri,istri/))
         index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
         zmea_era5(index)=var2d(ii,jj);enddo;enddo
       status = NF90_INQ_VARID(ncid,'zstd',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var2d,start=(/1,1/),count=(/iim,jjm/), &
                     stride=(/istri,istri/))
         index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
         zstd_era5(index)=var2d(ii,jj);enddo;enddo
       status = NF90_INQ_VARID(ncid,'zsig',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var2d,start=(/1,1/),count=(/iim,jjm/), &
                     stride=(/istri,istri/))
         index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
         zsig_era5(index)=var2d(ii,jj);enddo;enddo
       status = NF90_INQ_VARID(ncid,'zgam',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var2d,start=(/1,1/),count=(/iim,jjm/), &
                     stride=(/istri,istri/))
         index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
         zgam_era5(index)=var2d(ii,jj);enddo;enddo
       status = NF90_INQ_VARID(ncid,'zthe',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var2d,start=(/1,1/),count=(/iim,jjm/), &
                     stride=(/istri,istri/))
         index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
         zthe_era5(index)=var2d(ii,jj);enddo;enddo
       status = NF90_INQ_VARID(ncid,'zpic',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var2d,start=(/1,1/),count=(/iim,jjm/), &
                     stride=(/istri,istri/))
         index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
         zpic_era5(index)=var2d(ii,jj);enddo;enddo
       status = NF90_INQ_VARID(ncid,'zval',varid)
      status = &
        NF90_GET_VAR(ncid, varid, var2d,start=(/1,1/),count=(/iim,jjm/), &
                     stride=(/istri,istri/))
         index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
         zval_era5(index)=var2d(ii,jj);enddo;enddo
         status = NF90_CLOSE(ncid)

         print *,'zmea:',minval(zmea_era5),maxval(zmea_era5) 
         print *,'zstd:',minval(zstd_era5),maxval(zstd_era5) 
         print *,'zsig',minval(zsig_era5),maxval(zsig_era5) 
         print *,'zgam',minval(zgam_era5),maxval(zgam_era5) 
         print *,'zthe',minval(zthe_era5),maxval(zthe_era5) 
         print *,'zpic',minval(zpic_era5),maxval(zpic_era5) 
         print *,'zval',minval(zval_era5),maxval(zval_era5) 

!==============================IT======================================

!! IT: Mapping the SSO-PARAMETERS into the healpix grids using the closest point

        !call era52healp_closest(rlat_era5,rlat)   
        call era52healp_closest(zmea_era5,zmea)
        call era52healp_closest(zstd_era5,zstd)
        call era52healp_closest(zsig_era5,zsig)
        call era52healp_closest(zgam_era5,zgam)
        call era52healp_closest(zthe_era5,zthe)
        call era52healp_closest(zpic_era5,zpic)
        call era52healp_closest(zval_era5,zval)
        call era52healp_closest(rlat_era5,rlat)   
!======================================================================



!   Handle the vertical dimension of ERA5    
!   Read the ap and bp needed to calculate pressure from psol

       status = NF90_OPEN('VITUV_PRESS_TEMP.nc',NF90_NOWRITE,ncid)

       status = NF90_INQ_VARID(ncid,'level_full',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      status = NF90_INQ_DIMID(ncid,'level_full',dimid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
     status = NF90_INQUIRE_DIMENSION(ncid,dimid,len=ll)
      if(ll.ne.klev) then;print *,'pb z-dimension Healpix';stop;endif
       allocate(var2D_healp(npix,llm))
       allocate(var3D_healp(npix,1,llm))
       allocate(var1d_vert_healp(klev))
       status = &
              NF90_GET_VAR(ncid, varid, var1d_vert_healp, &
              start=(/1/),count=(/klev/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      print *,'klev=',klev
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       level(:)=int(var1d_vert_healp(:));print *,'level=',level
!       open(30,file='data/era5_levels.txt')
!       read(30,*)entete !; print *,entete
!       do ll=1,klev;read(30,*)ii,amoins,bmoins;read(30,*)ii,aplus,bplus
!       ap(klev+1-ll)=(amoins+aplus)/2.;bp(klev+1-ll)=(bmoins+bplus)/2.
!       enddo
!       do ll=1,klev;print *,'pf(',ll,')=',(ap(ll)+bp(ll)*101325)/100.
!       enddo
!       do ll=1,klev;zlev(ll)=-7.*log(bp(ll)+ap(ll)/101325.);enddo
!       close(30) 
!       !deallocate(var1d)
!  Temporal dimension of input:
       status = NF90_INQ_VARID(ncid,'time',varid)
       status = NF90_INQ_DIMID(ncid,'time',dimid)
       status = NF90_INQUIRE_DIMENSION(ncid,dimid,len=ll)
       allocate(var1d(ll))
       status = NF90_CLOSE(ncid)
       

!  Prepare outputs at grads and netcdf format:

! This prepare the storage at netdf and grads format:

       open(20,file='gwd_healpix_temp.ctl',form='formatted')
       nvar=0

! Open an healpix file and prepare the output accordingly

! Creating Netcdf file:
       status = nf90_create("gwd_healpix_temp.nc", &
                             nf90_clobber,ncio)
       print *,'Erreur create?',status

include "store_init.h"


!==============================IT======================================

! Creating Netcdf file:
       status = nf90_create("gwd_healpix_temp_lonlat.nc", &
                             nf90_clobber,ncio)
       print *,'Erreur create?',status

include "store_init_lonlat.h"


!==============================IT======================================
! This program handles on month of data only
!==============================IT======================================
!  geometric height at full level center :

      print *,' Geometric height begins:'
      status = NF90_OPEN('VITUV_PRESS_TEMP.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'zg',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
             NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1/),count=(/npix,klev/), &
              stride=(/istri,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
         do ii=1,npix; do ll=1,klev;
      zg(ii,ll)=var2d_healp(ii,klev+1-ll);enddo;enddo;
       status = NF90_CLOSE(ncid)

      do ll=1,klev; zlev(ll) = zg(10398,ll)/1000.; enddo;

      status = NF90_OPEN('VITUV_PRESS_TEMP.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'zghalf',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
             NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1/),count=(/npix,klev/), &
              stride=(/istri,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
         do ii=1,npix; do ll=1,klev+1;
      zghalf(ii,ll)=var2d_healp(ii,klev+1-ll);enddo;enddo;
       status = NF90_CLOSE(ncid)

!==============================IT======================================

!LOOP ON THE DAYS

      iter=0
      do it=itdeb,itfin,itjum
      iter=iter+1

! READ THE Healpix grids PRECIPS DATAS:
     
      status = NF90_OPEN('PREC.nc',NF90_NOWRITE,ncid)
! Precips:
    status = NF90_INQ_VARID(ncid,'pr',varid)
      status = &
              NF90_GET_VAR(ncid, varid, var1d_hor_healp, &
              start=(/1,it/),count=(/npix,1/), &
                     stride=(/istri,1/))
       do ii=1,npix;
! IT: In Healpix, precips are in mm accumulated over s:
       prec(ii)=var1d_hor_healp(ii)    
      enddo 
      if(it.eq.itdeb)print *,'Min and Max de prec:',minval(prec),maxval(prec)

! Time
      status = NF90_INQ_VARID(ncid,'time',varid)
      status = NF90_GET_VAR(ncid, varid, var1d, &
              start=(/it/),count=(/1/))

      print *,'Time in precips:',var1d(1)      
 
      status = NF90_CLOSE(ncid)


!  Read the Healpix vitu file: ! We get ICON data in reverse order vertically because of the way the data stored in ICON
   
      print *,' Vitu begins:'
      status = NF90_OPEN('VITUV_PRESS_TEMP.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'ua',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
              NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1,it/),count=(/npix,klev,1/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      do ii=1,npix; do ll=1,klev;
      u(ii,ll)=var2d_healp(ii,klev+1-ll);enddo;enddo;
      if(it.eq.itdeb)print *,'Min and Max de vitu:',minval(u),maxval(u)
      print *,'Une  de vitu:',u(klon/2,klev/2)
      status = NF90_CLOSE(ncid)
      print *,' Vitu ends:'


!  Read the Healpix vitv file:
   
      print *,' Vitv begins:'
      status = NF90_OPEN('VITUV_PRESS_TEMP.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'va',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
              NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1,it/),count=(/npix,klev,1/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      do ii=1,npix; do ll=1,klev;
      v(ii,ll)=var2d_healp(ii,klev+1-ll);enddo;enddo;
      if(it.eq.itdeb)print *,'Min and Max de vitv:',minval(v),maxval(v)
      print *,'Une  de vitu:',v(klon/2,klev/2)
      status = NF90_CLOSE(ncid)
      print *,' Vitv ends:'


!  Read the Healpix temp file

      print *,' Temp begins:'
      status = NF90_OPEN('VITUV_PRESS_TEMP.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'ta',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
             NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1,it/),count=(/npix,klev,1/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
         do ii=1,npix; do ll=1,klev;
      t(ii,ll)=var2d_healp(ii,klev+1-ll);enddo;enddo;
      if(it.eq.itdeb)print *,'Min and Max de Temp:',minval(t),maxval(t)
      print *,'Une  de Temp:',t(klon/2,klev/2)
       status = NF90_CLOSE(ncid)
      print *,' Temp ends:'
!  Build the pressures

!  We read the pressure data from the healpix grids

      print *,' Press begins:'
      status = NF90_OPEN('VITUV_PRESS_TEMP.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'pfull',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
             NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1,it/),count=(/npix,klev,1/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
         do ii=1,npix; do ll=1,klev;
      pplay(ii,ll)=var2d_healp(ii,klev+1-ll);enddo;enddo;
       status = NF90_CLOSE(ncid)

      print *,' pplay:',pplay(10398,1),pplay(10398,klev/2), &
                         pplay(10398,klev)

!  Pressure at 1/2 Levels:

       paprs(:,klev+1)=0.
       do l=2,klev; paprs(:,l)=0.5*(pplay(:,l)+pplay(:,l-1)); enddo
       paprs(:,1)=pplay(:,1)+0.5*(pplay(:,1)-pplay(:,2))
       paprs(:,klev+1)=pplay(:,klev)+0.5*(pplay(:,klev)-pplay(:,klev-1))

       print *,' paprs:',paprs(10398,2),paprs(10398,klev/2)


      print *,' Press ends:'

 ! Rotationnal       ! The vorticity stored in ICON is the correction direction vertically

         print *,' Vort begins:'
      status = NF90_OPEN('VORT.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'vor_z',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
              NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,it,1/),count=(/npix,1,klev/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      do ii=1,npix; do ll=1,klev;
      rot(ii,ll)=var2d_healp(ii,ll)-2.*ROMEGA*SIN(rlat(ii)*PI/180.);enddo;enddo;
      if(it.eq.itdeb)print *,'Min and Max de rot:',minval(rot),maxval(rot)
      if(it.eq.itdeb)print *,'Une  de rot:',rot(klon/2,klev/2)
      if(it.eq.itdeb)print *,'2  de rot:',rot(klon/2+1,klev/2)
      if(it.eq.itdeb)print *,'3  de rot:',rot(klon/2-1,klev/2)
      status = NF90_CLOSE(ncid)

           print *,' Vort ends:'

!==============================IT======================================      

      print *,' Product winds begins:'

     !!!! Initialize upwp, upvp and vptp
      upwp=0; upvp=0;vptp=0; 
      status = NF90_OPEN('products_winds.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'bar(udiv_prime*w_prime)',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
              NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1,it/),count=(/npix,80,1/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      do ii=1,npix; do ll=1,80;
      upwp(ii,ll)=var2d_healp(ii,80+1-ll);enddo;enddo;
      ! bar{upvp}
      status = NF90_OPEN('products_winds.nc',NF90_NOWRITE,ncid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
       status = NF90_INQ_VARID(ncid,'bar(udiv_prime*vdiv_prime)',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
              NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1,it/),count=(/npix,80,1/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      do ii=1,npix; do ll=1,80;
      upvp(ii,ll)=var2d_healp(ii,80+1-ll);enddo;enddo;
      ! bar{vptp}
      status = NF90_INQ_VARID(ncid,'bar(vdiv_prime*t_prime)',varid)
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
        status = &
              NF90_GET_VAR(ncid, varid, var2d_healp, &
              start=(/1,1,it/),count=(/npix,80,1/), &
              stride=(/istri,1,1/))
      if(status .ne. nf90_noerr) print *, trim(nf90_strerror(status))
      do ii=1,npix; do ll=1,80;
      vptp(ii,ll)=var2d_healp(ii,80+1-ll);enddo;enddo;
      print *,' Product winds ens:'
!==============================IT======================================      

!OFF-END

!  INITIALIZE CONSTANT FOR THE GWD SCHEME
!  ON-LINE SHOULD ONLY BE CALLED AT THE MODEL START UP.

        IF(iter.EQ.1) CALL SUGWD_strato(klon,klev,paprs,pplay)
        
       d_u_lif=0; d_v_lif=0; d_t_lif=0; zustrli=0; zvstrli=0
       d_u_oro=0; d_v_oro=0; d_t_oro=0; zustrdr=0; zvstrdr=0
       d_u_hin=0; d_v_hin=0; d_t_hin=0; zustrhi=0; zvstrhi=0
       d_u_lot=0; d_v_lot=0; d_t_lot=0; zustrlo=0; zvstrlo=0
       east_flott=0.;  west_flott=0.;  east_acama=0.;  west_acama=0.

!      goto 158

!  SSO Parameterization from Lott and Miller (1997)

!  SELECTION  POINTS WHERE THE SCHEME IS ACTIVE
      

        igwd=0
        DO i=1,klon
        itest(i)=0
        IF (((zpic(i)-zmea(i)).GT.100.).AND.(zstd(i).GT.10.0)) THEN
          itest(i)=1
          igwd=igwd+1
          idx(igwd)=i
        ENDIF
        ENDDO
        igwdim=MAX(1,igwd)
 
            print *,'Entree drag'

        CALL drag_noro_strato(1,klon,klev,dtime,paprs,pplay, &
                         zmea,zstd, zsig, zgam,  &
                         zthe,zpic,zval, &
!                        igwd,igwdim,idx,itest, &
                         igwd,idx,itest, &
                         t, u, v, &
                         zulow, zvlow, zustrdr, zvstrdr, &
                         d_t_oro, d_u_oro, d_v_oro)


       

!
!  LIFT parameterization from Lott (1999):
!
        igwd=0
        DO i=1,klon
        itest(i)=0
        IF ((zpic(i)-zmea(i)).GT.100.) THEN
          itest(i)=1
          igwd=igwd+1
          idx(igwd)=i
        ENDIF
        ENDDO
        igwdim=MAX(1,igwd)
!
            print *,'Entree Lift'

        CALL lift_noro_strato(klon,klev,dtime,paprs,pplay, &
                         rlat,zmea,zstd, zsig,  &
                         zgam, zthe,zpic,zval, &
!                        igwd,igwdim,idx,itest, &
                         igwd,idx,itest, &
                         t, u, v, &
                         zulow, zvlow, zustrli, zvstrli, &
                         d_t_lif, d_u_lif, d_v_lif)
!

!
!  Non-orographic GWS due to rain:
!
158    Continue

            print *,'Entree Lott'

 
             call flott_gwd_rando(dtime,pplay, &
                    t, u, v,prec, &
                    zustrlo,zvstrlo, &
                    d_u_lot,d_v_lot,east_flott,west_flott)


            print *,'Entree Acama'

!  Non-orographic GWS due to fronts:


       call acama_gwd_rando(dtime,pplay, rlat, &
                      t, u, v, rot, &
                     zustrhi,zvstrhi, &
                     d_u_hin,d_v_hin,east_acama,west_acama)

          print *,' Sortie Acama'

!OFF-BEG        

!==============================IT======================================   
        print *,'Entree Stress calculation'

        ustror_z(:,klev+1)=0.
        ustrpr_z(:,klev+1)=0.
        ustrfr_z(:,klev+1)=0.
        DO ll=klev,1,-1
             rho_z(:,ll) = -1/rg*(paprs(:,ll+1)-paprs(:,ll))/ &
                           (zghalf(:,ll+1)-zghalf(:,ll))
            
             ustror_z(:,ll) = ustror_z(:,ll+1) - d_u_oro(:,ll)/ &
                                dtime*1/rg*(paprs(:,ll+1)-paprs(:,ll))
             ustrpr_z(:,ll) = ustrpr_z(:,ll+1) - d_u_lot(:,ll)/ &
                                dtime*1/rg*(paprs(:,ll+1)-paprs(:,ll))
             ustrfr_z(:,ll) = ustrfr_z(:,ll+1) - d_u_hin(:,ll)/ &
                                dtime*1/rg*(paprs(:,ll+1)-paprs(:,ll))
	     rupwp(:,ll) = upwp(:,ll)  * rho_z(:,ll)
	     rvptp(:,ll) = vptp(:,ll) / t(:,ll) * (2.*ROMEGA*SIN(rlat(:)*PI/180.)) * rho_z(:,ll)
             EPF(:,ll) = rupwp(:,ll) - rvptp(:,ll)
             EPF_phi(:,ll) = - upvp(:,ll)* rho_z(:,ll)
        ENDDO

        DO ll=klev,1,-1
             dustrordz(:,ll) =-1/rho_z(:,ll)* (ustror_z(:,ll+1)-ustror_z(:,ll))/ (zghalf(:,ll+1)-zghalf(:,ll)) *24.*3600.
             dustrprdz(:,ll) =-1/rho_z(:,ll)* (ustrpr_z(:,ll+1)-ustrpr_z(:,ll))/ (zghalf(:,ll+1)-zghalf(:,ll)) *24.*3600.
             dustrfrdz(:,ll) =-1/rho_z(:,ll)* (ustrfr_z(:,ll+1)-ustrfr_z(:,ll))/ (zghalf(:,ll+1)-zghalf(:,ll)) *24.*3600.
	     dEPFdz(:,ll) = -1/rho_z(:,ll)*(EPF(:,ll+1)-EPF(:,ll))/ (zghalf(:,ll+1)-zghalf(:,ll)) *24.*3600.
	     drupwpdz(:,ll) = -1/rho_z(:,ll)*(rupwp(:,ll+1)-rupwp(:,ll))/ (zghalf(:,ll+1)-zghalf(:,ll)) *24.*3600.
        ENDDO
          print *,' Sortie Stress calculation'
      if(it.eq.itdeb)print *,'Min and Max de uptp:',minval(rupwp),maxval(rupwp)
      if(it.eq.itdeb)print *,'Min and Max de vptp:',minval(rvptp),maxval(rvptp)
!==============================IT======================================


         print *,'Sortie Acama....'

!  TENDENCIES IN M/S/DAY

       do i=1,klon
       do l=1,klev
          d_u_oro(i,l)=d_u_oro(i,l)*24.*3600./dtime
          d_v_oro(i,l)=d_v_oro(i,l)*24.*3600./dtime
          d_u_lif(i,l)=d_u_lif(i,l)*24.*3600./dtime
          d_v_lif(i,l)=d_v_lif(i,l)*24.*3600./dtime
          d_u_hin(i,l)=d_u_hin(i,l)*24.*3600./dtime
          d_v_hin(i,l)=d_v_hin(i,l)*24.*3600./dtime
          d_u_lot(i,l)=d_u_lot(i,l)*24.*3600./dtime
          d_v_lot(i,l)=d_v_lot(i,l)*24.*3600./dtime

       enddo
       enddo





!  NETCDF STORAGE:      

!     Here storage
      
      status = nf90_open('gwd_healpix_temp.nc',NF90_WRITE,ncio)

      print *,'Erreur open?',status


      print *,'storage begins....'
include  "store_resul.h"

      status = nf90_close(ncio)
!==============================IT======================================



!! IT: Mapping from the healpix grids into the longitude-latitude grids
!!     using the closest point

        call healp2lonlat_closest_2D(prec,prec_lonlat)   
        call healp2lonlat_closest_2D(zustrhi,zustrhi_lonlat)
        call healp2lonlat_closest_2D(zvstrhi,zvstrhi_lonlat)
        call healp2lonlat_closest(u,u_lonlat)
        call healp2lonlat_closest(v,v_lonlat)
        call healp2lonlat_closest(w,w_lonlat)
        call healp2lonlat_closest(t,t_lonlat)
        call healp2lonlat_closest(rot,rot_lonlat)
        call healp2lonlat_closest(pplay,pplay_lonlat)
        call healp2lonlat_closest(rupwp,rupwp_lonlat)
        call healp2lonlat_closest(rvptp,rvptp_lonlat)
        call healp2lonlat_closest(upwp,upwp_lonlat)
        call healp2lonlat_closest(vptp,vptp_lonlat)
        call healp2lonlat_closest(EPF,EPF_lonlat)
        call healp2lonlat_closest(EPF_phi,EPF_phi_lonlat)
        call healp2lonlat_closest(rho_z,rho_z_lonlat)
        call healp2lonlat_closest(d_u_oro,d_u_oro_lonlat)
        call healp2lonlat_closest(d_v_oro,d_v_oro_lonlat)
        call healp2lonlat_closest(d_u_lot,d_u_lot_lonlat)
        call healp2lonlat_closest(d_u_hin,d_u_hin_lonlat)
        call healp2lonlat_closest(east_acama,east_acama_lonlat)
        call healp2lonlat_closest(east_flott,east_flott_lonlat)
        call healp2lonlat_closest(west_acama,west_acama_lonlat)
        call healp2lonlat_closest(west_flott,west_flott_lonlat)
        call healp2lonlat_closest(ustror_z,ustror_z_lonlat)
        call healp2lonlat_closest(ustrpr_z,ustrpr_z_lonlat)
        call healp2lonlat_closest(ustrfr_z,ustrfr_z_lonlat)
        call healp2lonlat_closest(dustrordz,dustrordz_lonlat)
        call healp2lonlat_closest(dustrprdz,dustrprdz_lonlat)
        call healp2lonlat_closest(dustrfrdz,dustrfrdz_lonlat)
        call healp2lonlat_closest(dEPFdz,dEPFdz_lonlat)
        call healp2lonlat_closest(drupwpdz,drupwpdz_lonlat)
        call healp2lonlat_closest_2D(rlat,rlat_lonlat)
!==============================IT======================================
      print *,'Meridional derivative of EPF begins....'
         call mer_deriv(ylat,EPF_phi_lonlat,dEPFdphi_lonlat)
      print *,'Meridional derivative of EPF ended....'
      rad_e =  6371000.0
      do ii = 1, iim;do jj = 1, jjm;ij = ii + (jj-1) * iim
          dEPFdphi_lonlat(ij,:) = 1/(rho_z_lonlat(ij,:)*rad_e*COS(ylat(jj)*PI/180.))*dEPFdphi_lonlat(ij,:)*24.*3600.
      enddo;enddo;
!==============================IT======================================

 ! IT: We also store the data in lonlat grids
      status = nf90_open('gwd_healpix_temp_lonlat.nc',NF90_WRITE,ncio)

      print *,'Erreur open?',status


      print *,'storage in lonlat begins....'

include  "store_resul_lonlat.h"     

      status = nf90_close(ncio)


!  END OF LOOP ON DAYS

        ENDDO




!  Preparing ctl file for grads:

        open(20,file='gwd_healpix_temp.ctl', &
                form='formatted')
        write(20,220)'FORMAT ','sequential'
        write(20,220)'UNDEF  ','1.0E30'
220     format(a7,a28)

        write(20,221) 'XDEF',iim,' LINEAR ',xlon(1), &
                                    xlon(2)-xlon(1)
        write(20,2210) 'YDEF',jjm,'LINEAR ',ylat(1),ylat(2)-ylat(1)
        !write(20,224) 'ZDEF',klev,' LEVELS ', (zlev(l),l=1,klev)
        write(20,222) 'TDEF',iter,' LINEAR','00Z01JAN1960','24hr'
        write(20,223) 'VARS      ',nvar
        write(20,223) 'prec ',   0,'t,y,x','prec'
        write(20,223) 'zustrhi ',   0,'t,y,x','zustrhi'
        write(20,223) 'zvstrhi ',   0,'t,y,x','zvstrhi'
        write(20,223) 'vitu ',klev,'t,z,y,x','vitu'
        write(20,223) 'vitv ',klev,'t,z,y,x','vitv'
        write(20,223) 'temp ',klev,'t,z,y,x','temp'
        write(20,223) 'rot ',klev,'t,z,y,x','rot'
        write(20,223) 'd_u_mount ',klev,'t,z,y,x','du_mount'
        write(20,223) 'd_u_preci ',klev,'t,z,y,x','du_preci'
        write(20,223) 'd_u_front ',klev,'t,z,y,x','du_front'
        write(20,223) 'east_gwstr',klev,'t,z,y,x','east_gwstr'
        write(20,223) 'west_gwstr',klev,'t,z,y,x','west_gwstr'
        write(20,223) 'ENDVARS   ' 

221     format(a4,1x,i3,a8,f6.0,1x,f7.2)
2210    format(a4,1x,i3,a8,96(1x,f7.3)) 
222     format(a4,1x,i3,a8,1x,a12,1x,a6)
223     format(a12,3x,i3,3x,a7,1x,a12)
!224     format(a4,1x,i3,a8,19(1x,i5))
224     format(a4,1x,i3,a8,70(1x,f6.2))

        Print *,'Finishes right?'
!OFF-END
        
!       STOP 
        END


!==============================IT======================================


subroutine era52healp_mean(varin,varout)
	use healpix_modules
	include "dimensions.h" 

	double precision :: theta, phi
	real :: i, j
	real, intent(in) :: varin(klon_era5)
	real, intent(out) :: varout(npix)

	integer :: ii, i1, i2, j1, j2, val1, val2, val3, val4
	real :: delta_i1, delta_i2, delta_j1, delta_j2

        DO ii = 1, npix

		call pix2ang_nest(nside, ii-1,theta, phi)
		j = theta / PI * 180
		i = phi / PI * 180

		if ( anint(i) .gt. i ) then; i2 = anint(i); i1 = i2 - 1;
        	else; i1 = anint(i); i2 = i1 + 1; endif;
		if ( anint(j) .gt. j ) then; j2 = anint(j); j1 = j2 - 1;
        	else; j1 = anint(j); j2 = j1 + 1; endif;
	
		delta_i2 = i2 - i
		delta_i1 = i - i1
		delta_j2 = j2 - j
		delta_j1 = j - j1

		val1 = iim * (j1) + i1; 
		val2 = iim * (j1) + i2; 
		val3 = iim * (j2) + i1; 
		val4 = iim * (j2) + i2; 


		varout(ii) =  delta_j1*(delta_i2 * varin (val3) &
		+ delta_i1* varin (val4) ) + delta_j2*(delta_i2* varin (val1) & 
		+ delta_i1* varin (val2) )

	ENDDO

end subroutine era52healp_mean

subroutine era52healp_closest(varin,varout)
	use healpix_modules
	include "dimensions.h" 
        include "dimphy.h"  
	double precision :: theta, phi
	real :: i, j
	real, intent(in) :: varin(klon_era5)
	real, intent(out) :: varout(klon)

	integer :: ii, i1, i2, j1, j2, val1, val2, val3, val4
	real :: delta_i1, delta_i2, delta_j1, delta_j2
	real :: r1, r2, r3, r4, r(4)
	integer :: rm,  val(4)

        DO ii = 1, npix

		call pix2ang_nest(nside, ii-1,theta, phi)
		j = theta / PI * 180
		i = phi / PI * 180

		if ( anint(i) .gt. i ) then; i2 = anint(i); i1 = i2 - 1;
        	else; i1 = anint(i); i2 = i1 + 1; endif;
		if ( anint(j) .gt. j ) then; j2 = anint(j); j1 = j2 - 1;
        	else; j1 = anint(j); j2 = j1 + 1; endif;
	
		delta_i2 = i2 - i
		delta_i1 = i - i1
		delta_j2 = j2 - j
		delta_j1 = j - j1


		val1 = iim * (j1) + i1; 
		val2 = iim * (j1) + i2; 
		val3 = iim * (j2) + i1; 
		val4 = iim * (j2) + i2; 
		val = [val1, val2, val3, val4]

	
		r1 = delta_i1**2 + delta_j1**2
		r2 = delta_i2**2 + delta_j1**2
		r3 = delta_i1**2 + delta_j2**2 
		r4 = delta_i2**2 + delta_j2**2 
		r = [r1, r2, r3, r4]
		rm = MINLOC(r, DIM=1)



		varout(ii) =  varin (val(rm))

	ENDDO

end subroutine era52healp_closest

subroutine era52healp_closest_reversed(varin,varout)
	use healpix_modules
	include "dimensions.h" 
        include "dimphy.h"  
	double precision :: theta, phi
	real :: i, j
	real, intent(in) :: varin(klon_era5)
	real, intent(out) :: varout(klon)

	integer :: ii, i1, i2, j1, j2, val1, val2, val3, val4
	real :: delta_i1, delta_i2, delta_j1, delta_j2
	real :: r1, r2, r3, r4, r(4)
	integer :: rm,  val(4)

        DO ii = 1, npix

		call pix2ang_nest(nside, ii-1,theta, phi)
		j = 180. - theta / PI * 180
		i = phi / PI * 180

		if ( anint(i) .gt. i ) then; i2 = anint(i); i1 = i2 - 1;
        	else; i1 = anint(i); i2 = i1 + 1; endif;
		if ( anint(j) .gt. j ) then; j2 = anint(j); j1 = j2 - 1;
        	else; j1 = anint(j); j2 = j1 + 1; endif;
	
		delta_i2 = i2 - i
		delta_i1 = i - i1
		delta_j2 = j2 - j
		delta_j1 = j - j1


		val1 = iim * (j1) + i1; 
		val2 = iim * (j1) + i2; 
		val3 = iim * (j2) + i1; 
		val4 = iim * (j2) + i2; 
		val = [val1, val2, val3, val4]

	
		r1 = delta_i1**2 + delta_j1**2
		r2 = delta_i2**2 + delta_j1**2
		r3 = delta_i1**2 + delta_j2**2 
		r4 = delta_i2**2 + delta_j2**2 
		r = [r1, r2, r3, r4]
		rm = MINLOC(r, DIM=1)



		varout(ii) =  varin (val(rm))

	ENDDO

end subroutine era52healp_closest_reversed

subroutine healp2lonlat_closest_2D(varin,varout)
	use healpix_modules
	include "dimensions.h" 
        include "dimphy.h"   
 
	double precision :: theta(klon_era5), phi(klon_era5)
	real, intent(in) :: varin(klon)
	real, intent(out) :: varout(klon_era5)

	integer :: ii, jj, ll
	integer :: index, ipnest

        index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
 	theta(index) = (jj-1)*PI/180  ; phi(index) = (ii-1)*PI/180 ; enddo;enddo;


	DO index = 1,klon_era5
		call ang2pix_nest(nside,theta(index), phi(index), ipnest)

		varout(index) =  varin (ipnest+1)
        ENDDO

end subroutine healp2lonlat_closest_2D

subroutine healp2lonlat_closest(varin,varout)
	use healpix_modules
	include "dimensions.h" 
        include "dimphy.h"   
 
	double precision :: theta(klon_era5), phi(klon_era5)
	real, intent(in) :: varin(klon,klev)
	real, intent(out) :: varout(klon_era5,klev)

	integer :: ii, jj, ll
	integer :: index, ipnest

        index=0;do jj=jjm,1,-1;do ii=1,iim;index=index+1
 	theta(index) = (jj-1)*PI/180  ; phi(index) = (ii-1)*PI/180 ; enddo;enddo;

     DO ll=1, klev
	DO index = 1,klon_era5
		call ang2pix_nest(nside,theta(index), phi(index), ipnest)

		varout(index,ll) =  varin (ipnest+1,ll)
	ENDDO
     ENDDO

end subroutine healp2lonlat_closest


subroutine mer_deriv(ylat, varin, dvarindphi)

include "dimensions.h" ! IT: New dimension.h for the healpix grids
include "dimphy.h"   ! IT: New dimphy.h for the healpix grids
include "YOEGWD.h"
      
      integer ii, jj, ij
!  INPUT DYNAMICAL FIELDS
!  PLACED ON THE PHYSICAL GRID:


      real, intent(in) :: ylat(jjm) !!! latitude grids [degree]
      REAL, intent(in)  :: varin(klon_era5,klev)   
      REAL, intent(out) :: dvarindphi(klon_era5,klev)  
      REAL              :: varin_ii_jj(iim,jjm,klev), rlat_ii_jj(iim,jjm,klev)  
      REAL              :: dy
      REAL :: ylat_rad(jjm)
      REAL, PARAMETER :: pi = 3.141592653589793

      ! Converting varin into longitude (ii) and latitude (jj) grids
      do ii=1,iim
          do jj=1,jjm
              ij = ii + (jj-1) * iim
              varin_ii_jj(ii,jj,:) = varin(ij,:)
          end do
      end do

      ! Converting ylat to radians
      ylat_rad=0.0
      do jj = 1, jjm
          ylat_rad(jj) = ylat(jj) * PI / 180.0
      end do


      ! Calculating meridional derivative using finite differences
      do ii = 1, iim;do jj = 2, jjm-1
         dy = ylat_rad(jj+1) - ylat_rad(jj-1)
         dvarindphi(ii+(jj-1)*iim, :) = (COS(ylat_rad(jj+1))*COS(ylat_rad(jj+1))* varin_ii_jj(ii,jj+1,:) - &
                                COS(ylat_rad(jj-1))*COS(ylat_rad(jj-1))* varin_ii_jj(ii,jj-1,:))/dy/COS(ylat_rad(jj))

      end do;end do

      ! Boundary conditions (simple one-sided differences)
      do ii = 1, iim
          ! Southern boundary
         dy = ylat_rad(2) - ylat_rad(1)
         dvarindphi(ii+(1-1)*iim, :) = (COS(ylat_rad(2))*COS(ylat_rad(2))* varin_ii_jj(ii,2,:) - &
                                COS(ylat_rad(1))*COS(ylat_rad(1))* varin_ii_jj(ii,1,:))/dy/COS(ylat_rad(1))

          ! Northern boundary
         dy = ylat_rad(jjm) - ylat_rad(jjm-1)
         dvarindphi(ii+(jjm-1)*iim, :) = (COS(ylat_rad(jjm))*COS(ylat_rad(jjm))* varin_ii_jj(ii,jjm,:) - &
                               COS(ylat_rad(jjm-1))*COS(ylat_rad(jjm-1))* varin_ii_jj(ii,jjm-1,:))/dy/COS(ylat_rad(jjm))
      end do

end subroutine mer_deriv


!======================================================================
