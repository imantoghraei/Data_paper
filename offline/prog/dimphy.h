!-----------------------------------------------------------------------
      INTEGER KIDIA, KFDIA, KLON, KLEV
      INTEGER KIDIA_era5, KFDIA_era5, KLON_era5, KLEV_era5
!LOTT PARAMETER (KIDIA_era5=1,KFDIA_era5=iim*(jjm-1)+2,  | IN THE ON-LINE LMD-GCM
!    IN THE OFFLINE VERSION:
      PARAMETER (KIDIA_era5=1,KFDIA_era5=iim*jjm,                                 &
     &           KLON_era5=KFDIA_era5-KIDIA_era5+1,KLEV_era5=llm_era5)

!!!!! IT !!!!!!!!!!!!
  PARAMETER (KIDIA=1,KFDIA=npix,                                 &
     &           KLON=KFDIA-KIDIA+1,KLEV=llm) 
!-----------------------------------------------------------------------
      INTEGER nbtr ! nombre de vrais traceurs
      PARAMETER (nbtr=nqmx-2+1/(nqmx-1))
!-----------------------------------------------------------------------
