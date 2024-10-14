      rupwp=0.
      ustror=0.
      ustrpr=0.
      ustrfr=0.

      rupwp_or=0.
      ustror_or=0.
      ustrpr_or=0.
      ustrfr_or=0.

      rupwp_nor_land=0.
      ustror_nor_land=0.
      ustrpr_nor_land=0.
      ustrfr_nor_land=0.

      rupwp_nor_sea=0.
      ustror_nor_sea=0.
      ustrpr_nor_sea=0.
      ustrfr_nor_sea=0.

      rupwp_bar_z=0.  
      ustror_bar_z=0.
      ustrpr_bar_z=0.
      ustrfr_bar_z=0.

      rupwp_zonal=0.  
      ustror_zonal=0.
      ustrpr_zonal=0.
      ustrfr_zonal=0.

      rupwp_or_bar_z  =0.
      ustror_or_bar_z=0.
      ustrpr_or_bar_z=0.
      ustrfr_or_bar_z=0.

      rupwp_or_zonal=0.  
      ustror_or_zonal=0.
      ustrpr_or_zonal=0.
      ustrfr_or_zonal=0.

      rupwp_nor_land_bar_z  =0.
      ustror_nor_land_bar_z=0.
      ustrpr_nor_land_bar_z=0.
      ustrfr_nor_land_bar_z=0.

      rupwp_nor_land_zonal  =0.
      ustror_nor_land_zonal=0.
      ustrpr_nor_land_zonal=0.
      ustrfr_nor_land_zonal=0.

      rupwp_nor_sea_bar_z  =0.
      ustror_nor_sea_bar_z=0.
      ustrpr_nor_sea_bar_z=0.
      ustrfr_nor_sea_bar_z=0.

      rupwp_nor_sea_zonal  =0.
      ustror_nor_sea_zonal=0.
      ustrpr_nor_sea_zonal=0.
      ustrfr_nor_sea_zonal=0.

      rupwp_mean  =0.
      ustror_mean=0.
      ustrpr_mean=0.
      ustrfr_mean=0.

      rupwp_f2_1=0.;rupwp_f2_2=0.;rupwp_f2_3=0.;rupwp_f2_4=0.;rupwp_f2_5=0.;rupwp_f2_6=0.;
!======================================================================

             do ll=1,llm;do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;
             if (rupwp_lonlat(ij,ll) .GT. 0) then
                rupwp(ii,jj,ll) = rupwp_lonlat(ij,ll)
                rupwp_f2_1(ii,jj,ll) = rupwp_lonlat(ij,ll);rupwp_f2_2(ii,jj,ll) = rupwp_lonlat(ij,ll);
                rupwp_f2_3(ii,jj,ll) = rupwp_lonlat(ij,ll);rupwp_f2_4(ii,jj,ll) = rupwp_lonlat(ij,ll);
                rupwp_f2_5(ii,jj,ll) = rupwp_lonlat(ij,ll);rupwp_f2_6(ii,jj,ll) = rupwp_lonlat(ij,ll);
             endif
             if (zustror_lonlat(ij,ll) .GT. 0) then
                ustror(ii,jj,ll) = zustror_lonlat(ij,ll)
             endif
             if (zustrpr_lonlat(ij,ll) .GT. 0) then
                ustrpr(ii,jj,ll) = zustrpr_lonlat(ij,ll)
             endif
             if (zustrfr_lonlat(ij,ll) .GT. 0) then
                ustrfr(ii,jj,ll) = zustrfr_lonlat(ij,ll)
             endif
             if (prec12(ii,jj) .GT. prec_cr .AND. filtering.EQ.1) then
		rupwp(ii,jj,ll)=0;
                ustror(ii,jj,ll)=0;ustrpr(ii,jj,ll)=0;ustrfr(ii,jj,ll)=0 
             endif

            if (filtering2 .EQ. 1) then

             if (prec12(ii,jj) .GE. 0 .AND. prec12(ii,jj) .LT. 0.00000005 ) then
		rupwp_f2_2(ii,jj,ll)=0;rupwp_f2_3(ii,jj,ll)=0;rupwp_f2_4(ii,jj,ll)=0;rupwp_f2_5(ii,jj,ll)=0;rupwp_f2_6(ii,jj,ll)=0;
            endif
             if (prec12(ii,jj) .GT. 0.00000005 .AND. prec12(ii,jj) .LT. 0.0000001) then
		rupwp_f2_1(ii,jj,ll)=0;rupwp_f2_3(ii,jj,ll)=0;rupwp_f2_4(ii,jj,ll)=0;rupwp_f2_5(ii,jj,ll)=0;rupwp_f2_6(ii,jj,ll)=0;
            endif
             if (prec12(ii,jj) .GT. 0.0000001 .AND. prec12(ii,jj) .LT. 0.0000005 ) then
		rupwp_f2_1(ii,jj,ll)=0;rupwp_f2_2(ii,jj,ll)=0;rupwp_f2_4(ii,jj,ll)=0;rupwp_f2_5(ii,jj,ll)=0;rupwp_f2_6(ii,jj,ll)=0;
            endif
             if (prec12(ii,jj) .GT. 0.0000005 .AND. prec12(ii,jj) .LT. 0.000001) then
		rupwp_f2_1(ii,jj,ll)=0;rupwp_f2_2(ii,jj,ll)=0;rupwp_f2_3(ii,jj,ll)=0;rupwp_f2_5(ii,jj,ll)=0;rupwp_f2_6(ii,jj,ll)=0;
            endif
             if (prec12(ii,jj) .GT. 0.000001 .AND. prec12(ii,jj) .LT. 0.000005) then
		rupwp_f2_1(ii,jj,ll)=0;rupwp_f2_2(ii,jj,ll)=0;rupwp_f2_3(ii,jj,ll)=0;rupwp_f2_4(ii,jj,ll)=0;rupwp_f2_6(ii,jj,ll)=0;
            endif
             if (prec12(ii,jj) .GT. 0.000005) then
		rupwp_f2_1(ii,jj,ll)=0;rupwp_f2_2(ii,jj,ll)=0;rupwp_f2_3(ii,jj,ll)=0;rupwp_f2_4(ii,jj,ll)=0;rupwp_f2_5(ii,jj,ll)=0;
            endif

            endif
             enddo;enddo;enddo
