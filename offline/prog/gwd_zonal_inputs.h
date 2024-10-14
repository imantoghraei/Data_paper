


   !LOOP ON THE DAYS

	if (iter .GT. 1 .AND. iter .LT. itfin-1) then;
            it1=iter-1*itjum; it2=iter; it3=iter+1*itjum; it4=iter+2*itjum   
        elseif (iter .EQ. 1) then
            it1=iter; it2=iter+1*itjum; it3=iter+2*itjum; it4=iter+3*itjum   
	elseif(iter .GT. itfin-2) then
            it1=iter; it2=iter-1*itjum; it3=iter-2*itjum; it4=iter-3*itjum   
	endif

      do ii=1,iim;do jj=1,jjm;ij=ii+(jj-1)*iim;

             rupwp(ii,jj,:) = rupwp_lonlat(ij,:)
             ustror(ii,jj,:) = zustror_lonlat(ij,:)
             ustrpr(ii,jj,:) = zustrpr_lonlat(ij,:)
             ustrfr(ii,jj,:) = zustrfr_lonlat(ij,:)

             rupwp_f2_1(ii,jj,:) = rupwp_lonlat(ij,:);rupwp_f2_2(ii,jj,:) = rupwp_lonlat(ij,:);
             rupwp_f2_3(ii,jj,:) = rupwp_lonlat(ij,:);rupwp_f2_4(ii,jj,:) = rupwp_lonlat(ij,:);
             rupwp_f2_5(ii,jj,:) = rupwp_lonlat(ij,:);rupwp_f2_6(ii,jj,:) = rupwp_lonlat(ij,:);

             prec12(ii,jj) = (prec(it1,ii,jj)+prec(it2,ii,jj)+prec(it3,ii,jj)+prec(it4,ii,jj))/1000*3600*3
             !prec12(ii,jj) = (prec(iter,ii,jj))
             if (prec12(ii,jj) .GT. prec_cr .AND. filtering .EQ. 1) then
		rupwp(ii,jj,:)=0
                ustror(ii,jj,:)=0; ustrpr(ii,jj,:)=0;ustrfr(ii,jj,:)=0 
            endif

            if (filtering2 .EQ. 1) then
           
             if (prec12(ii,jj) .GE. 0 .AND. prec12(ii,jj) .LT. 0.00000005 ) then
		rupwp_f2_2(ii,jj,:)=0;rupwp_f2_3(ii,jj,:)=0;rupwp_f2_4(ii,jj,:)=0;rupwp_f2_5(ii,jj,:)=0;rupwp_f2_6(ii,jj,:)=0;
            endif
             if (prec12(ii,jj) .GT. 0.00000005 .AND. prec12(ii,jj) .LT. 0.0000001) then
		rupwp_f2_1(ii,jj,:)=0;rupwp_f2_3(ii,jj,:)=0;rupwp_f2_4(ii,jj,:)=0;rupwp_f2_5(ii,jj,:)=0;rupwp_f2_6(ii,jj,:)=0;
            endif
             if (prec12(ii,jj) .GT. 0.0000001 .AND. prec12(ii,jj) .LT. 0.0000005 ) then
		rupwp_f2_1(ii,jj,:)=0;rupwp_f2_2(ii,jj,:)=0;rupwp_f2_4(ii,jj,:)=0;rupwp_f2_5(ii,jj,:)=0;rupwp_f2_6(ii,jj,:)=0;
            endif
             if (prec12(ii,jj) .GT. 0.0000005 .AND. prec12(ii,jj) .LT. 0.000001) then
		rupwp_f2_1(ii,jj,:)=0;rupwp_f2_2(ii,jj,:)=0;rupwp_f2_3(ii,jj,:)=0;rupwp_f2_5(ii,jj,:)=0;rupwp_f2_6(ii,jj,:)=0;
            endif
             if (prec12(ii,jj) .GT. 0.000001 .AND. prec12(ii,jj) .LT. 0.000005) then
		rupwp_f2_1(ii,jj,:)=0;rupwp_f2_2(ii,jj,:)=0;rupwp_f2_3(ii,jj,:)=0;rupwp_f2_4(ii,jj,:)=0;rupwp_f2_6(ii,jj,:)=0;
            endif
             if (prec12(ii,jj) .GT. 0.000005) then
		rupwp_f2_1(ii,jj,:)=0;rupwp_f2_2(ii,jj,:)=0;rupwp_f2_3(ii,jj,:)=0;rupwp_f2_4(ii,jj,:)=0;rupwp_f2_5(ii,jj,:)=0;
            endif

            endif


      enddo;enddo


