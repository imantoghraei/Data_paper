!     -----------------------------------------------------------------
!*    *COMMON* *YOEGWD* - PARAMETERS FOR GRAVITY WAVE DRAG CALCULATIONS
!     FOR DETAILS SEE THE SUBROUTINE SUGWD
!     ----------------------------------------------------------------
!     IT ALSO CONTAINS 5 CONSTANTS WHICH ARE USED BY THE SCHEME
!     BUT WHICH ARE NOT SPECIFIC OF IT (LIKE RG THE GRAVIT. CTE)
!     -----------------------------------------------------------------
!
      REAL GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT                   &
     &    ,GHMAX,GRAHILO,GSIGCR,GSSEC,GTSEC,GVSEC
      INTEGER NKTOPG,NSTRA
!    
!     PHYSICAL CONSTANT NOT SPECIFIC OF THE SSO SCHEME,
!     BUT WHICH ARE USED BY IT:
      
      REAL RD,RG,RCPD,ROMEGA,RPI,RLVTT

      COMMON/YOEGWD/ GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT         &
     &        ,GHMAX,GRAHILO,GSIGCR,NKTOPG,NSTRA,GSSEC,GTSEC,GVSEC      &
!
     &        ,RD,RG,RCPD,RLVTT,ROMEGA,RPI
!


