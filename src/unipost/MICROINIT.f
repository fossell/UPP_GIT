      SUBROUTINE MICROINIT(imp_physics)
!
!-- ABSTRACT:
!     Initializes arrays for new cloud microphysics
!
!-- Program History Log:
!     02-02-08  B. Ferrier
!     04-11-19 H CHUANG - WRF VERSION
!
!-- Input argument list:
!     None
!
!-- Output argument list:
!     None
!
!-- Subprograms called:
!     Function FPVS
!
!-- Common blocks:
!     CMASSI
!     RMASS_TABLES
!     MAPOT
!     CRHgrd
!
!-- Attributes:
!     Language: FORTRAN 90
!     Machine : IBM SP
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use rhgrd_mod
     use params_mod, only: tfrz, pi
     use cmassi_mod, only: dmrmax, t_ice, nlimax, flarge2, xmrmax, &
                           mdrmax, mdrmin, trad_ice, massi, &
                           rqr_drmin, n0r0, rqr_drmax, cn0r0, &
                           cn0r_dmrmin, cn0r_dmrmax, dmrmin
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      REAL, PARAMETER :: RHOL=1000.
      real ax,C_N0r0
      integer i
      integer, intent(in):: imp_physics
      real, allocatable:: MASSR(:)
      character filename*80
!
!------------------------ START EXECUTION ------------------------
!
!---  READ IN MASSI FROM LOOKUP TABLES 
!
      if(imp_physics==5 .or. imp_physics==85)then
        if(imp_physics==5)then
           RHgrd=0.98
           filename="eta_micro_lookup.dat"
        else
           RHgrd=1.
           filename="nam_micro_lookup.dat"
        endif
        DMRmax=1.E-3
      else  !-- Should be imp_physics==95
        RHgrd=1.
        DMRmax=.45E-3
        NLImax=5.E3
        FLARGE2=0.03
        filename="nam_micro_lookup.dat"
      end if 
      T_ICE=-40.
      XMRmax=1.E6*DMRmax 
      MDRmax=XMRmax
      allocate(MASSR(MDRmin:MDRmax))
      TRAD_ice=0.5*T_ICE+TFRZ
      
      OPEN (UNIT=1,FILE=filename,convert='big_endian',FORM="UNFORMATTED")
      DO I=1,3
        READ(1)
      ENDDO
      READ(1) MASSR
      DO I=1,5
        READ(1)
      ENDDO
      READ(1) MASSI
      CLOSE(1)
      RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
      RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
!      PI=ACOS(-1.) ! defined in params now
      C_N0r0=PI*RHOL*N0r0
      CN0r0=1.E6/SQRT(SQRT(C_N0r0))
      CN0r_DMRmin=1./(PI*RHOL*DMRmin*DMRmin*DMRmin*DMRmin)
      CN0r_DMRmax=1./(PI*RHOL*DMRmax*DMRmax*DMRmax*DMRmax)
      print *,'MICROINIT: MDRmin, MASSR(MDRmin)=',MDRmin,MASSR(MDRmin)
      print *,'MICROINIT: MDRmax, MASSR(MDRmax)=',MDRmax,MASSR(MDRmax)
!      print *,  'ETA2P:MASSI(50)= ', MASSI(50)
!      print *,  'ETA2P:MASSI(450)= ', MASSI(450)
!      print *,  'ETA2P:MASSI(1000)= ', MASSI(1000)
!
!--- Initialize saturation vapor pressure lookup tables (functions FPVS, FPVS0)
!
      CALL GPVS
!
      deallocate(MASSR)
!--- 
      RETURN
      END
