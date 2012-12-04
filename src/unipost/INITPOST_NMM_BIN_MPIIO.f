      SUBROUTINE INITPOST_NMM_BIN_MPIIO
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    INITPOST    INITIALIZE POST FOR RUN
!   PRGRMMR: RUSS TREADON    ORG: W/NP2      DATE: 93-11-10
!     
! ABSTRACT:  THIS ROUTINE INITIALIZES CONSTANTS AND
!   VARIABLES AT THE START OF AN ETA MODEL OR POST 
!   PROCESSOR RUN.
!
!   THIS ROUTINE ASSUMES THAT INTEGERS AND REALS ARE THE SAME SIZE
!   .     
!     
! PROGRAM HISTORY LOG:
!   93-11-10  RUSS TREADON - ADDED DOCBLOC
!   98-05-29  BLACK - CONVERSION OF POST CODE FROM 1-D TO 2-D
!   99-01 20  TUCCILLO - MPI VERSION
!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!   02-06-19  MIKE BALDWIN - WRF VERSION
!   02-08-15  H CHUANG - UNIT CORRECTION AND GENERALIZE PROJECTION OPTIONS
!   02-10-31  H CHUANG - MODIFY TO READ WRF BINARY OUTPUT
!   05-12-05  H CHUANG - ADD CAPABILITY TO OUTPUT OFF-HOUR FORECAST WHICH HAS
!               NO INPACTS ON ON-HOUR FORECAST
!     
! USAGE:    CALL INIT
!   INPUT ARGUMENT LIST:
!     NONE     
!
!   OUTPUT ARGUMENT LIST: 
!     NONE
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!       NONE
!     LIBRARY:
!       COMMON   - CTLBLK
!                  LOOKUP
!                  SOILDEPTH
!
!    
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
      use vrbls3d
      use vrbls2d
      use soil
      use masks
      use kinds, only             : i_llong
      use wrf_io_flags_mod
      use params_mod
      use lookup_mod
      use ctlblk_mod
      use gridspec_mod
      use module_io_int_idx, only: io_int_index, r_info
      use module_io_int_read, only: io_int_fetch_data
!
!     INCLUDE/SET PARAMETERS.
!     
      INCLUDE "mpif.h"

      character(len=31) :: VarName
      integer :: Status
      character startdate*19,SysDepInfo*80,cgar*1
      character startdate2(19)*4
! 
!     NOTE: SOME INTEGER VARIABLES ARE READ INTO DUMMY ( A REAL ). THIS IS OK
!     AS LONG AS REALS AND INTEGERS ARE THE SAME SIZE.
!
!     ALSO, EXTRACT IS CALLED WITH DUMMY ( A REAL ) EVEN WHEN THE NUMBERS ARE
!     INTEGERS - THIS IS OK AS LONG AS INTEGERS AND REALS ARE THE SAME SIZE.
      LOGICAL RUNB,SINGLRST,SUBPOST,NEST,HYDRO
      LOGICAL IOOMG,IOALL
      CHARACTER*32 LABEL
      CHARACTER*40 CONTRL,FILALL,FILMST,FILTMP,FILTKE,FILUNV                  &
         , FILCLD,FILRAD,FILSFC
      CHARACTER*4 RESTHR
      CHARACTER FNAME*80,ENVAR*50,BLANK*4
      INTEGER IDATB(3),IDATE(8),JDATE(8)
!     
!     DECLARE VARIABLES.
!     
      REAL SLDPTH2(NSOIL)
      REAL RINC(5)
      REAL ETA1(LM), ETA2(LM)
      REAL DUM1D (LM+1)
      REAL DUMMY ( IM, JM )
      REAL DUMMY2 ( IM, JM )
      REAL FI(IM,JM,2)
      INTEGER IDUMMY ( IM, JM )
      REAL DUM3D ( IM, LM, JM )
      REAL DUM3D2 ( IM, LM+1, JM ),DUMSOIL ( IM, NSOIL, JM )

      integer ibuf(im,jsta_2l:jend_2u)
      real buf(im,jsta_2l:jend_2u),bufsoil(im,nsoil,jsta_2l:jend_2u)   &
        ,buf3d(im,jm,lm),buf3d2(im,jm,lp1),buf3dx(im,lm,jm)
!jw
      integer ii,jj,js,je,jev,iyear,imn,iday,itmp,ioutcount,istatus,   &
              nsrfc,nrdlw,nrdsw,nheat,nclod,                           &
              iunit,nrecs,I,J,L

      character*80        :: titlestring

      type(r_info), pointer :: r(:) => NULL()

!
      DATA BLANK/'    '/
!
!***********************************************************************
!     START INIT HERE.
!
      WRITE(6,*)'INITPOST:  ENTER INITPOST'
!     
!     
!     STEP 1.  READ MODEL OUTPUT FILE
!
!***
! LMH always = LM for sigma-type vert coord
! LMV always = LM for sigma-type vert coord

       do j = jsta_2l, jend_2u
        do i = 1, im
            LMV ( i, j ) = lm
            LMH ( i, j ) = lm
        end do
       end do

! HTM VTM all 1 for sigma-type vert coord

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTM ( i, j, l ) = 1.0
            VTM ( i, j, l ) = 1.0
        end do
       end do
      end do

!  The end j row is going to be jend_2u for all variables except for V.
      JS=JSTA_2L
      JE=JEND_2U
      IF (JEND_2U.EQ.JM) THEN
       JEV=JEND_2U+1
      ELSE
       JEV=JEND_2U
      ENDIF

      ! Get an index of the file
      call io_int_index(filename, r, ierr)
      if (ierr /= 0) then
       print*,'Error obtinaing index of: ', trim(filename)
       stop
      end if

      ! MPI Open the file
      call mpi_file_open(mpi_comm_world, filename,                      &
                         mpi_mode_rdonly, mpi_info_null, iunit, ierr)
      if (ierr /= 0) then
       print*,"Error opening file with mpi io"
       stop
      end if

      call io_int_fetch_data(iunit, r, 'TITLE', titlestring, ierr)

!  OK, since all of the variables are dimensioned/allocated to be
!  the same size, this means we have to be careful int getVariable
!  to not try to get too much data.  For example, 
!  DUM3D is dimensioned IM+1,JM+1,LM+1 but there might actually
!  only be im,jm,lm points of data available for a particular variable.  
! get metadata

      call io_int_fetch_data(iunit, r, 'MP_PHYSICS', imp_physics, ierr)
      if (ierr /= 0) then
         imp_physics=5        ! assume ferrier if nothing specified
      endif
      if(imp_physics==85) imp_physics=5  ! HWRF scheme = Ferrier scheme
      print*,'MP_PHYSICS= ',imp_physics

      call io_int_fetch_data(iunit, r,'CU_PHYSICS', icu_physics, ierr)
      if (ierr /= 0) then
         icu_physics=4        ! assume SAS if nothing specified
      endif
      if(icu_physics==84) icu_physics=4  ! HWRF SAS = SAS
      print*,'CU_PHYSICS= ',icu_physics

      call io_int_fetch_data(iunit,r,'SF_SURFACE_PHYSICS',isf_physics,ierr)
      print*,'SF_PHYSICS= ',isf_physics

      call io_int_fetch_data(iunit, r, 'START_DATE', startdate, ierr)
      if (ierr /= 0) then
        print*,"Error reading START_DATE using MPIIO"
      else
        print*,'START_DATE from MPIIO READ= ', startdate
      end if

      jdate=0
      idate=0
      read(startdate,15)iyear,imn,iday,ihrst,imin       
 15   format(i4,1x,i2,1x,i2,1x,i2,1x,i2)
      print*,'start yr mo day hr min =',iyear,imn,iday,ihrst,imin
      print*,'processing yr mo day hr min='                             &
         ,idat(3),idat(1),idat(2),idat(4),idat(5)
      idate(1)=iyear
      idate(2)=imn
      idate(3)=iday
      idate(5)=ihrst
      idate(6)=imin
      SDAT(1)=imn
      SDAT(2)=iday
      SDAT(3)=iyear
      jdate(1)=idat(3)
      jdate(2)=idat(1)
      jdate(3)=idat(2)
      jdate(5)=idat(4)
      jdate(6)=idat(5)

      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
      ifhr=nint(rinc(2)+rinc(1)*24.)
      ifmin=nint(rinc(3))
      print*,' in INITPOST ifhr ifmin fileName=',ifhr,ifmin,fileName

! Getting tstart
      tstart=0.
      call io_int_fetch_data(iunit, r, 'TSTART', tstart, ierr)
      print*,'tstart= ',tstart

      IF(tstart .GT. 1.0E-2)THEN
       ifhr=ifhr+NINT(tstart)
       rinc=0
       idate=0
       rinc(2)=-1.0*ifhr
       call w3movdat(rinc,jdate,idate)
       SDAT(1)=idate(2)
       SDAT(2)=idate(3)
       SDAT(3)=idate(1)
       IHRST=idate(5)       
       print*,'new forecast hours for restrt run= ',ifhr
       print*,'new start yr mo day hr min =',sdat(3),sdat(1)               &
             ,sdat(2),ihrst,imin
      END IF

      RESTRT=.TRUE.  ! set RESTRT as default

      call io_int_fetch_data(iunit, r, 'DX', garb, ierr)
      print*,'DX from MPIIO READ= ',garb
      dxval=nint(garb*1000.) ! E-grid dlamda in degree
      write(6,*) 'dxval= ', dxval

      call io_int_fetch_data(iunit, r, 'DY', garb, ierr)
      print*,'DY from MPIIO READ= ',garb
      dyval=nint(garb*1000.) ! E-grid dlamda in degree
      write(6,*) 'dyval= ', dyval

      call io_int_fetch_data(iunit, r, 'DT', dt, ierr)
      write(6,*) 'DT= ', DT

      call io_int_fetch_data(iunit, r, 'CEN_LAT', garb, ierr)
      print*,'CEN_LAT from MPIIO READ= ',garb
      cenlat=nint(garb*1000.)
      write(6,*) 'cenlat= ', cenlat

      call io_int_fetch_data(iunit, r, 'CEN_LON', garb, ierr)
      print*,'CEN_LON from MPIIO READ= ',garb
      cenlon=nint(garb*1000.)
      write(6,*) 'cenlon= ', cenlon

      ! Does HWRF (NMM) use TRUELAT1 and TRUELAT2?
!      call io_int_fetch_data(iunit, r, 'TRUELAT1', garb, ierr)
!      print*,'TRUELAT1 from MPIIO READ= ',garb
!      TRUELAT1=nint(garb*1000.)
!      write(6,*) 'truelat1= ', TRUELAT1
!
!      call io_int_fetch_data(iunit, r, 'TRUELAT2', garb, ierr)
!      print*,'TRUELAT2 from MPIIO READ= ',garb
!      TRUELAT2=nint(garb*1000.)
!      write(6,*) 'truelat2= ', TRUELAT2

      call io_int_fetch_data(iunit, r, 'MAP_PROJ', maptype, ierr)
      write(6,*) 'maptype is ', maptype

      call io_int_fetch_data(iunit, r, 'GRIDTYPE', gridtype, ierr)
      write(6,*) 'gridtype is ', gridtype

      call io_int_fetch_data(iunit, r, 'HBM2', hbm2, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading HBM2: Assigned missing values"
          HBM2=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SM', sm, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SM: Assigned missing values"
          sm=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SICE', sice, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SICE: Assigned missing values"
          sice=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'PD', pd, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PD: Assigned missing values"
          pd=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'FIS', fis, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading FIS: Assigned missing values"
          fis=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'T', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading T: Assigned missing values"
          t=SPVAL
      else
          do l = 1, lm
              ll=lm-l+1
              do j = jsta_2l, jend_2u
                  do i = 1, im
                      T ( i, j, l ) = buf3d ( i, j, ll )
                  end do
              end do
          end do
      endif

      call io_int_fetch_data(iunit, r, 'Q', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading Q: Assigned missing values"
          q=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  Q ( i, j, l ) = buf3d ( i, j, ll )
              end do
          end do
      end do
      ii=im/2
      jj=(jsta+jend)/2

      call io_int_fetch_data(iunit, r, 'U', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading U: Assigned missing values"
          u=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  U ( i, j, l ) = buf3d ( i, j, ll )
                  UH( i, j, l ) = U( i, j, l )
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'V', vwbuf3dsice, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading V: Assigned missing values"
          v=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  V ( i, j, l ) = buf3d ( i, j, ll )
                  VH( i, j, l ) = V( i, j, l )
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'DX_NMM', dx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading DX_NMM: Assigned missing values"
          dx=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ETA1', eta1, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ETA1: Assigned missing values"
          eta1=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ETA2', eta2, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ETA2: Assigned missing values"
          eta2=SPVAL
      endif

      open(75,file='ETAPROFILE.txt',form='formatted',status='unknown')
      DO L=1,lm+1
        IF(L .EQ. 1)THEN
          write(75,1020)L, 0., 0.
        ELSE
          write(75,1020)L, ETA1(lm+2-l), ETA2(lm+2-l)
        END IF
      END DO
1020   format(I3,2E17.10)
      close (75)

      call io_int_fetch_data(iunit, r, 'PDTOP', pdtop, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PDTOP: Assigned missing values"
          pdtop=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'PT', pt, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PT: Assigned missing values"
          pt=SPVAL
      endif

      print*,'PT, PDTOP= ',PT,PDTOP

      call io_int_fetch_data(iunit, r, 'PBLH', pblh, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PBLH: Assigned missing values"
          pblh=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'MIXHT', mixht, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading MIXHT: Assigned missing values"
          mixht=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'USTAR', ustar, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SICE: Assigned missing values"
          sice=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'Z0', z0, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading Z0: Assigned missing values"
          sice=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'THS', ths, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading THS: Assigned missing values"
          ths=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'QS', qs, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading QS: Assigned missing values"
          qs=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'TWBS', twbs, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TWBS: Assigned missing values"
          twbs=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'QWBS', qwbs, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading QWBS: Assigned missing values"
          qwbs=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'PREC', prec, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PREC: Assigned missing values"
          prec=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ACPREC', acprec, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ACPREC: Assigned missing values"
          acprec=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CUPREC', cuprec, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CUPREC: Assigned missing values"
          cuprec=SPVAL
      endif
      do j = jsta_2l, jend_2u
        do i = 1, im
          ANCPRC(I,J)=ACPREC(I,J)-CUPREC(I,J)
        enddo
      enddo

      call io_int_fetch_data(iunit, r, 'LSPA', lspa, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading LSPA: Assigned missing values"
          lspa=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SNO', sno, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SNO: Assigned missing values"
          sno=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SI', si, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SI: Assigned missing values"
          si=SPVAL
      endif
      call io_int_fetch_data(iunit, r, 'CLDEFI', sice, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CLDEFI: Assigned missing values"
          cldefi=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'TH10', th10, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TH10: Assigned missing values"
          th10=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'Q10', q10, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading Q10: Assigned missing values"
          q10=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'PSHLTR', pshltr, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PSHLTR: Assigned missing values"
          pshltr=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'TSHLTR', tshltr, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TSHLTR: Assigned missing values"
          tshltr=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'QSHLTR', qshltr, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading QSHLTR: Assigned missing values"
          qshltr=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'Q2', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading Q2: Assigned missing values"
          q2=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  Q2 ( i, j, l ) = buf3d ( i, j, ll )
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'AKHS_OUT', akhs, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading AKHS_OUT: Assigned missing values"
          akhs=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'AKMS_OUT', akms, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading AKMS_OUT: Assigned missing values"
          akms=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ALBASE', albase, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ALBASE: Assigned missing values"
          albase=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ALBEDO', albedo, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ALBEDO: Assigned missing values"
          albedo=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CZEN', czen, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SICE: Assigned missing values"
          sice=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CZMEAN', czmean, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CZMEAN: Assigned missing values"
          czmean=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'GLAT', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading GLAT: Assigned missing values"
          sice=SPVAL
      endif
      do j = jsta_2l, jend_2u
        do i = 1, im
          F(I,J)=1.454441e-4*sin(buf(I,J))   ! 2*omeg*sin(phi)
          GDLAT(I,J)=buf(I,J)*RTD
        enddo
      enddo

      call io_int_fetch_data(iunit, r, 'GLON', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SICE: Assigned missing values"
          gdlon=SPVAL
      endif
      do j = jsta_2l, jend_2u
        do i = 1, im
          GDLON(I,J)=buf(I,J)*RTD
        enddo
      enddo

      call io_int_fetch_data(iunit, r, 'MXSNAL', mxsnal, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading MXSNAL: Assigned missing values"
          mxsnal=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'RADOT', radot, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RADOT: Assigned missing values"
          radot=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SIGT4', sigt4, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SIGT4: Assigned missing values"
          sigt4=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'TGROUND', tg, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TGROUND: Assigned missing values"
          tg=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CWM', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CWM: Assigned missing values"
          cwm=SPVAL
      endif
      do l = 1, lm
        ll=lm-l+1
        do j = jsta_2l, jend_2u
          do i = 1, im
            CWM( i, j, l ) = buf3d ( i, j, ll )
          end do
        end do
      end do

      call io_int_fetch_data(iunit, r, 'F_ICE', buf3dx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading F_ICE: Assigned missing values"
          f_ice=SPVAL
      endif
      do l = 1, lm
        ll=lm-l+1
        do j = jsta_2l, jend_2u
          do i = 1, im
            F_ice( i, j, l ) = buf3dx ( i, ll, j )
          end do
         end do
      end do

      call io_int_fetch_data(iunit, r, 'F_RAIN', buf3dx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading F_RAIN: Assigned missing values"
          f_rain=SPVAL
      endif
      do l = 1, lm
        ll=lm-l+1
        do j = jsta_2l, jend_2u
          do i = 1, im
            F_rain( i, j, l ) = buf3dx ( i, ll, j )
          end do
         end do
      end do

      call io_int_fetch_data(iunit, r, 'F_RIMEF', buf3dx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading F_RIMEF: Assigned missing values"
          f_rimef=SPVAL
      endif
      do l = 1, lm
        ll=lm-l+1
        do j = jsta_2l, jend_2u
          do i = 1, im
            F_rimef( i, j, l ) = buf3dx ( i, ll, j )
          end do
         end do
      end do

      call io_int_fetch_data(iunit, r, 'CLDFRA', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CLDFRA: Assigned missing values"
          cfr=SPVAL
      endif
      do l = 1, lm
        ll=lm-l+1
        do j = jsta_2l, jend_2u
          do i = 1, im
            CFR( i, j, l ) = buf3d ( i, j, ll )
          end do
        end do
      end do

      call io_int_fetch_data(iunit, r, 'SR', sr, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SR: Assigned missing values"
          sr=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CFRACH', cfrach, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CFRACH: Assigned missing values"
          cfrach=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CFRACL', cfracl, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CFRACL: Assigned missing values"
          cfracl=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CFRACM', cfracm, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CFRACM: Assigned missing values"
          cfracm=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ISLOPE', islope, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ISLOPE: Assigned missing values"
          islope=nint(SPVAL)
      endif

! either assign SLDPTH to be the same as eta (which is original
! setup in WRF LSM) or extract thickness of soil layers from wrf
! output

! assign SLDPTH to be the same as eta
! jkw comment out because Pleim Xiu only has 2 layers
! jkw         SLDPTH(1)=0.10
! jkw         SLDPTH(2)=0.3
! jkw         SLDPTH(3)=0.6
! jkw         SLDPTH(4)=1.0
! Initialize soil depth to some bogus value
! to alert user if not found in wrfout file
       do I=1,NSOIL
        SLDPTH(I) = 0.0
       end do

      if (isf_PHYSICS == 3) then
! get SLDPTH from wrf output
        call io_int_fetch_data(iunit, r, 'SLDPTH', SLDPTH2, ierr)
        if (ierr .ne. 0) then
            print*,"Error reading ISLOPE: Assigned missing values"
            SLDPTH2=SPVAL
        endif

        DUMCST=0.0
        DO N=1,NSOIL
          DUMCST=DUMCST+SLDPTH2(N)
        END DO
        IF(ABS(DUMCST-0.).GT.1.0E-2)THEN
          DO N=1,NSOIL
            SLLEVEL(N)=SLDPTH2(N)
          END DO
        END IF
!        print*,'SLLEVEL ',(SLLEVEL(N),N=1,NSOIL)

      else ! isf_PHYSICS /= 3
        call io_int_fetch_data(iunit, r, 'DZSOIL', SLDPTH2, ierr)
        if (ierr .ne. 0) then
            print*,"Error reading ISLOPE: Assigned missing values"
            SLDPTH2=SPVAL
        endif

        DUMCST=0.0
        DO N=1,NSOIL
          DUMCST=DUMCST+SLDPTH2(N)
        END DO
        IF(ABS(DUMCST-0.).GT.1.0E-2)THEN
          DO N=1,NSOIL
            SLDPTH(N)=SLDPTH2(N)
          END DO
        END IF
!        print*,'SLDPTH= ',(SLDPTH(N),N=1,NSOIL)
      end if   ! if (isf_PHYSICS==3)

      call io_int_fetch_data(iunit, r, 'CMC', cmc, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CMC: Assigned missing values"
          cmc=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'GRNFLX', grnflx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading GRNFLX: Assigned missing values"
          grnflx=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'PCTSNO', pctsno, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PCTSNO: Assigned missing values"
          pctsno=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SOILTB', soiltb, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SOILTB: Assigned missing values"
          soiltb=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'VEGFRC', vegfrc, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading VEGFRC: Assigned missing values"
          vegfrc=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SH2O', bufsoil, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SH2O: Assigned missing values"
          sh20=SPVAL
      endif
      do l = 1, nsoil
          do j = jsta_2l, jend_2u
              do i = 1, im
                  SH2O(I,J,L)=bufSOIL(I,L,J)
              enddo
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'SMC', bufsoil, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SMC: Assigned missing values"
          smc=SPVAL
      endif
      DO L = 1, NSOIL
          DO J = JSTA_2L, JEND_2U
              DO I = 1, IM
                  SMC(I,J,L)=BUFSOIL(I,L,J)
              ENDDO
          ENDDO
      ENDDO

      call io_int_fetch_data(iunit, r, 'STC', bufsoil, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading STC: Assigned missing values"
          stc=SPVAL
      endif
      DO L = 1, NSOIL
          DO J = JSTA_2L, JEND_2U
              DO I = 1, IM
                  STC(I,J,L)=BUFSOIL(I,L,J)
              ENDDO
          ENDDO
      ENDDO

      call io_int_fetch_data(iunit, r, 'PINT', buf3d2, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading PINT: Assigned missing values"
          pint=SPVAL
      endif
      DO L = 1, lp1
          ll=lp1-l+1
          DO J = JSTA_2L, JEND_2U
              DO I = 1, IM
                  PINT(I,J,L) = buf3d2(I,J,LL)
                  ALPINT(I,J,L)=ALOG(PINT(I,J,L))
              ENDDO
          ENDDO
      ENDDO

      do l = 2, lm+1
       do j = jsta_2l, jend_2u
        do i = 1, im
!         PMID ( i, j, l-1 ) = EXP((ALOG(PINT(I,J,L-1))+
!     &               ALOG(PINT(I,J,L)))*0.5)
         PMID ( i, j, l-1 ) = (PINT(I,J,L-1)+                              &
                     PINT(I,J,L))*0.5 ! representative of what model does
        end do
       end do
      end do
!      write(0,*)' after PMID'

      do l = 1, lm
       do j = jsta, jend
        do i = 1, im-MOD(J,2) 
	 IF(J .EQ. 1 .AND. I .LT. IM)THEN   !SOUTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
         ELSE IF(J.EQ.JM .AND. I.LT.IM)THEN   !NORTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
         ELSE IF(I .EQ. 1 .AND. MOD(J,2) .EQ. 0) THEN   !WESTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))
	 ELSE IF(I .EQ. IM .AND. MOD(J,2) .EQ. 0                             &  
      	 .AND. J .LT. JM) THEN   !EASTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))  
         ELSE IF (MOD(J,2) .LT. 1) THEN
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I-1,J,L)                       &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
         ELSE
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I+1,J,L)                       &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
         END IF  
        end do
       end do
      end do
      write(0,*)' after PMIDV'


!!!!! COMPUTE Z
       do j = jsta_2l, jend_2u
        do i = 1, im
            ZINT(I,J,LM+1)=FIS(I,J)/G
	if (I .eq. 1 .and. J .eq. jsta_2l) then
                   write(6,*) 'G,ZINT: ', G,ZINT(I,J,LM+1)
	endif
            FI(I,J,1)=FIS(I,J)
        end do
       end do
      write(0,*)' after FI'

! SECOND, INTEGRATE HEIGHT HYDROSTATICLY
      DO L=LM,1,-1
       do j = jsta_2l, jend_2u
        do i = 1, im
         FI(I,J,2)=HTM(I,J,L)*T(I,J,L)*(Q(I,J,L)*D608+1.0)*RD*                &
                   (ALPINT(I,J,L+1)-ALPINT(I,J,L))+FI(I,J,1)
         ZINT(I,J,L)=FI(I,J,2)/G
      if(i==1.and.j==250.and.l<=2)then
        write(0,*)' zint=',zint(i,j,l),' fi=',fi(i,j,2)                       &
      ,           fi(i,j,1)
        write(0,*)' t=',t(i,j,l),' q=',q(i,j,l)
        write(0,*)' alpint=',alpint(i,j,l+1),alpint(i,j,l)
      endif
         if(i.eq.ii.and.j.eq.jj)                                              &
        print*,'L,sample HTM,T,Q,ALPINT(L+1),ALPINT(l),ZINT= '                &
        ,l,HTM(I,J,L),T(I,J,L),Q(I,J,L),ALPINT(I,J,L+1),                      &
        ALPINT(I,J,L),ZINT(I,J,L)
         FI(I,J,1)=FI(I,J,2)
        ENDDO
       ENDDO
      END DO
      print*,'finish deriving geopotential in nmm'
      write(0,*)' after ZINT lm=',lm,' js=',js,' je=',je,' im=',im
      write(0,*)' zmid lbounds=',lbound(zmid),' ubounds=',ubound(zmid)
      write(0,*)' zint lbounds=',lbound(zint),' ubounds=',ubound(zint)
      write(0,*)' pmid lbounds=',lbound(pmid),' ubounds=',ubound(pmid)
      write(0,*)' pint lbounds=',lbound(pint),' ubounds=',ubound(pint)
!
      DO L=1,LM
!      write(0,*)' zmid l=',l
       DO J=JS,JE
!      write(0,*)' zmid j=',j
        DO I=1,IM
!      write(0,*)' zmid i=',i
!         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
!      write(0,*)' pmid=',pmid(i,j,l)
!      write(0,*)' pint=',pint(i,j,l),pint(i,j,l+1)
!      write(0,*)' zint=',zint(i,j,l),zint(i,j,l+1)
         FACT=(ALOG(PMID(I,J,L))-ALOG(PINT(I,J,L)))/                      &
               (ALOG(PINT(I,J,L+1))-ALOG(PINT(I,J,L)))	 
         ZMID(I,J,L)=ZINT(I,J,L)+(ZINT(I,J,L+1)-ZINT(I,J,L))*FACT
        ENDDO
       ENDDO
      ENDDO

      call io_int_fetch_data(iunit, r, 'W', buf3d2, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading W: Assigned missing values"
          WH=SPVAL
      endif
      DO L = 1, lm
          ll=lm-l+1
          DO J = JSTA_2L, JEND_2U
              DO I = 1, IM
                  WH(I,J,L) = buf3d2(I,J,LL)
              ENDDO
          ENDDO
      ENDDO

      call io_int_fetch_data(iunit, r, 'ACFRCV', acfrcv, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ACFRCV: Assigned missing values"
          acfrcv=SPVAL
      endif
      write(6,*) 'MAX ACFRCV: ', maxval(ACFRCV)


      call io_int_fetch_data(iunit, r, 'ACFRST', acfrst, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ACFRST: Assigned missing values"
          acfrst=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SSROFF', ssroff, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SSROFF: Assigned missing values"
          ssroff=SPVAL
      endif

! reading UNDERGROUND RUNOFF
      call io_int_fetch_data(iunit, r, 'BGROFF', bgroff, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading BGROFF: Assigned missing values"
          bgroff=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'RLWIN', rlwin, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RLWIN: Assigned missing values"
          rlwin=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'RLWTOA', rlwtoa, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RLWTOA: Assigned missing values"
          rlwtoa=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ALWIN', alwin, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ALWIN: Assigned missing values"
          alwin=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ALWOUT', alwout, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ALWOUT: Assigned missing values"
          alwout=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ALWTOA', alwtoa, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ALWTOA: Assigned missing values"
          alwtoa=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'RSWIN', rswin, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RSWIN: Assigned missing values"
          rswin=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'RSWINC', rswinc, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RSWINC: Assigned missing values"
          rswinc=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'RSWOUT', rswout, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RSWOUT: Assigned missing values"
          rswout=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ASWIN', aswin, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ASWIN: Assigned missing values"
          aswin=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ASWOUT', aswout, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ASWOUT: Assigned missing values"
          aswout=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ASWTOA', aswtoa, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ASWTOA: Assigned missing values"
          aswtoa=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SFCSHX', sfcshx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SFCSHX: Assigned missing values"
          sfcshx=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SFCLHX', sfclhx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SFCLHX: Assigned missing values"
          sfclhx=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SUBSHX', subshx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SUBSHX: Assigned missing values"
          subshx=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SNOPCX', snopcx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SNOPCX: Assigned missing values"
          snopcx=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SFCUVX', sfcuvx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SFCUVX: Assigned missing values"
          sfcuvx=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'POTEVP', potevp, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading POTEVP: Assigned missing values"
          potevp=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'RLWTT', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RLWTT: Assigned missing values"
          rlwtt=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  RLWTT( i, j, l ) = buf3d ( i, j, ll )
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'RSWTT', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading RSWTT: Assigned missing values"
          rswtt=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  RSWTT( i, j, l ) = buf3d ( i, j, ll )
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'TCUCN', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TCUCN: Assigned missing values"
          tcucn=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  TCUCN( i, j, l ) = buf3d ( i, j, ll )
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'TRAIN', buf3d, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TRAIN: Assigned missing values"
          train=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  TRAIN( i, j, l ) = buf3d ( i, j, ll )
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'NCFRCV', ibuf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NCFRCV: Assigned missing values"
          ncfrcv=SPVAL
      endif
      ! Is this needed?....
      do j = jsta_2l, jend_2u
          do i = 1, im
              NCFRCV(I,J)=FLOAT(ibuf(I,J))
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'NCFRST', ibuf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NCFRST: Assigned missing values"
          ncfrst=SPVAL
      endif
      ! Is this needed?....
      do j = jsta_2l, jend_2u
          do i = 1, im
              ncfrst(I,J)=FLOAT(ibuf(I,J))
          enddo
      enddo

! set default to not empty buket
      NSRFC=0
      NRDLW=0
      NRDSW=0
      NHEAT=0
      NCLOD=0
      NPREC=0

      call io_int_fetch_data(iunit, r, 'NPHS0', nphs, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NPHS0: Assigned missing values"
          nphs=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'NPREC', nprec, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NPREC: Assigned missing values"
          nprec=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'NCLOD', nclod, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NCLOD: Assigned missing values"
          nclod=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'NHEAT', nheat, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NHEAT: Assigned missing values"
          nheat=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'NRDLW', nrdlw, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NRDLW: Assigned missing values"
          nrdlw=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'NRDSW', nrdsw, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NRDSW: Assigned missing values"
          nrdsw=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'NSRFC', nsrfc, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading NSRFC: Assigned missing values"
          nsrfc=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'AVRAIN', avrain, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading AVRAIN: Assigned missing values"
          avrain=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'AVCNVC', avcnvc, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading AVCNVC: Assigned missing values"
          avcnvc=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ARDLW', ardlw, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ARDLW: Assigned missing values"
          ardlw=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ARDSW', ardsw, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ARDSW: Assigned missing values"
          ardsw=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ASRFC', asrfc, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ASRFC: Assigned missing values"
          asrfc=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'APHTIM', aphtim, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading APHTIM: Assigned missing values"
          aphtim=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'U10', u10, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading U10: Assigned missing values"
          u10=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'V10', v10, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading V10: Assigned missing values"
          v10=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SMSTAV', smstav, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SMSTAV: Assigned missing values"
          smstav=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SMSTOT', smstot, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SMSTOT: Assigned missing values"
          smstot=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'IVGTYP', ivgtyp, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading IVGTYP: Assigned missing values"
          ivgtyp=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ISLTYP', isltyp, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ISLTYP: Assigned missing values"
          isltyp=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SFCEVP', sfcevp, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SFCEVP: Assigned missing values"
          sfcevp=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SFCEXC', sfcexc, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SFCEXC: Assigned missing values"
          sfcexc=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ACSNOW', acsnow, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ACSNOW: Assigned missing values"
          acsnow=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'ACSNOM', acsnom, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading ACSNOM: Assigned missing values"
          acsnom=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'SST', sst, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading SST: Assigned missing values"
          sst=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'TAUX', mdltaux, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TAUX: Assigned missing values"
          mdltaux=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'TAUY', mdltauy, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading TAUY: Assigned missing values"
          mdltauy=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'EL_PBL', buf3dx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading EL_PBL: Assigned missing values"
          el_pbl=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  EL_PBL( i, j, l ) = buf3dx ( i, j ,ll)
              end do
          end do
      end do

      call io_int_fetch_data(iunit, r, 'EXCH_H', buf3dx, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading EXCH_H: Assigned missing values"
          EXCH_H=SPVAL
      endif
      do l = 1, lm
          ll=lm-l+1
          do j = jsta_2l, jend_2u
              do i = 1, im
                  EXCH_H( i, j, l ) = buf3dx ( i, j ,ll)
              end do
          end do
      end do 

      call io_int_fetch_data(iunit, r, 'THZ0', thz0, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading THZ0: Assigned missing values"
          thz0=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'QZ0', qz0, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading QZ0: Assigned missing values"
          qz0=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'UZ0', uz0, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading UZ0: Assigned missing values"
          uz0=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'VZ0', vz0, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading VZ0: Assigned missing values"
          vz0=SPVAL
      endif

!
! Very confusing story ...
!
! Retrieve htop and hbot => They are named CNVTOP, CNVBOT in the model and
! with HBOTS,HTOPS (shallow conv) and HBOTD,HTOPD (deep conv) represent
! the 3 sets of convective cloud base/top arrays tied to the frequency
! that history files are written.
!
! IN THE *MODEL*, arrays HBOT,HTOP are similar to CNVTOP,CNVBOT but are
! used in radiation and are tied to the frequency of radiation updates.
!
! For historical reasons model arrays CNVTOP,CNVBOT are renamed HBOT,HTOP
! and manipulated throughout the post. 

      call io_int_fetch_data(iunit, r, 'CNVTOP', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CNVTOP: Assigned missing values"
          htop=SPVAL
      endif
      do j = jsta_2l, jend_2u
          do i = 1, im
              HTOP ( i, j ) = float(LM)-buf(i,j)+1.0
              HTOP ( i, j ) = max(1.0,min(HTOP(I,J),float(LM)))
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'CNVBOT', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CNVBOT: Assigned missing values"
          hbot=SPVAL
      endif
      do j = jsta_2l, jend_2u
          do i = 1, im
              HBOT ( i, j ) = float(LM)-buf(i,j)+1.0
              HBOT ( i, j ) = max(1.0,min(HBOT(I,J),float(LM)))
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'HTOPD', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading HTOPD: Assigned missing values"
          htopd=SPVAL
      endif
      do j = jsta_2l, jend_2u
          do i = 1, im
              HTOPD ( i, j ) = float(LM)-buf(i,j)+1.0
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'HBOTD', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading HBOTD: Assigned missing values"
          hbotd=SPVAL
      endif
      do j = jsta_2l, jend_2u
          do i = 1, im
              HBOTD ( i, j ) = float(LM)-buf(i,j)+1.0
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'HTOPS', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading HTOPS: Assigned missing values"
          htops=SPVAL
      endif
      do j = jsta_2l, jend_2u
          do i = 1, im
              HTOPS ( i, j ) = float(LM)-buf(i,j)+1.0
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'HBOTS', buf, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading HBOTS: Assigned missing values"
          hbots=SPVAL
      endif
      do j = jsta_2l, jend_2u
          do i = 1, im
              HBOTS ( i, j ) = float(LM)-buf(i,j)+1.0
          enddo
      enddo

      call io_int_fetch_data(iunit, r, 'CUPPT', cuppt, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CUPPT: Assigned missing values"
          cuppt=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'CPRATE', cprate, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading CPRATE: Assigned missing values"
          cprate=SPVAL
      endif

      call io_int_fetch_data(iunit, r, 'HBM2', hbm2, ierr)
      if (ierr .ne. 0) then
          print*,"Error reading HBM2: Assigned missing values"
          hbm2=SPVAL
      endif

!!!! DONE GETTING

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            IF(ABS(T(I,J,L)).GT.1.0E-3)                                &
              OMGA(I,J,L) = -WH(I,J,L)*PMID(I,J,L)*G/                   &
                       (RD*T(I,J,L)*(1.+D608*Q(I,J,L)))

        end do
       end do
      end do
      write(0,*)' after OMGA'

! pos east
      call collect_loc(gdlat,dummy)
      if(me.eq.0)then
        latstart=nint(dummy(1,1)*1000.)
        latlast=nint(dummy(im,jm)*1000.)
! temporary patch for nmm wrf for moving nest. gopal's doing
! jkw changed if statement as per MP's suggestion
! jkw        if(mod(im,2).ne.0) then
! chuang: test
        icen=(im+1)/2
        jcen=(jm+1)/2
        if(mod(im,2).ne.0)then !per Pyle, jm is always odd
         if(mod(jm+1,4).ne.0)then
          cenlat=nint(dummy(icen,jcen)*1000.)
         else
          cenlat=nint(0.5*(dummy(icen-1,jcen)+dummy(icen,jcen))*1000.)
         end if
        else
         if(mod(jm+1,4).ne.0)then
          cenlat=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
         else
          cenlat=nint(dummy(icen,jcen)*1000.)
         end if
        end if

!        if(mod(im,2).eq.0) then
!           icen=(im+1)/2
!           jcen=(jm+1)/2
!           cenlat=nint(dummy(icen,jcen)*1000.)
!        else
!           icen=im/2
!           jcen=(jm+1)/2
!           cenlat=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
!        end if
        
      end if
      write(6,*) 'laststart,latlast,cenlat B calling bcast= ', &
     &    latstart,latlast,cenlat
      call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
      call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
      call mpi_bcast(cenlat,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
      write(6,*) 'laststart,latlast,cenlat A calling bcast= ', &
     &    latstart,latlast,cenlat

      call collect_loc(gdlon,dummy)
      if(me.eq.0)then
        lonstart=nint(dummy(1,1)*1000.)
        lonlast=nint(dummy(im,jm)*1000.)
! temporary patch for nmm wrf for moving nest. gopal's doing
!lrb changed if statement as per MP's suggestion
!lrb        if(mod(im,2).ne.0) then
!Chuang: test
        icen=(im+1)/2
        jcen=(jm+1)/2
        if(mod(im,2).ne.0)then !per Pyle, jm is always odd
         if(mod(jm+1,4).ne.0)then
          cenlon=nint(dummy(icen,jcen)*1000.)
         else
          cenlon=nint(0.5*(dummy(icen-1,jcen)+dummy(icen,jcen))*1000.)
         end if
        else
         if(mod(jm+1,4).ne.0)then
          cenlon=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
         else
          cenlon=nint(dummy(icen,jcen)*1000.)
         end if
        end if

!        if(mod(im,2).eq.0) then
!           icen=(im+1)/2
!           jcen=(jm+1)/2
!           cenlon=nint(dummy(icen,jcen)*1000.)
!        else
!           icen=im/2
!           jcen=(jm+1)/2
!           cenlon=nint(0.5*(dummy(icen,jcen)+dummy(icen+1,jcen))*1000.)
!        end if
       end if

       write(6,*)'lonstart,lonlast,cenlon B calling bcast= ', &
     &      lonstart,lonlast,cenlon
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       call mpi_bcast(cenlon,1,MPI_INTEGER,0,mpi_comm_comp,irtn)
       write(6,*)'lonstart,lonlast,cenlon A calling bcast= ', &
     &      lonstart,lonlast,cenlon
!
        write(6,*) 'filename in INITPOST=', filename


!MEB not sure how to get these 
       do j = jsta_2l, jend_2u
        do i = 1, im
!            DX ( i, j ) = dxval  !MEB ???
!            DY ( i, j ) = dyval*DTR*ERAD  

!!!!!!!!!!!!!!!!!!!!! DY ????

            DY ( i, j ) =   0.001*ERAD*DYVAL*DTR  ! like A*DPH
        end do
       end do
!MEB not sure how to get these 
      write(0,*)' after DX DY'

! close up shop
!      call ext_int_ioclose ( DataHandle, Status )

! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,                                       &
                RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

      CALL TABLEQ(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)
      write(0,*)' after TABLEQ'


!     
!     
      IF(ME.EQ.0)THEN
        WRITE(6,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        WRITE(6,51) (SPL(L),L=1,LSM)
   50   FORMAT(14(F4.1,1X))
   51   FORMAT(8(F8.1,1X))
      ENDIF
!     
!     COMPUTE DERIVED TIME STEPPING CONSTANTS.
!
!MEB need to get DT
!      DT = 120. !MEB need to get DT
!      NPHS = 4  !MEB need to get physics DT
      DTQ2 = DT * NPHS  !MEB need to get physics DT
      TSPH = 3600./DT   !MEB need to get DT

      TSRFC=float(NSRFC)/TSPH
      IF(NSRFC.EQ.0)TSRFC=float(ifhr)  !in case buket does not get emptied
      TRDLW=float(NRDLW)/TSPH
      IF(NRDLW.EQ.0)TRDLW=float(ifhr)  !in case buket does not get emptied
      TRDSW=float(NRDSW)/TSPH
      IF(NRDSW.EQ.0)TRDSW=float(ifhr)  !in case buket does not get emptied
      THEAT=float(NHEAT)/TSPH
      IF(NHEAT.EQ.0)THEAT=float(ifhr)  !in case buket does not get emptied
      TCLOD=float(NCLOD)/TSPH
      IF(NCLOD.EQ.0)TCLOD=float(ifhr)  !in case buket does not get emptied
      TPREC=float(NPREC)/TSPH
      IF(NPREC.EQ.0)TPREC=float(ifhr)  !in case buket does not get emptied
!       TPREC=float(ifhr)
      print*,'TSRFC TRDLW TRDSW= ',TSRFC, TRDLW, TRDSW
!MEB need to get DT

!how am i going to get this information?
!      NPREC  = INT(TPREC *TSPH+D50)
!      NHEAT  = INT(THEAT *TSPH+D50)
!      NCLOD  = INT(TCLOD *TSPH+D50)
!      NRDSW  = INT(TRDSW *TSPH+D50)
!      NRDLW  = INT(TRDLW *TSPH+D50)
!      NSRFC  = INT(TSRFC *TSPH+D50)
!how am i going to get this information?
!     
!     IF(ME.EQ.0)THEN
!       WRITE(6,*)' '
!       WRITE(6,*)'DERIVED TIME STEPPING CONSTANTS'
!       WRITE(6,*)' NPREC,NHEAT,NSRFC :  ',NPREC,NHEAT,NSRFC
!       WRITE(6,*)' NCLOD,NRDSW,NRDLW :  ',NCLOD,NRDSW,NRDLW
!     ENDIF
!
!     COMPUTE DERIVED MAP OUTPUT CONSTANTS.
      DO L = 1,LSM
         ALSL(L) = ALOG(SPL(L))
      END DO
      write(0,*)' after ALSL'
!
!HC WRITE IGDS OUT FOR WEIGHTMAKER TO READ IN AS KGDSIN
        if(me.eq.0)then
        print*,'writing out igds'
        igdout=110
!        open(igdout,file='griddef.out',form='unformatted'
!     +  ,status='unknown')
        if(maptype .eq. 1)THEN  ! Lambert conformal
          WRITE(igdout)3
          WRITE(6,*)'igd(1)=',3
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        ELSE IF(MAPTYPE .EQ. 2)THEN  !Polar stereographic
          WRITE(igdout)5
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2  !Assume projection at +-90
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        ELSE IF(MAPTYPE .EQ. 3)THEN  !Mercator
          WRITE(igdout)1
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)latlast
          WRITE(igdout)lonlast
          WRITE(igdout)TRUELAT1
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)255
        ELSE IF(MAPTYPE.EQ.0 .OR. MAPTYPE.EQ.203)THEN  !A STAGGERED E-GRID
          WRITE(igdout)203
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)136
          WRITE(igdout)CENLAT
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)64
          WRITE(igdout)0
          WRITE(igdout)0
          WRITE(igdout)0
        END IF
        end if
      write(0,*)' after writes'

        call mpi_file_close(iunit,ierr)

      RETURN
      END
