      SUBROUTINE W3FI74 (IGDS,ICOMP,GDS,LENGDS,NPTS,IGERR)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    W3FI74      CONSTRUCT GRID DEFINITION SECTION (GDS)
C   PRGMMR: FARLEY           ORG: W/NMC42    DATE: 93-08-24
C
C ABSTRACT: THIS SUBROUTINE CONSTRUCTS A GRIB GRID DEFINITION
C   SECTION.
C
C PROGRAM HISTORY LOG:
C   92-07-07  M. FARLEY   ORIGINAL AUTHOR
C   92-10-16  R.E.JONES   ADD CODE TO LAT/LON SECTION TO DO
C                         GAUSSIAN GRIDS.
C   93-03-29  R.E.JONES   ADD SAVE STATEMENT
C   93-08-24  R.E.JONES   CHANGES FOR GRIB GRIDS 37-44
C   93-09-29  R.E.JONES   CHANGES FOR GAUSSIAN GRID FOR DOCUMENT
C                         CHANGE IN W3FI71.
C   94-02-15  R.E.JONES   CHANGES FOR ETA MODEL GRIDS 90-93
C   95-04-20  R.E.JONES   CHANGE 200 AND 201 TO 201 AND 202
C   95-10-31  IREDELL     REMOVED SAVES AND PRINTS
C   98-08-20  BALDWIN     ADD TYPE 203
C   07-03-20  VUONG       ADD TYPE 204
C   10-01-21  GAYNO       ADD GRID 205 - ROTATED LAT/LON A,B,C,D STAGGERS
C
C
C USAGE:    CALL W3FI74 (IGDS, ICOMP, GDS, LENGDS, NPTS, IGERR)
C   INPUT ARGUMENT LIST:
C     IGDS        - INTEGER ARRAY SUPPLIED BY W3FI71
C     ICOMP       - TABLE 7- RESOLUTION & COMPONENT FLAG (BIT 5)
C                   FOR GDS(17) WIND COMPONENTS
C
C   OUTPUT ARGUMENT LIST:
C     GDS       - COMPLETED GRIB GRID DEFINITION SECTION
C     LENGDS    - LENGTH OF GDS
C     NPTS      - NUMBER OF POINTS IN GRID
C     IGERR     - 1, GRID REPRESENTATION TYPE NOT VALID
C
C REMARKS: SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C ATTRIBUTES:
C   LANGUAGE: CRAY CFT77 FORTRAN 77, IBM370 VS FORTRAN
C   MACHINE:  CRAY C916-128, CRAY Y-MP8/864, CRAY Y-MP EL2/256, HDS
C
C$$$
C
      INTEGER       IGDS  (*)
      CHARACTER*1   GDS   (*)
      INTEGER       ITEMP
C
      ISUM  = 0
      IGERR = 0
C
C       PRINT *,' '
C       PRINT *,'(W3FI74-IGDS = )'
C       PRINT *,(IGDS(I),I=1,18)
C       PRINT *,' '
C
C     COMPUTE LENGTH OF GDS IN OCTETS (OCTETS 1-3)
C       LENGTH =  32 FOR LAT/LON, GNOMIC, GAUSIAN LAT/LON,
C                    POLAR STEREOGRAPHIC, SPHERICAL HARMONICS
C                    ROTATED LAT/LON E-STAGGER
C       LENGDS =  34 ROTATED LAT/LON A,B,C,D STAGGERS
C       LENGTH =  42 FOR MERCATOR, LAMBERT, TANGENT CONE
C       LENGTH = 178 FOR MERCATOR, LAMBERT, TANGENT CONE
C
      IF (IGDS(3) .EQ. 0  .OR.  IGDS(3) .EQ. 2  .OR.
     &    IGDS(3) .EQ. 4  .OR.  IGDS(3) .EQ. 5  .OR.
     &    IGDS(3) .EQ. 50 .OR.  IGDS(3) .EQ. 201.OR.
     &    IGDS(3) .EQ. 202.OR.  IGDS(3) .EQ. 203.OR.
     &    IGDS(3) .EQ. 204) THEN
          LENGDS = 32
C
C       CORRECTION FOR GRIDS 37-44
C
        IF (IGDS(3).EQ.0.AND.IGDS(1).EQ.0.AND.IGDS(2).NE.
     &  255) THEN
          LENGDS = IGDS(5) * 2 + 32
        ENDIF
      ELSE IF (IGDS(3) .EQ. 1  .OR.  IGDS(3) .EQ. 3  .OR.
     &         IGDS(3) .EQ. 13) THEN
        LENGDS = 42
      ELSE IF (IGDS(3) .EQ. 205) THEN
        LENGDS = 34
      ELSE
C       PRINT *,' W3FI74 ERROR, GRID REPRESENTATION TYPE NOT VALID'
        IGERR = 1
        RETURN
      ENDIF
C
C     PUT LENGTH OF GDS SECTION IN BYTES 1,2,3
C
      GDS(1) = CHAR(MOD(LENGDS/65536,256))
      GDS(2) = CHAR(MOD(LENGDS/  256,256))
      GDS(3) = CHAR(MOD(LENGDS      ,256))
C
C     OCTET 4 = NV, NUMBER OF VERTICAL COORDINATE PARAMETERS
C     OCTET 5 = PV, PL OR 255
C     OCTET 6 = DATA REPRESENTATION TYPE (TABLE 6)
C
      GDS(4) = CHAR(IGDS(1))
      GDS(5) = CHAR(IGDS(2))
      GDS(6) = CHAR(IGDS(3))
C
C     FILL OCTET THE REST OF THE GDS BASED ON DATA REPRESENTATION
C     TYPE (TABLE 6)
C
C$$
C     PROCESS ROTATED LAT/LON A,B,C,D STAGGERS
C
      IF (IGDS(3).EQ.205) THEN
        GDS( 7) = CHAR(MOD(IGDS(4)/256,256))
        GDS( 8) = CHAR(MOD(IGDS(4)    ,256))
        GDS( 9) = CHAR(MOD(IGDS(5)/256,256))
        GDS(10) = CHAR(MOD(IGDS(5)    ,256))
        LATO    = IGDS(6) ! LAT OF FIRST POINT
        IF (LATO .LT. 0) THEN
          LATO = -LATO
          LATO = IOR(LATO,8388608)
        ENDIF
        GDS(11) = CHAR(MOD(LATO/65536,256))
        GDS(12) = CHAR(MOD(LATO/  256,256))
        GDS(13) = CHAR(MOD(LATO      ,256))
        LONO    = IGDS(7) ! LON OF FIRST POINT
        IF (LONO .LT. 0) THEN
          LONO = -LONO
          LONO = IOR(LONO,8388608)
        ENDIF
        GDS(14) = CHAR(MOD(LONO/65536,256))
        GDS(15) = CHAR(MOD(LONO/  256,256))
        GDS(16) = CHAR(MOD(LONO      ,256))
        LATEXT  = IGDS(9)  ! CENTER LAT
        IF (LATEXT .LT. 0) THEN
          LATEXT = -LATEXT
          LATEXT = IOR(LATEXT,8388608)
        ENDIF
        GDS(18) = CHAR(MOD(LATEXT/65536,256))
        GDS(19) = CHAR(MOD(LATEXT/  256,256))
        GDS(20) = CHAR(MOD(LATEXT      ,256))
        LONEXT  = IGDS(10) ! CENTER LON
        IF (LONEXT .LT. 0) THEN
          LONEXT = -LONEXT
          LONEXT = IOR(LONEXT,8388608)
        ENDIF
        GDS(21) = CHAR(MOD(LONEXT/65536,256))
        GDS(22) = CHAR(MOD(LONEXT/  256,256))
        GDS(23) = CHAR(MOD(LONEXT      ,256))
        GDS(24) = CHAR(MOD(IGDS(11)/256,256))
        GDS(25) = CHAR(MOD(IGDS(11)    ,256))
        GDS(26) = CHAR(MOD(IGDS(12)/256,256))
        GDS(27) = CHAR(MOD(IGDS(12)    ,256))
        GDS(28) = CHAR(IGDS(13))
        LATO    = IGDS(14) ! LAT OF LAST POINT
        IF (LATO .LT. 0) THEN
          LATO = -LATO
          LATO = IOR(LATO,8388608)
        ENDIF
        GDS(29) = CHAR(MOD(LATO/65536,256))
        GDS(30) = CHAR(MOD(LATO/  256,256))
        GDS(31) = CHAR(MOD(LATO      ,256))
        LONO    = IGDS(15) ! LON OF LAST POINT
        IF (LONO .LT. 0) THEN
          LONO = -LONO
          LONO = IOR(LONO,8388608)
        ENDIF
        GDS(32) = CHAR(MOD(LONO/65536,256))
        GDS(33) = CHAR(MOD(LONO/  256,256))
        GDS(34) = CHAR(MOD(LONO      ,256))
C
C     PROCESS LAT/LON GRID TYPES OR GAUSSIAN GRID OR ARAKAWA
C     STAGGERED, SEMI-STAGGERED, OR FILLED E-GRIDS
C
      ELSEIF (IGDS(3).EQ.0.OR.IGDS(3).EQ.4.OR.
     &    IGDS(3).EQ.201.OR.IGDS(3).EQ.202.OR.
     &    IGDS(3).EQ.203.OR.IGDS(3).EQ.204) THEN
        GDS( 7) = CHAR(MOD(IGDS(4)/256,256))
        GDS( 8) = CHAR(MOD(IGDS(4)    ,256))
        GDS( 9) = CHAR(MOD(IGDS(5)/256,256))
        GDS(10) = CHAR(MOD(IGDS(5)    ,256))
        LATO    = IGDS(6)
        IF (LATO .LT. 0) THEN
          LATO = -LATO
          LATO = IOR(LATO,8388608)
        ENDIF
        GDS(11) = CHAR(MOD(LATO/65536,256))
        GDS(12) = CHAR(MOD(LATO/  256,256))
        GDS(13) = CHAR(MOD(LATO      ,256))
        LONO    = IGDS(7)
        IF (LONO .LT. 0) THEN
          LONO = -LONO
          LONO = IOR(LONO,8388608)
        ENDIF
        GDS(14) = CHAR(MOD(LONO/65536,256))
        GDS(15) = CHAR(MOD(LONO/  256,256))
        GDS(16) = CHAR(MOD(LONO      ,256))
        LATEXT  = IGDS(9)
        IF (LATEXT .LT. 0) THEN
          LATEXT = -LATEXT
          LATEXT = IOR(LATEXT,8388608)
        ENDIF
        GDS(18) = CHAR(MOD(LATEXT/65536,256))
        GDS(19) = CHAR(MOD(LATEXT/  256,256))
        GDS(20) = CHAR(MOD(LATEXT      ,256))
        LONEXT  = IGDS(10)
        IF (LONEXT .LT. 0) THEN
          LONEXT = -LONEXT
          LONEXT = IOR(LONEXT,8388608)
        ENDIF
        GDS(21) = CHAR(MOD(LONEXT/65536,256))
        GDS(22) = CHAR(MOD(LONEXT/  256,256))
        GDS(23) = CHAR(MOD(LONEXT      ,256))
        IRES    = IAND(IGDS(8),128)
        IF (IGDS(3).EQ.201.OR.IGDS(3).EQ.202.OR.
     &      IGDS(3).EQ.203.OR.IGDS(3).EQ.204) THEN
          GDS(24) = CHAR(MOD(IGDS(11)/256,256))
          GDS(25) = CHAR(MOD(IGDS(11)    ,256))
        ELSE IF (IRES.EQ.0) THEN
          GDS(24) = CHAR(255)
          GDS(25) = CHAR(255)
        ELSE
          GDS(24) = CHAR(MOD(IGDS(12)/256,256))
          GDS(25) = CHAR(MOD(IGDS(12)    ,256))
        END IF
        IF (IGDS(3).EQ.4) THEN
          GDS(26) = CHAR(MOD(IGDS(11)/256,256))
          GDS(27) = CHAR(MOD(IGDS(11)    ,256))
        ELSE IF (IGDS(3).EQ.201.OR.IGDS(3).EQ.202.OR.
     &           IGDS(3).EQ.203.OR.IGDS(3).EQ.204) THEN
          GDS(26) = CHAR(MOD(IGDS(12)/256,256))
          GDS(27) = CHAR(MOD(IGDS(12)    ,256))
        ELSE IF (IRES.EQ.0) THEN
          GDS(26) = CHAR(255)
          GDS(27) = CHAR(255)
        ELSE
          GDS(26) = CHAR(MOD(IGDS(11)/256,256))
          GDS(27) = CHAR(MOD(IGDS(11)    ,256))
        END IF
        GDS(28) = CHAR(IGDS(13))
        GDS(29) = CHAR(0)
        GDS(30) = CHAR(0)
        GDS(31) = CHAR(0)
        GDS(32) = CHAR(0)
        IF (LENGDS.GT.32) THEN
          ISUM = 0
          I    = 19
          DO 10 J = 33,LENGDS,2
            ISUM     = ISUM + IGDS(I)
            GDS(J)   = CHAR(MOD(IGDS(I)/256,256))
            GDS(J+1) = CHAR(MOD(IGDS(I)    ,256))
            I        = I + 1
 10       CONTINUE
        END IF
C
C$$     PROCESS MERCATOR GRID TYPES
C
      ELSE IF (IGDS(3) .EQ. 1) THEN
        GDS( 7) = CHAR(MOD(IGDS(4)/256,256))
        GDS( 8) = CHAR(MOD(IGDS(4)    ,256))
        GDS( 9) = CHAR(MOD(IGDS(5)/256,256))
        GDS(10) = CHAR(MOD(IGDS(5)    ,256))
        LATO = IGDS(6)
        IF (LATO .LT. 0) THEN
          LATO = -LATO
          LATO = IOR(LATO,8388608)
        ENDIF
        GDS(11) = CHAR(MOD(LATO/65536,256))
        GDS(12) = CHAR(MOD(LATO/  256,256))
        GDS(13) = CHAR(MOD(LATO      ,256))
        LONO = IGDS(7)
        IF (LONO .LT. 0) THEN
          LONO = -LONO
          LONO = IOR(LONO,8388608)
        ENDIF
        GDS(14) = CHAR(MOD(LONO/65536,256))
        GDS(15) = CHAR(MOD(LONO/  256,256))
        GDS(16) = CHAR(MOD(LONO      ,256))
        LATEXT = IGDS(9)
        IF (LATEXT .LT. 0) THEN
          LATEXT = -LATEXT
          LATEXT = IOR(LATEXT,8388608)
        ENDIF
        GDS(18) = CHAR(MOD(LATEXT/65536,256))
        GDS(19) = CHAR(MOD(LATEXT/  256,256))
        GDS(20) = CHAR(MOD(LATEXT      ,256))
        LONEXT  = IGDS(10)
        IF (LONEXT .LT. 0) THEN
          LONEXT = -LONEXT
          LONEXT = IOR(LONEXT,8388608)
        ENDIF
        GDS(21) = CHAR(MOD(LONEXT/65536,256))
        GDS(22) = CHAR(MOD(LONEXT/  256,256))
        GDS(23) = CHAR(MOD(LONEXT      ,256))
        GDS(24) = CHAR(MOD(IGDS(13)/65536,256))
        GDS(25) = CHAR(MOD(IGDS(13)/  256,256))
        GDS(26) = CHAR(MOD(IGDS(13)      ,256))
        GDS(27) = CHAR(0)
        GDS(28) = CHAR(IGDS(14))
        GDS(29) = CHAR(MOD(IGDS(12)/65536,256))
        GDS(30) = CHAR(MOD(IGDS(12)/  256,256))
        GDS(31) = CHAR(MOD(IGDS(12)      ,256))
        GDS(32) = CHAR(MOD(IGDS(11)/65536,256))
        GDS(33) = CHAR(MOD(IGDS(11)/  256,256))
        GDS(34) = CHAR(MOD(IGDS(11)      ,256))
        GDS(35) = CHAR(0)
        GDS(36) = CHAR(0)
        GDS(37) = CHAR(0)
        GDS(38) = CHAR(0)
        GDS(39) = CHAR(0)
        GDS(40) = CHAR(0)
        GDS(41) = CHAR(0)
        GDS(42) = CHAR(0)
C$$     PROCESS LAMBERT CONFORMAL GRID TYPES

      ELSE IF (IGDS(3) .EQ. 3) THEN
        GDS( 7) = CHAR(MOD(IGDS(4)/256,256))
        GDS( 8) = CHAR(MOD(IGDS(4)    ,256))
        GDS( 9) = CHAR(MOD(IGDS(5)/256,256))
        GDS(10) = CHAR(MOD(IGDS(5)    ,256))
        LATO = IGDS(6)
        IF (LATO .LT. 0) THEN
          LATO = -LATO
          LATO = IOR(LATO,8388608)
        ENDIF
        GDS(11) = CHAR(MOD(LATO/65536,256))
        GDS(12) = CHAR(MOD(LATO/  256,256))
        GDS(13) = CHAR(MOD(LATO      ,256))
        LONO = IGDS(7)
        IF (LONO .LT. 0) THEN
          LONO = -LONO
          LONO = IOR(LONO,8388608)
        ENDIF
        GDS(14) = CHAR(MOD(LONO/65536,256))
        GDS(15) = CHAR(MOD(LONO/  256,256))
        GDS(16) = CHAR(MOD(LONO      ,256))
        LONM = IGDS(9)
        IF (LONM .LT. 0) THEN
          LONM = -LONM
          LONM = IOR(LONM,8388608)
        ENDIF
        GDS(18) = CHAR(MOD(LONM/65536,256))
        GDS(19) = CHAR(MOD(LONM/  256,256))
        GDS(20) = CHAR(MOD(LONM      ,256))
        GDS(21) = CHAR(MOD(IGDS(10)/65536,256))
        GDS(22) = CHAR(MOD(IGDS(10)/  256,256))
        GDS(23) = CHAR(MOD(IGDS(10)      ,256))
        GDS(24) = CHAR(MOD(IGDS(11)/65536,256))
        GDS(25) = CHAR(MOD(IGDS(11)/  256,256))
        GDS(26) = CHAR(MOD(IGDS(11)      ,256))
        GDS(27) = CHAR(IGDS(12))
        GDS(28) = CHAR(IGDS(13))

! Modified by Robert Rozumalski to account for negative
! Lat/Lon values in IGDS[15] through IGDS[18]
!
        LAT2 = IGDS(15)
        IF (LAT2 .LT. 0) THEN
          LAT2 = -LAT2
          LAT2 = IOR(LAT2,8388608)
        ENDIF
        GDS(29) = CHAR(MOD(LAT2/65536,256))
        GDS(30) = CHAR(MOD(LAT2/  256,256))
        GDS(31) = CHAR(MOD(LAT2      ,256))
!       GDS(29) = CHAR(MOD(IGDS(15)/65536,256))
!       GDS(30) = CHAR(MOD(IGDS(15)/  256,256))
!       GDS(31) = CHAR(MOD(IGDS(15)      ,256))
        LAT1 = IGDS(16)
        IF (LAT1 .LT. 0) THEN
          LAT1 = -LAT1
          LAT1 = IOR(LAT1,8388608)
        ENDIF
        GDS(32) = CHAR(MOD(LAT1/65536,256))
        GDS(33) = CHAR(MOD(LAT1/  256,256))
        GDS(34) = CHAR(MOD(LAT1      ,256))
!       GDS(32) = CHAR(MOD(IGDS(16)/65536,256))
!       GDS(33) = CHAR(MOD(IGDS(16)/  256,256))
!       GDS(34) = CHAR(MOD(IGDS(16)      ,256))
        LATP = IGDS(17)
        IF (LATP .LT. 0) THEN
          LATP = -LATP
          LATP = IOR(LATP,8388608)
        ENDIF
        GDS(35) = CHAR(MOD(LATP/65536,256))
        GDS(36) = CHAR(MOD(LATP/  256,256))
        GDS(37) = CHAR(MOD(LATP      ,256))
!       GDS(35) = CHAR(MOD(IGDS(17)/65536,256))
!       GDS(36) = CHAR(MOD(IGDS(17)/  256,256))
!       GDS(37) = CHAR(MOD(IGDS(17)      ,256))
        LONP = IGDS(18)
        IF (LONP .LT. 0) THEN
          LONP = -LONP
          LONP = IOR(LONP,8388608)
        ENDIF
        GDS(38) = CHAR(MOD(LONP/65536,256))
        GDS(39) = CHAR(MOD(LONP/  256,256))
        GDS(40) = CHAR(MOD(LONP      ,256))
!       GDS(38) = CHAR(MOD(IGDS(18)/65536,256))
!       GDS(39) = CHAR(MOD(IGDS(18)/  256,256))
!       GDS(40) = CHAR(MOD(IGDS(18)      ,256))

        GDS(41) = CHAR(0)
        GDS(42) = CHAR(0)

C$$     PROCESS POLAR STEREOGRAPHIC GRID TYPES
      ELSE IF (IGDS(3) .EQ. 5) THEN
        GDS( 7) = CHAR(MOD(IGDS(4)/256,256))
        GDS( 8) = CHAR(MOD(IGDS(4)    ,256))
        GDS( 9) = CHAR(MOD(IGDS(5)/256,256))
        GDS(10) = CHAR(MOD(IGDS(5)    ,256))
        LATO = IGDS(6)
        IF (LATO .LT. 0) THEN
          LATO = -LATO
          LATO = IOR(LATO,8388608)
        ENDIF
        GDS(11) = CHAR(MOD(LATO/65536,256))
        GDS(12) = CHAR(MOD(LATO/  256,256))
        GDS(13) = CHAR(MOD(LATO      ,256))
        LONO = IGDS(7)
        IF (LONO .LT. 0) THEN
          LONO = -LONO
          LONO = IOR(LONO,8388608)
        ENDIF
        GDS(14) = CHAR(MOD(LONO/65536,256))
        GDS(15) = CHAR(MOD(LONO/  256,256))
        GDS(16) = CHAR(MOD(LONO      ,256))
        LONM = IGDS(9)
        IF (LONM .LT. 0) THEN
          LONM = -LONM
          LONM = IOR(LONM,8388608)
        ENDIF
        GDS(18) = CHAR(MOD(LONM/65536,256))
        GDS(19) = CHAR(MOD(LONM/   256,256))
        GDS(20) = CHAR(MOD(LONM       ,256))
        GDS(21) = CHAR(MOD(IGDS(10)/65536,256))
        GDS(22) = CHAR(MOD(IGDS(10)/  256,256))
        GDS(23) = CHAR(MOD(IGDS(10)      ,256))
        GDS(24) = CHAR(MOD(IGDS(11)/65536,256))
        GDS(25) = CHAR(MOD(IGDS(11)/  256,256))
        GDS(26) = CHAR(MOD(IGDS(11)      ,256))
        GDS(27) = CHAR(IGDS(12))
        GDS(28) = CHAR(IGDS(13))
        GDS(29) = CHAR(0)
        GDS(30) = CHAR(0)
        GDS(31) = CHAR(0)
        GDS(32) = CHAR(0)
      ENDIF
C       PRINT 10,(GDS(IG),IG=1,32)
C10     FORMAT ('  GDS= ',32(1X,Z2.2))
C
C     COMPUTE NUMBER OF POINTS IN GRID BY MULTIPLYING
C       IGDS(4) AND IGDS(5) ... NEEDED FOR PACKER
C
      IF (IGDS(3).EQ.0.AND.IGDS(1).EQ.0.AND.IGDS(2).NE.
     & 255) THEN
        NPTS = ISUM
      ELSE
        NPTS = IGDS(4) * IGDS(5)
      ENDIF
C
C     'IOR' ICOMP-BIT 5 RESOLUTION & COMPONENT FLAG FOR WINDS
C       WITH IGDS(8) INFO (REST OF RESOLUTION & COMPONENT FLAG DATA)
C       ITEMP used to avoid setting an input parameter
C
      ITEMP   = ISHFT(ICOMP,3)
      GDS(17) = CHAR(IOR(IGDS(8),ITEMP))
C
      RETURN
      END
