C-----------------------------------------------------------------------
      SUBROUTINE POLATES0(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
     &                    NO,RLAT,RLON,IBO,LO,GO,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  POLATES0   INTERPOLATE SCALAR FIELDS (BILINEAR)
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
C
C ABSTRACT: THIS SUBPROGRAM PERFORMS BILINEAR INTERPOLATION
C           FROM ANY GRID TO ANY GRID FOR SCALAR FIELDS.
C           OPTIONS ALLOW VARYING THE MINIMUM PERCENTAGE FOR MASK,
C           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
C           (IPOPT(1)) WHICH DEFAULTS TO 50 (IF IPOPT(1)=-1).
C           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
C           IF NO INPUT DATA IS FOUND NEAR THE OUTPUT POINT, A SPIRAL
C           SEARCH MAY BE INVOKED BY SETTING IPOPT(2)> 0.
C           NO SEARCHING IS DONE IF OUTPUT POINT IS OUTSIDE THE INPUT GRID.
C           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
C           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
C           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
C             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
C             (KGDS(1)=001) MERCATOR CYLINDRICAL
C             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
C             (KGDS(1)=004) GAUSSIAN CYLINDRICAL (SPECTRAL NATIVE)
C             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
C             (KGDS(1)=202) ROTATED EQUIDISTANT CYLINDRICAL (ETA NATIVE)
C           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
C           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
C           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
C           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
C           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
C           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
C           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
C           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
C           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
C           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
C        
C PROGRAM HISTORY LOG:
C   96-04-10  IREDELL
C 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
C 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
C 2007-05-22  IREDELL  EXTRAPOLATE UP TO HALF A GRID CELL
C 2008-06-04  GAYNO    ADDED SPIRAL SEARCH OPTION
C 2009-10-19  IREDELL  SAVE WEIGHTS AND THREAD FOR PERFORMANCE
C 2012-06-26  GAYNO    FIX OUT-OF-BOUNDS ERROR.  SEE NCEPLIBS
C                      TICKET #9.
C
C USAGE:    CALL POLATES0(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,
C    &                    NO,RLAT,RLON,IBO,LO,GO,IRET)
C
C   INPUT ARGUMENT LIST:
C     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
C                IPOPT(1) IS MINIMUM PERCENTAGE FOR MASK
C                (DEFAULTS TO 50 IF IPOPT(1)=-1)
C                IPOPT(2) IS WIDTH OF SQUARE TO EXAMINE IN SPIRAL SEARCH
C                (DEFAULTS TO NO SEARCH IF IPOPT(2)=-1)
C     KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
C     KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
C                (KGDSO(1)<0 IMPLIES RANDOM STATION POINTS)
C     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
C                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
C     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
C                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
C     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
C     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
C     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
C     GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
C     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)<0)
C     RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
C     RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
C
C   OUTPUT ARGUMENT LIST:
C     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)>=0)
C     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)>=0)
C     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)>=0)
C     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
C     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
C     GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
C     IRET     - INTEGER RETURN CODE
C                0    SUCCESSFUL INTERPOLATION
C                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
C                3    UNRECOGNIZED OUTPUT GRID
C
C SUBPROGRAMS CALLED:
C   GDSWIZ       GRID DESCRIPTION SECTION WIZARD
C   IJKGDS0      SET UP PARAMETERS FOR IJKGDS1
C   (IJKGDS1)    RETURN FIELD POSITION FOR A GIVEN GRID POINT
C   POLFIXS      MAKE MULTIPLE POLE SCALAR VALUES CONSISTENT
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      IMPLICIT NONE
      INTEGER,INTENT(IN):: IPOPT(20),KGDSI(200),KGDSO(200),MI,MO,KM
      INTEGER,INTENT(IN):: IBI(KM)
      LOGICAL*1,INTENT(IN):: LI(MI,KM)
      REAL,INTENT(IN):: GI(MI,KM)
      INTEGER,INTENT(INOUT):: NO
      REAL,INTENT(INOUT):: RLAT(MO),RLON(MO)
      INTEGER,INTENT(OUT):: IBO(KM)
      LOGICAL*1,INTENT(OUT):: LO(MO,KM)
      REAL,INTENT(OUT):: GO(MO,KM)
      INTEGER,INTENT(OUT):: IRET
      REAL XPTS(MO),YPTS(MO)
      INTEGER IJX(2),IJY(2)
      REAL WX(2),WY(2)
      INTEGER IJKGDSA(20)
      REAL,PARAMETER:: FILL=-9999.
      INTEGER MP,N,I,J,K,NK,NV,IJKGDS1
      INTEGER MSPIRAL,I1,J1,IXS,JXS,MX,KXS,KXT,IX,JX,NX
      REAL PMP,XIJ,YIJ,XF,YF,G,W
      INTEGER,SAVE:: KGDSIX(200)=-1,KGDSOX(200)=-1,NOX=-1,IRETX=-1
      INTEGER,ALLOCATABLE,SAVE:: NXY(:,:,:)
      REAL,ALLOCATABLE,SAVE:: RLATX(:),RLONX(:),WXY(:,:,:)
      REAL,ALLOCATABLE:: CROT(:),SROT(:)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SET PARAMETERS
      IRET=0
      MP=IPOPT(1)
      IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
      IF(MP.LT.0.OR.MP.GT.100) IRET=32
      PMP=MP*0.01
      MSPIRAL=MAX(IPOPT(2),0)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SAVE OR SKIP WEIGHT COMPUTATION
      IF(IRET.EQ.0.AND.(KGDSO(1).LT.0.OR.
     &    ANY(KGDSI.NE.KGDSIX).OR.ANY(KGDSO.NE.KGDSOX))) THEN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        IF(KGDSO(1).GE.0) THEN
          ALLOCATE (CROT(MO))
          ALLOCATE (SROT(MO))
          CALL GDSWIZ(KGDSO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO,0,
     &                CROT,SROT)
          DEALLOCATE (CROT,SROT)
          IF(NO.EQ.0) IRET=3
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  LOCATE INPUT POINTS
        ALLOCATE (CROT(NO))
        ALLOCATE (SROT(NO))
        CALL GDSWIZ(KGDSI,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV,0,
     &              CROT,SROT)
        DEALLOCATE (CROT,SROT)
        IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ALLOCATE AND SAVE GRID DATA
        KGDSIX=KGDSI
        KGDSOX=KGDSO
        IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,NXY,WXY)
          ALLOCATE(RLATX(NO),RLONX(NO),NXY(2,2,NO),WXY(2,2,NO))
          NOX=NO
        ENDIF
        IRETX=IRET
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE WEIGHTS
        IF(IRET.EQ.0) THEN
          CALL IJKGDS0(KGDSI,IJKGDSA)
C$OMP PARALLEL DO
C$OMP&PRIVATE(N,XIJ,YIJ,IJX,IJY,XF,YF,J,I,WX,WY)
          DO N=1,NO
            RLONX(N)=RLON(N)
            RLATX(N)=RLAT(N)
            XIJ=XPTS(N)
            YIJ=YPTS(N)
            IF(XIJ.NE.FILL.AND.YIJ.NE.FILL) THEN
              IJX(1:2)=FLOOR(XIJ)+(/0,1/)
              IJY(1:2)=FLOOR(YIJ)+(/0,1/)
              XF=XIJ-IJX(1)
              YF=YIJ-IJY(1)
              WX(1)=(1-XF)
              WX(2)=XF
              WY(1)=(1-YF)
              WY(2)=YF
              DO J=1,2
                DO I=1,2
                  NXY(I,J,N)=IJKGDS1(IJX(I),IJY(J),IJKGDSA)
                  WXY(I,J,N)=WX(I)*WY(J)
                ENDDO
              ENDDO
            ELSE
              NXY(:,:,N)=0
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INTERPOLATE OVER ALL FIELDS
      IF(IRET.EQ.0.AND.IRETX.EQ.0) THEN
        IF(KGDSO(1).GE.0) THEN
          NO=NOX
          DO N=1,NO
            RLON(N)=RLONX(N)
            RLAT(N)=RLATX(N)
          ENDDO
        ENDIF
C$OMP PARALLEL DO
C$OMP&PRIVATE(NK,K,N,G,W,J,I)
C$OMP&PRIVATE(I1,J1,IXS,JXS,MX,KXS,KXT,IX,JX,NX)
        DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          G=0
          W=0
          DO J=1,2
            DO I=1,2
              IF(NXY(I,J,N).GT.0)THEN
                IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K)) THEN
                  G=G+WXY(I,J,N)*GI(NXY(I,J,N),K)
                  W=W+WXY(I,J,N)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          LO(N,K)=W.GE.PMP
          IF(LO(N,K)) THEN
            GO(N,K)=G/W
          ELSEIF(MSPIRAL.GT.0.AND.XPTS(N).NE.FILL.AND.
     &                            YPTS(N).NE.FILL) THEN
            I1=NINT(XPTS(N))
            J1=NINT(YPTS(N))
            IXS=SIGN(1.,XPTS(N)-I1)
            JXS=SIGN(1.,YPTS(N)-J1)
            SPIRAL : DO MX=1,MSPIRAL**2
              KXS=SQRT(4*MX-2.5)
              KXT=MX-(KXS**2/4+1)
              SELECT CASE(MOD(KXS,4))
              CASE(1)
                IX=I1-IXS*(KXS/4-KXT)
                JX=J1-JXS*KXS/4
              CASE(2)
                IX=I1+IXS*(1+KXS/4)
                JX=J1-JXS*(KXS/4-KXT)
              CASE(3)
                IX=I1+IXS*(1+KXS/4-KXT)
                JX=J1+JXS*(1+KXS/4)
              CASE DEFAULT
                IX=I1-IXS*KXS/4
                JX=J1+JXS*(KXS/4-KXT)
              END SELECT
              NX=IJKGDS1(IX,JX,IJKGDSA)
              IF(NX.GT.0.)THEN
                IF(LI(NX,K).OR.IBI(K).EQ.0)THEN
                  GO(N,K)=GI(NX,K)
                  LO(N,K)=.TRUE.
                  EXIT SPIRAL
                ENDIF
              ENDIF
            ENDDO SPIRAL
            IF(.NOT.LO(N,K))THEN
              IBO(K)=1
              GO(N,K)=0.
            ENDIF
          ELSE
            GO(N,K)=0.
          ENDIF
        ENDDO
        DO K=1,KM
          IBO(K)=IBI(K)
          IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
        ENDDO
        IF(KGDSO(1).EQ.0) CALL POLFIXS(NO,MO,KM,RLAT,RLON,IBO,LO,GO)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        IF(IRET.EQ.0) IRET=IRETX
        IF(KGDSO(1).GE.0) NO=0
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
