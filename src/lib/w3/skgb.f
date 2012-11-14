C-----------------------------------------------------------------------
      SUBROUTINE SKGB(LUGB,ISEEK,MSEEK,LSKIP,LGRIB)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: SKGB           SEARCH FOR NEXT GRIB MESSAGE
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 93-11-22
C
C ABSTRACT: THIS SUBPROGRAM SEARCHES A FILE FOR THE NEXT GRIB 1 MESSAGE.
C   A GRIB 1 MESSAGE IS IDENTIFIED BY ITS INDICATOR SECTION, I.E.
C   AN 8-BYTE SEQUENCE WITH 'GRIB' IN BYTES 1-4 AND 1 IN BYTE 8.
C   IF FOUND, THE LENGTH OF THE MESSAGE IS DECODED FROM BYTES 5-7.
C   THE SEARCH IS DONE OVER A GIVEN SECTION OF THE FILE.
C   THE SEARCH IS TERMINATED IF AN EOF OR I/O ERROR IS ENCOUNTERED.
C
C PROGRAM HISTORY LOG:
C   93-11-22  IREDELL
C   95-10-31  IREDELL   ADD CALL TO BAREAD 
C   97-03-14  IREDELL   CHECK FOR '7777'
C 2001-06-05  IREDELL   APPLY LINUX PORT BY EBISUZAKI
C
C USAGE:    CALL SKGB(LUGB,ISEEK,MSEEK,LSKIP,LGRIB)
C   INPUT ARGUMENTS:
C     LUGB         INTEGER LOGICAL UNIT OF INPUT GRIB FILE
C     ISEEK        INTEGER NUMBER OF BYTES TO SKIP BEFORE SEARCH
C     MSEEK        INTEGER MAXIMUM NUMBER OF BYTES TO SEARCH
C   OUTPUT ARGUMENTS:
C     LSKIP        INTEGER NUMBER OF BYTES TO SKIP BEFORE MESSAGE
C     LGRIB        INTEGER NUMBER OF BYTES IN MESSAGE (0 IF NOT FOUND)
C
C SUBPROGRAMS CALLED:
C   BAREAD       BYTE-ADDRESSABLE READ
C   GBYTEC       GET INTEGER DATA FROM BYTES
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      PARAMETER(LSEEK=128)
      CHARACTER Z(LSEEK)
      CHARACTER Z4(4)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LGRIB=0
      KS=ISEEK
      KN=MIN(LSEEK,MSEEK)
      KZ=LSEEK
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  LOOP UNTIL GRIB MESSAGE IS FOUND
      DOWHILE(LGRIB.EQ.0.AND.KN.GE.8.AND.KZ.EQ.LSEEK)
C  READ PARTIAL SECTION
        CALL BAREAD(LUGB,KS,KN,KZ,Z)
        KM=KZ-8+1
        K=0
C  LOOK FOR 'GRIB...1' IN PARTIAL SECTION
        DOWHILE(LGRIB.EQ.0.AND.K.LT.KM)
          CALL GBYTEC(Z,I4,(K+0)*8,4*8)
          CALL GBYTEC(Z,I1,(K+7)*8,1*8)
          IF(I4.EQ.1196575042.AND.I1.EQ.1) THEN
C  LOOK FOR '7777' AT END OF GRIB MESSAGE
            CALL GBYTEC(Z,KG,(K+4)*8,3*8)
            CALL BAREAD(LUGB,KS+K+KG-4,4,K4,Z4)
            IF(K4.EQ.4) THEN
              CALL GBYTEC(Z4,I4,0,4*8)
              IF(I4.EQ.926365495) THEN
C  GRIB MESSAGE FOUND
                LSKIP=KS+K
                LGRIB=KG
              ENDIF
            ENDIF
          ENDIF
          K=K+1
        ENDDO
        KS=KS+KM
        KN=MIN(LSEEK,ISEEK+MSEEK-KS)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
