      SUBROUTINE ETAMP_Q2F(QRIMEF)
      ! This subroutine is to be used with the WRF "advected Ferrier
      ! scheme" to calculate the F_ICE, F_RIMEF and F_RAIN arrays from
      ! the QQW, QQR, QQI and the input array QRIMEF.
        use CTLBLK_mod, only: lm,im,jsta,jend,jsta_2l,jend_2u
        use VRBLS3D, only: QQW,QQR,QQI,CWM, f_rain,f_ice,f_rimef, T

        implicit none

        real, intent(in) :: QRIMEF(im,jsta_2l:jend_2u,lm)
        integer :: i,j,l
        real :: qt

        ! NOTE: these parameters must match the WRF Ferrier scheme.
        ! They're wrong elsewhere in the post:
        real, parameter :: t_ice=-40., t0c=273.15, t_icek=t0c+t_ice
        real, parameter :: epsq=1.e-12

      bigl: do l=1,lm
         bigj: do j=jsta,jend
            bigi: do i=1,im
               QT=QQW(I,J,L)+QQR(I,J,L)+QQI(I,J,L)
               CWM(i,j,l)=QT
               if(QQI(i,j,l)<=EPSQ) then
                  f_ice(i,j,l)=0.
                  f_rimef(i,j,l)=0.
                  if(T(i,j,l)<T_ICEK) f_ice(i,j,l)=1.
               else
                  f_ice(i,j,l)=max(0.,min(1.,QQI(i,j,l)/QT))
                  f_rimef(i,j,l)=QRIMEF(i,j,l)/QQI(i,j,l)
               endif
               if(QQR(i,j,l) <= EPSQ) then
                  F_RAIN(i,j,l)=0.
               else
                  F_RAIN(i,j,l)=QQR(I,J,L)/(QQR(I,J,L)+QQW(I,J,L))
               endif
            enddo bigi
         enddo bigj
      enddo bigl
      END SUBROUTINE ETAMP_Q2F
