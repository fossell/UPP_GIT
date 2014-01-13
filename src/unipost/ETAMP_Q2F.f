      SUBROUTINE ETAMP_Q2F(QRIMEF)
      ! This subroutine is to be used with the WRF "advected Ferrier
      ! scheme" to calculate the F_ICE, F_RIMEF and F_RAIN arrays from
      ! the QQW, QQR, QQI and the input array QRIMEF.
        use CTLBLK_mod, only: lm,im,jsta,jend,jsta_2l,jend_2u
        use VRBLS3D_mod, only: QQW,QQR,QQI
        real, parameter :: t_ice=-40., t0c=273.15, t_icek=t0c+t_ice
        real, parameter :: epsq=1.e-12
        real, intent(in) :: QRIMEF(im,jsta_2l:jend_2u,lm)

      bigl: do l=1,lm
         bigj: do j=jsta,jend
            bigi: do i=1,im
               QT(I,J,K)=QQW(I,J,K)+QQR(I,J,K)+QQI(I,J,K)
               if(QQI(i,j,k)<=EPSQ) then
                  f_ice(i,j,k)=0.
                  f_rimef(i,j,k)=0.
                  if(T(i,j,k)<T_ICEK) f_ice(i,j,k)=1.
               else
                  f_ice(i,j,k)=max(o.,min(1.,QQI(i,j,k)/QT(i,j,k)))
                  f_rimef(i,j,k)=rimeq(i,j,k)/QQI(i,j,k)
               endif
               if(QQR(i,j,k) <= EPSQ) then
                  F_RAIN(i,j,k)=0.
               else
                  F_RAIN=QQR(I,J,K)/(QQR(I,J,K)+QQW(I,J,K))
               endif
            enddo bigi
         enddo bigj
      enddo bigl
      END SUBROUTINE ETAMP_Q2F
