       SUBROUTINE CALRAD
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALRAD      
!   PRGRMMR: CHUANG        ORG: EMC      DATE: 07-01-17       
!     
! ABSTRACT:
!     THIS ROUTINE COMPUTES MODEL DERIVED BRIGHTNESS TEMPERATURE
!     USING CRTM. IT IS PATTERNED AFTER GSI SETUPRAD WITH TREADON'S HELP     
! PROGRAM HISTORY LOG:
!
! USAGE:    CALL MDLFLD
!   INPUT ARGUMENT LIST:
!     NONE
!   OUTPUT ARGUMENT LIST: 
!     NONE
!
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!
!     LIBRARY:
!     /nwprod/lib/sorc/crtm2
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : IBM
!$$$  
       use vrbls3d
       use vrbls2d
       use masks
       use soil
      
       use kinds, only: r_kind,r_single,i_kind
       use crtm_module
       use crtm_parameters, only: limit_exp,toa_pressure,max_n_layers
       use error_handler, only: success

       use params_mod
       use ctlblk_mod
       use rqstfld_mod

      implicit none
!     
!     DECLARE VARIABLES.
!     
! Mapping land surface type of GFS to CRTM
!  Note: index 0 is water, and index 13 is ice. The two indices are not
!        used and just assigned to COMPACTED_SOIL.
        integer, parameter, dimension(0:13) :: gfs_to_crtm=(/COMPACTED_SOIL,    &
      &  BROADLEAF_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, PINE_FOREST,&  
      &  PINE_FOREST, BROADLEAF_BRUSH, SCRUB, SCRUB, SCRUB_SOIL, TUNDRA,        &
      &  COMPACTED_SOIL, TILLED_SOIL, COMPACTED_SOIL/)
!       integer gfs_to_crtm(0:13)
!       data gfs_to_crtm /COMPACTED_SOIL, 
!     &  BROADLEAF_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, 
!     &  PINE_FOREST, PINE_FOREST, BROADLEAF_BRUSH, SCRUB, SCRUB, 
!     &  SCRUB_SOIL, TUNDRA, 
!     &  COMPACTED_SOIL, TILLED_SOIL, COMPACTED_SOIL/

! Mapping land surface type of NMM to CRTM
!  Note: index 16 is water, and index 24 is ice. The two indices are not
!        used and just assigned to COMPACTED_SOIL.
       integer, parameter, dimension(24) :: nmm_to_crtm=(/URBAN_CONCRETE,     &
      &   COMPACTED_SOIL, IRRIGATED_LOW_VEGETATION, GRASS_SOIL, MEADOW_GRASS, &
      &   MEADOW_GRASS, MEADOW_GRASS, SCRUB, GRASS_SCRUB, MEADOW_GRASS,       &
      &   BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST,       &
      &   BROADLEAF_PINE_FOREST, COMPACTED_SOIL, WET_SOIL, WET_SOIL,          &
      &   IRRIGATED_LOW_VEGETATION, TUNDRA, TUNDRA, TUNDRA, TUNDRA,           &
      &   COMPACTED_SOIL/)

!      integer nmm_to_crtm(24)
!      data nmm_to_crtm /URBAN_CONCRETE, COMPACTED_SOIL, 
!     &   IRRIGATED_LOW_VEGETATION, GRASS_SOIL, MEADOW_GRASS, 
!     &   MEADOW_GRASS, MEADOW_GRASS, SCRUB, GRASS_SCRUB, MEADOW_GRASS, 
!     &   BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, PINE_FOREST, 
!     &   BROADLEAF_PINE_FOREST, COMPACTED_SOIL, WET_SOIL, WET_SOIL, 
!     &   IRRIGATED_LOW_VEGETATION, TUNDRA, TUNDRA, TUNDRA, TUNDRA, 
!     &   COMPACTED_SOIL/ 
       
      integer, parameter:: ndat=100
! CRTM structure variable declarations.
      integer,parameter::  n_absorbers = 2
      integer,parameter::  n_clouds = 0 ! may add cloud in the future
      integer,parameter::  n_aerosols = 0
      integer(i_kind) lunin,nobs,nchanl,nreal
      integer(i_kind) error_status,itype
      integer(i_kind) err1,err2,err3,err4
      integer(i_kind) i,j,k,msig
      integer jdn
  
      real(r_kind),parameter:: r100=100.0_r_kind
      real,parameter:: ozsmall = 1.e-10 ! to convert to mass mixing ratio
      real(r_kind) tsfc 
      real(r_kind),dimension(4):: sfcpct
      real(r_kind) snodepth,snoeqv,vegcover
      real snofrac
      real(r_kind),dimension(im,jsta:jend):: tb1,tb2,tb3,tb4
      real,dimension(im,jm):: grid1
      real sun_zenith,sun_azimuth
  
      character(10)::obstype
      character(20)::isis
  
      logical hirs2,msu,goessndr,hirs3,hirs4,hirs,amsua,amsub,airs,hsb  &
     &      ,goes_img,mhs
      logical avhrr,avhrr_navy,lextra,ssu
      logical ssmi,ssmis,amsre,amsre_low,amsre_mid,amsre_hig,change
      logical ssmis_las,ssmis_uas,ssmis_env,ssmis_img
      logical sea,mixed,land,ice,snow,toss
      logical micrim,microwave
!  logical,dimension(nobs):: luse
    
      type(crtm_atmosphere_type) :: atmosphere
      type(crtm_surface_type) :: surface
      type(crtm_geometryinfo_type) :: geometryinfo
      type(crtm_channelinfo_type):: channelinfo
      type(crtm_rtsolution_type),allocatable,dimension(:):: rtsolution     

      integer ii,jj
      integer, external :: iw3jdn
!     
!*****************************************************************************
!     START SUBROUTINE CALRAD.
      if (iget(327) > 0 .or. iget(328) > 0 .or. iget(329) > 0     &
     & .or. iget(330) > 0 ) then
! Initialize ozone to zeros for WRF NMM and ARW for now
       if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR')o3=0.0
! Compute solar zenith angle for GFS
       if (MODELNAME == 'GFS')then
        jdn=iw3jdn(idat(3),idat(1),idat(2))
	do j=jsta,jend
	 do i=1,im
	  call zensun(jdn,float(idat(4)),gdlat(i,j),gdlon(i,j)     &
     &	    ,pi,sun_zenith,sun_azimuth)
          czen(i,j)=cos(sun_zenith)
	 end do
	end do
!	ii=361
!	jj=278
	ii=1035
	jj=219
	if(jj>=jsta .and. jj<=jend)                                 &
     &    print*,'sample GFS zenith angle=',acos(czen(ii,jj))*rtd   
       end if	       
! Initialize CRTM.  Load satellite sensor array.
! The optional arguments Process_ID and Output_Process_ID limit
! generation of runtime informative output to mpi task
! Output_Process_ID (which here is set to be task 0)
       print*,'success in CALRAD= ',success
       error_status = crtm_init(channelinfo, Process_ID=me,         &
     &      Output_Process_ID=0)
       if (error_status /= 0_i_kind)                                &
     &   write(6,*)'ERROR*** crtm_init error_status=',error_status
  
!   lunin=1 ! will read data file in the future, only simulate GOES for now
!   open(lunin,file='obs_setup',form='unformatted') ! still need to find out filename
!   rewind lunin
  
! Loop over data types to process    
!  do is=1,ndat
!    read(lunin,end=125) obstype,isis,nreal,nchanl
! Test GOES for now.
       obstype='goes_img'
       isis='imgr_g12'
!       isis='amsua_n15'

! Initialize logical flags for satellite platform

       hirs2      = obstype == 'hirs2'
       hirs3      = obstype == 'hirs3'
       hirs4      = obstype == 'hirs4'
       hirs       = hirs2 .or. hirs3 .or. hirs4
       msu        = obstype == 'msu'
       ssu        = obstype == 'ssu'
       goessndr   = obstype == 'sndr'  .or. obstype == 'sndrd1' .or.   &
     &                obstype == 'sndrd2'.or. obstype == 'sndrd3' .or. &
     &                obstype == 'sndrd4'
       amsua      = obstype == 'amsua'
       amsub      = obstype == 'amsub'
       mhs        = obstype == 'mhs'
       airs       = obstype == 'airs'
       hsb        = obstype == 'hsb'
       goes_img   = obstype == 'goes_img'
       avhrr      = obstype == 'avhrr'
       avhrr_navy = obstype == 'avhrr_navy'
       ssmi       = obstype == 'ssmi'
       amsre_low  = obstype == 'amsre_low'
       amsre_mid  = obstype == 'amsre_mid'
       amsre_hig  = obstype == 'amsre_hig'
       amsre      = amsre_low .or. amsre_mid .or. amsre_hig
       ssmis      = obstype == 'ssmis'
       ssmis_las  = obstype == 'ssmis_las'
       ssmis_uas  = obstype == 'ssmis_uas'
       ssmis_img  = obstype == 'ssmis_img'
       ssmis_env  = obstype == 'ssmis_env'

       ssmis=ssmis_las.or.ssmis_uas.or.ssmis_img.or.ssmis_env

       micrim=ssmi .or. ssmis .or. amsre   ! only used for MW-imager-QC and id_qc(ch)

       microwave=amsua .or. amsub  .or. mhs .or. msu .or. hsb .or.   &
     &             micrim 

! Set CRTM to process given satellite/sensor
       Error_Status = CRTM_Set_ChannelInfo(isis, ChannelInfo)
       if (error_status /= success)                                  &
     &    write(6,*)'ERROR*** crtm_set_channelinfo error_status=',   &
     &    error_status,' for satsensor=',isis

! Allocate structures for radiative transfer
       print*,'channel number= ',channelinfo%n_channels
       allocate(rtsolution  (channelinfo%n_channels))
       err1=0; err2=0; err3=0; err4=0
       if(lm > max_n_layers)then
        write(6,*) 'CALRAD: lm > max_n_layers - '//                  &
     &	  'increase crtm max_n_layers ',lm,max_n_layers
        stop 2
       end if
       err1 = crtm_allocate_atmosphere(lm,n_absorbers,n_clouds       &
     &       ,n_aerosols,atmosphere)
       err2 = crtm_allocate_surface(channelinfo%n_channels,surface)
       err3 = crtm_allocate_rtsolution(lm,rtsolution)
       if (err1/=success) write(6,*)' ERROR** allocating atmosphere. &
     &       err=',err1
       if (err2/=success) write(6,*)' ***ERROR** allocating surface. &
     &       err=',err2
       if (err3/=success) write(6,*)' ***ERROR** allocating rtsolution. &
     &       err=',err3

       atmosphere%n_layers = lm
       atmosphere%level_temperature_input = 0
       atmosphere%absorber_id(1) = H2O_ID
       atmosphere%absorber_id(2) = O3_ID
       atmosphere%absorber_units(1) = MASS_MIXING_RATIO_UNITS
       atmosphere%absorber_units(2) = MASS_MIXING_RATIO_UNITS ! Use mass mixing ratio
       atmosphere%level_pressure(0) = TOA_PRESSURE

!    if(nchanl /= channelinfo%n_channels) write(6,*)'***ERROR** nchanl,n_channels ', &
!           nchanl,channelinfo%n_channels

! Load surface sensor data structure
       surface%sensordata%n_channels = channelinfo%n_channels
       surface%sensordata%sensor_id  = channelinfo%wmo_sensor_id(1)   

! Loop through all grid points
!      ii=361
!      jj=278
      do j=jsta,jend
       do i=1,im

!    Load geometry structure
!    geometryinfo%sensor_zenith_angle = zasat*rtd  ! local zenith angle ???????
        geometryinfo%sensor_zenith_angle = 0. ! 44.
        geometryinfo%source_zenith_angle = acos(czen(i,j))*rtd ! solar zenith angle
!        geometryinfo%source_zenith_angle = 0. 
        geometryinfo%sensor_scan_angle   = 0. ! scan angle, assuming nadir
        if(abs(geometryinfo%sensor_zenith_angle) > h1 ) then
          geometryinfo%distance_ratio =                       &
     &        abs( sin(geometryinfo%sensor_scan_angle*dtr)/   &
     &        sin(geometryinfo%sensor_zenith_angle*dtr) )
        endif
        if(i==ii.and.j==jj)print*,'sample geometry ',         &
     &  geometryinfo%sensor_zenith_angle                      &
     &  ,geometryinfo%source_zenith_angle                     &
     &  ,czen(i,j)*rtd 
!  Set land/sea, snow, ice percentages and flags
        if (MODELNAME == 'GFS')then ! GFS uses 13 veg types
         itype=IVGTYP(I,J)
         itype = min(max(0,ivgtyp(i,j)),13)
!         IF(itype <= 0 .or. itype > 13)itype=7 !use scrub for ocean point
         if(sno(i,j)/=spval)then
          snoeqv=sno(i,j)
         else
          snoeqv=0.
         end if
         if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp B cing snfrc = ', &
     &      snoeqv,itype,IVGTYP(I,J)
         if(sm(i,j) > 0.1)then
          sfcpct(4)=0.
         else 
          call snfrac_gfs(i,j,SNOeqv,itype,sfcpct(4))
         end if
         if(i==ii.and.j==jj)print*,'sno,itype,ivgtyp,sfcpct(4) = ',     &
     &      snoeqv,itype,IVGTYP(I,J),sfcpct(4)
        else          
         itype=IVGTYP(I,J)
         IF(itype == 0)itype=8
         CALL SNFRAC (SNO(I,J),itype,sfcpct(4))
        end if 
!	CALL SNFRAC (SNO(I,J),IVGTYP(I,J),snofrac)
!	sfcpct(4)=snofrac
	if(sm(i,j) > 0.1)then ! water
!	 tsfc=sst(i,j)
         tsfc = ths(i,j)*(pint(i,j,nint(lmh(i,j))+1)/p1000)**capa
         vegcover=0.0
	 if(sfcpct(4) > 0.0_r_kind)then ! snow and water
          sfcpct(1) = 1.0_r_kind-sfcpct(4)
          sfcpct(2) = 0.0_r_kind
	  sfcpct(3) = 0.0_r_kind
	 else ! pure water
	  sfcpct(1) = 1.0_r_kind
          sfcpct(2) = 0.0_r_kind
	  sfcpct(3) = 0.0_r_kind
	 end if  
        else ! land and sea ice
	 tsfc = ths(i,j)*(pint(i,j,nint(lmh(i,j))+1)/p1000)**capa
         vegcover=vegfrc(i,j)
	 if(sice(i,j) > 0.1)then ! sea ice
	  if(sfcpct(4) > 0.0_r_kind)then ! sea ice and snow
	   sfcpct(3) = 1.0_r_kind-sfcpct(4)
	   sfcpct(1) = 0.0_r_kind
           sfcpct(2) = 0.0_r_kind
	  else ! pure sea ice
	   sfcpct(3)= 1.0_r_kind
	   sfcpct(1) = 0.0_r_kind
           sfcpct(2) = 0.0_r_kind
	  end if
	 else ! land
	  if(sfcpct(4) > 0.0_r_kind)then ! land and snow
	   sfcpct(2)= 1.0_r_kind-sfcpct(4)
           sfcpct(1) = 0.0_r_kind
           sfcpct(3) = 0.0_r_kind
	  else ! pure land
	   sfcpct(2)= 1.0_r_kind
           sfcpct(1) = 0.0_r_kind
           sfcpct(3) = 0.0_r_kind
	  end if  
         end if
	end if 
        if(si(i,j)/=spval)then 	
	 snodepth = si(i,j)
        else
         snodepth = 0.
        end if

        sea  = sfcpct(1)  >= 0.99_r_kind
        land = sfcpct(2)  >= 0.99_r_kind
        ice  = sfcpct(3)  >= 0.99_r_kind
        snow = sfcpct(4)  >= 0.99_r_kind
        mixed = .not. sea  .and. .not. ice .and.          &
     &             .not. land .and. .not. snow
        if((sfcpct(1)+sfcpct(2)+sfcpct(3)+sfcpct(4))      &
     &    >1._r_kind)print*,'ERROR sfcpct ',i,j,sfcpct(1) &
     &    ,sfcpct(2),sfcpct(3),sfcpct(4)
!    Load surface structure

!    Define land characteristics

!    **NOTE:  The model surface type --> CRTM surface type
!             mapping below is specific to the versions NCEP
!             GFS and NNM as of September 2005
!    itype = ivgtyp(i,j)
       if (MODELNAME == 'NMM' .OR. MODELNAME == 'NCAR' .OR. MODELNAME == 'RAPR') then
        itype = min(max(1,ivgtyp(i,j)),24)
        surface%land_type = nmm_to_crtm(itype)
       else
        itype = min(max(0,ivgtyp(i,j)),13)
        surface%land_type = gfs_to_crtm(itype)
       end if

       surface%wind_speed            = sqrt(u10(i,j)*u10(i,j) &
     &                                   +v10(i,j)*v10(i,j))    
       surface%water_coverage        = sfcpct(1)
       surface%land_coverage         = sfcpct(2)
       surface%ice_coverage          = sfcpct(3)
       surface%snow_coverage         = sfcpct(4)
       
       surface%land_temperature      = tsfc
       surface%snow_temperature      = min(tsfc,280._r_kind)
       surface%water_temperature     = max(tsfc,270._r_kind)
       surface%ice_temperature       = min(tsfc,280._r_kind)
       surface%soil_moisture_content = smc(i,j,1) !???
       surface%vegetation_fraction   = vegcover
!       surface%vegetation_fraction   = vegfrc(i,j)
       surface%soil_temperature      = 283.
!       surface%soil_temperature      = stc(i,j,1)
       surface%snow_depth            = snodepth ! in mm
       if(i==ii.and.j==jj)print*,'sample surface in CALRAD=',  &
     &   i,j,surface%wind_speed,surface%water_coverage,        &
     &   surface%land_coverage,surface%ice_coverage,           &
     &   surface%snow_coverage,surface%land_temperature,       &
     &   surface%snow_temperature,surface%water_temperature,   &
     &   surface%ice_temperature,surface%vegetation_fraction,  &
     &   surface%soil_temperature,surface%snow_depth,          &
     &   surface%land_type,sm(i,j)

!       Load profiles into model layers

!       Load atmosphere profiles into RTM model layers
!       CRTM counts from top down just as post does
       if(i==ii.and.j==jj)print*,'TOA= ',atmosphere%level_pressure(0)
       do k = 1,lm
        atmosphere%level_pressure(k) = pint(i,j,k+1)/r100
        atmosphere%pressure(k)       = pmid(i,j,k)/r100
        atmosphere%temperature(k)    = t(i,j,k)
!        atmosphere%absorber(k,1)     = q(i,j,k) ! use specific humidity
	atmosphere%absorber(k,1)     = q(i,j,k)*h1000/(h1-q(i,j,k)) ! use mixing ratio like GSI
        atmosphere%absorber(k,2)     = max(ozsmall,o3(i,j,k)*h1000) ! convert to g/kg
! fill in cloud mixing ratio later      
        if(i==ii.and.j==jj)print*,'sample atmosphere in CALRAD=',      &
     &	  i,j,k,atmosphere%level_pressure(k),atmosphere%pressure(k),   &
     &    atmosphere%temperature(k),atmosphere%absorber(k,1),          &
     &    atmosphere%absorber(k,2)
       end do

!     call crtm forward model
       error_status = crtm_forward(atmosphere,surface                  &
     &         ,geometryinfo,channelinfo,rtsolution)
       if (error_status /=0) then
        print*,'***ERROR*** during crtm_forward call ',                &
     &       error_status
       else 	 
        tb1(i,j)=rtsolution(1)%brightness_temperature
        tb2(i,j)=rtsolution(2)%brightness_temperature
        tb3(i,j)=rtsolution(3)%brightness_temperature	 
        tb4(i,j)=rtsolution(4)%brightness_temperature
	if(i==ii.and.j==jj)print*,'sample rtsolution in CALRAD=',      &
     &    rtsolution(1)%brightness_temperature,                        &
     &    rtsolution(2)%brightness_temperature,                        &
     &    rtsolution(3)%brightness_temperature,                        &
     &    rtsolution(4)%brightness_temperature
	if(i==ii.and.j==jj)print*,'sample TB in CALRAD=',             &
     &    tb1(i,j),tb2(i,j),tb3(i,j),tb4(i,j)
!        if(tb1(i,j) < 400. )
!     &        print*,'good tb1 ',i,j,tb1(i,j),gdlat(i,j),gdlon(i,j)
!        if(tb2(i,j) > 400.)print*,'bad tb2 ',i,j,tb2(i,j)
!        if(tb3(i,j) > 400.)print*,'bad tb3 ',i,j,tb3(i,j)
!        if(tb4(i,j) > 400.)print*,'bad tb4 ',i,j,tb4(i,j)
       end if  	 
       end do ! end loop for i
      end do ! end loop for j 
  
      error_status = crtm_destroy(channelinfo)
      if (error_status /= success)                                      &
     &   print*,'ERROR*** crtm_destroy error_status=',error_status
  
      if( iget(327) > 0 ) then
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb1(i,j)
         if(grid1(i,j)>400. .or. grid1(i,j)/=grid1(i,j))grid1(i,j)=400.
        enddo
       enddo
       id(1:25) = 0
       id(02) = 130
       call gribit(iget(327),lvls(1,iget(327)), grid1,im,jm)
      endif
  
      if( iget(328) > 0 ) then
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb2(i,j)
         if(grid1(i,j)>400.)grid1(i,j)=400.
        enddo
       enddo
       id(1:25) = 0
       id(02) = 130
       call gribit(iget(328),lvls(1,iget(328)),  grid1,im,jm)
      endif 
    
      if( iget(329) > 0 ) then
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb3(i,j)
         if(grid1(i,j)>400.)grid1(i,j)=400.
        enddo
       enddo
       id(1:25) = 0
       id(02) = 130
       call gribit(iget(329),lvls(1,iget(329)),  grid1,im,jm)
      endif  
    
      if( iget(330) > 0 ) then
       do j=jsta,jend
        do i=1,im
         grid1(i,j)=tb4(i,j)
         if(grid1(i,j)>400.)grid1(i,j)=400.
        enddo
       enddo
       id(1:25) = 0
       id(02) = 130
       call gribit(iget(330),lvls(1,iget(330)), grid1,im,jm)
      endif
    
! Deallocate arrays
      err1=0; err2=0; err3=0; err4=0
      err1 = crtm_destroy_atmosphere(atmosphere)
      err2 = crtm_destroy_surface(surface)
      err3 = crtm_destroy_rtsolution(rtsolution)
      if (err1/=success) write(6,*)'ERROR** destroy atmosphere.err=', &
     &      err1
      if (err2/=success) write(6,*)'ERROR** destroy surface.err=', &
     &      err2
      if (err3/=success) write(6,*)'ERROR** destroy rtsolution.err=', &
     &      err3
      deallocate (rtsolution)
!     

      endif ! for all iget logical
      return
      end
