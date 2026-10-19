! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================

!=================================================================================================================
 module mpas_atmphys_driver_sf_urban
!=================================================================================================================

 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer,           only : mpas_timer_start,mpas_timer_stop
 use mpas_timekeeping,     only : mpas_get_timeInterval,mpas_get_clock_timestep
 use mpas_atmphys_manager, only : gmt,curr_julday,year,month,day 

 use mpas_atmphys_constants
 use mpas_atmphys_vars     ! mpas_noahmp is included here

 use NoahmpUrbanDriverMainMod, only : noahmp_urban

 implicit none
 private
 public:: driver_sf_urban
 
!
! MPAS driver for parameterization of urban canopy model (UCM).
! Cenlin He (cenlinhe@ucar.edu).
! 2026-10-05
!
! subroutines in mpas_atmphys_driver_sf_urban:
! ---------------------------------------
! driver_urb_slucm : main driver (called from subroutine physics_driver).
!
! WRF physics called from driver_urb_slucm:
! ------------------------ ----------
!    * noahmp_urban subroutine in NoahmpUrbanDriverMainMod : urban canopy scheme
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * slucm source code and driver changed to MPAS 1-D vector (i) from WRF 2-D grids (i,j)
!   bep and bem code inside noahmp_urban driver is commented out to deactivate
! Cenlin He (cenlinhe@ucar.edu) / 2026-10-05.

! DOCUMENTATION:
! ./physics_noahmp/drivers/mpas/NoahmpUrbanDriverMainMod.F90: main driver for the urban parameterization.
! In the argument list, "OPTIONAL" is added to the declaration of a few arrays to avoid compiling
! with the bep and bem options which need future coupling.
! Cenlin He (2026-10-05).
!

 contains

!=================================================================================================================
 subroutine driver_sf_urban(configs,mesh,diag,state,itimestep,time_lev,xtime_s, &
                            sfc_input,diag_physics,diag_physics_slucm,output_slucm, &
                            diag_physics_urban,output_urban,output_noahmp,its,ite)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in) :: configs
 type(mpas_pool_type),intent(in) :: mesh
 type(mpas_pool_type),intent(in) :: diag
 type(mpas_pool_type),intent(in) :: state

 integer,         intent(in) :: itimestep,its,ite,time_lev
 real(kind=RKIND),intent(in) :: xtime_s
 
!--- inout arguments:
 type(mpas_pool_type),intent(inout) :: diag_physics
 type(mpas_pool_type),intent(inout) :: diag_physics_slucm
 type(mpas_pool_type),intent(inout) :: output_slucm
 type(mpas_pool_type),intent(inout) :: sfc_input
 type(mpas_pool_type),intent(inout) :: diag_physics_urban
 type(mpas_pool_type),intent(inout) :: output_urban
 type(mpas_pool_type),intent(inout) :: output_noahmp

! local variables:
 integer:: i,k,ierr
 character(len=StrKIND),pointer :: microp_scheme,convection_scheme
 real(kind=RKIND) :: fzm,fzp,mult,totm,totp
 real(kind=RKIND) :: w1,w2,z0,z1,z2
 real(kind=RKIND) :: radt,xtime_m,xt24
 integer, pointer :: index_qv

 real(kind=RKIND),dimension(:),    pointer :: coszr,latCell,lonCell,swddir,swddif,glw,gsw
 real(kind=RKIND),dimension(:,:),  pointer :: zgrid
 real(kind=RKIND),dimension(:,:),  pointer :: qv,theta_m,u,v
 real(kind=RKIND),dimension(:,:),  pointer :: exner,pressure_b,pressure_p
 real(kind=RKIND),dimension(:,:),  pointer :: dz8w,t_phy,p8w
 real(kind=RKIND),dimension(:,:,:),pointer :: scalars

 integer,pointer :: sf_urban_physics,nSoilLevels,config_distributed_ahe
 integer,pointer :: num_roof_layers,num_wall_layers,num_road_layers
 integer, dimension(:), pointer :: ivgtyp,utype_urb2d

 real(kind=RKIND),pointer :: dt

 real(kind=RKIND),dimension(:),pointer :: swdown,rainbl,xlat,xlong,hrang,vegfra
 real(kind=RKIND),dimension(:),pointer :: raincv,rainncv,znt,skintemp,qsfc
 real(kind=RKIND),dimension(:),pointer :: hfx,qfx,lh,grdflx,sfc_emiss,sfc_albedo
 real(kind=RKIND),dimension(:),pointer :: cmr_sfcdif,chr_sfcdif,cmc_sfcdif
 real(kind=RKIND),dimension(:),pointer :: chc_sfcdif,cmgr_sfcdif,chgr_sfcdif
 real(kind=RKIND),dimension(:),pointer :: tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d
 real(kind=RKIND),dimension(:),pointer :: qc_urb2d,uc_urb2d,q2_urb2d,ust_urb2d
 real(kind=RKIND),dimension(:),pointer :: xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d
 real(kind=RKIND),dimension(:),pointer :: sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d
 real(kind=RKIND),dimension(:),pointer :: psim_urb2d,psih_urb2d,u10_urb2d,v10_urb2d 
 real(kind=RKIND),dimension(:),pointer :: gz1oz0_urb2d,akms_urb2d,th2_urb2d
 real(kind=RKIND),dimension(:),pointer :: dzr,dzb,dzg,flxhumr_urb2d,flxhumb_urb2d
 real(kind=RKIND),dimension(:),pointer :: cmcr_urb2d,tgr_urb2d,flxhumg_urb2d
 real(kind=RKIND),dimension(:),pointer :: drelr_urb2d,drelb_urb2d,drelg_urb2d
 real(kind=RKIND),dimension(:),pointer :: tvg_urb2d,xxxvg_urb2d,tt_urb2d
 real(kind=RKIND),dimension(:),pointer :: cmcg_urb2d,frc_urb2d,mh_urb2d
 real(kind=RKIND),dimension(:),pointer :: flxhumvg_urb2d,flxhumt_urb2d
 real(kind=RKIND),dimension(:),pointer :: hgt_urb2d,lb_urb2d,lp_urb2d
 real(kind=RKIND),dimension(:),pointer :: stdh_urb2d,lf_urb2d_s,z0_urb2d
 real(kind=RKIND),dimension(:),pointer :: t2mxy,q2mxy,psfc,th2m,u10,v10,psim,psih,gz1oz0
 real(kind=RKIND),dimension(:),pointer :: chs,chs2,cqs2,ust
                    
 real(kind=RKIND),dimension(:,:),pointer :: trl_urb3d,tgrl_urb3d,smr_urb3d
 real(kind=RKIND),dimension(:,:),pointer :: tbl_urb3d
 real(kind=RKIND),dimension(:,:),pointer :: tgl_urb3d,tvgl_urb3d,smg_urb3d
 real(kind=RKIND),dimension(:,:),pointer :: lf_urb2d,ahe_urb2d

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_sf_urban:')

 call mpas_pool_get_config(configs,'config_urban_physics',sf_urban_physics)

 if ( sf_urban_physics > 0 ) then

     call mpas_timer_start('sf_urban')

     !--- compute atmos forcing (no direct corresponding mpas vars)
     !--- consistent with Noah-MP forcing processing in mpas_atmphys_driver_lsm_noahmp.F
     call mpas_pool_get_config(configs, 'config_convection_scheme',convection_scheme)
     call mpas_pool_get_config(configs, 'config_microp_scheme'    ,microp_scheme    )
     call mpas_pool_get_array(mesh,     'zgrid'                   ,zgrid            )
     call mpas_pool_get_array(mesh,     'latCell'                 ,latCell          )
     call mpas_pool_get_array(mesh,     'lonCell'                 ,lonCell          )
     call mpas_pool_get_array(diag,     'exner'                   ,exner            )
     call mpas_pool_get_array(diag,     'pressure_base'           ,pressure_b       )
     call mpas_pool_get_array(diag,     'pressure_p'              ,pressure_p       )
     call mpas_pool_get_array(diag,     'uReconstructZonal'       ,u                )
     call mpas_pool_get_array(diag,     'uReconstructMeridional'  ,v                )
     call mpas_pool_get_array(state,    'theta_m'                 ,theta_m,time_lev )
     call mpas_pool_get_array(state,    'scalars'                 ,scalars,time_lev )
     call mpas_pool_get_dimension(state,'index_qv'                ,index_qv         )
     call mpas_pool_get_array(sfc_input,'xlat'                    ,xlat             )
     call mpas_pool_get_array(sfc_input,'xlong'                   ,xlong            )
     call mpas_pool_get_array(sfc_input,'hrang'                   ,hrang            )
     call mpas_pool_get_array(diag_physics,'gsw'                  ,gsw              )
     call mpas_pool_get_array(diag_physics,'sfc_albedo'           ,sfc_albedo       )
     call mpas_pool_get_array(diag_physics,'raincv'               ,raincv           )
     call mpas_pool_get_array(diag_physics,'rainncv'              ,rainncv          )
     call mpas_pool_get_array(output_urban,'dz8w'                 ,dz8w             )
     call mpas_pool_get_array(output_urban,'p8w'                  ,p8w              )
     call mpas_pool_get_array(output_urban,'t_phy'                ,t_phy            )
     call mpas_pool_get_array(output_urban,'swdown'               ,swdown           )
     call mpas_pool_get_array(output_urban,'rainbl'               ,rainbl           )

     ! initialize layer thickness, temperature, radiation, humidity
     qv => scalars(index_qv,:,:)

     do i = its,ite
        do k = kts,kte
           dz8w(k,i)  = zgrid(k+1,i)-zgrid(k,i)
           t_phy(k,i) = (theta_m(k,i) / (1.0 + R_v/R_d*qv(k,i))) * exner(k,i)
        enddo
        swdown(i) = gsw(i) / (1.0 - sfc_albedo(i))
     enddo

     ! initialization of pressure at interface between layers:
     do i = its,ite
        k = kts
        z0 = zgrid(k,i)
        z1 = 0.5*(zgrid(k,i)+zgrid(k+1,i))
        z2 = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
        w1 = (z0-z2)/(z1-z2)
        w2 = 1.-w1
        totm = pressure_p(k,i)+pressure_b(k,i)
        totp = pressure_p(k+1,i)+pressure_b(k+1,i)
        p8w(k,i) = w1*totm + w2*totp
        do k = kts+1,kte
           totm = pressure_p(k-1,i)+pressure_b(k-1,i)
           totp = pressure_p(k,i)+pressure_b(k,i)
           mult = 1.0 / (zgrid(k+1,i)-zgrid(k-1,i))
           fzm  = mult*(zgrid(k,i)-zgrid(k-1,i))
           fzp  = mult*(zgrid(k+1,i)-zgrid(k,i))
           p8w(k,i) = fzm*totp + fzp*totm
        enddo
     enddo

     ! compute total precipitation
     do i = its,ite
        rainbl(i) = 0.0
     enddo
     if (microp_scheme .ne. 'off') then
        do i = its,ite
           rainbl(i) = rainbl(i) + rainncv(i)
        enddo
     endif
     if (convection_scheme .ne. 'off') then
        do i = its,ite
           rainbl(i) = rainbl(i) + raincv(i)
        enddo
     endif

     ! compute radiation variables follow module_ra_rrtmg_sw.F
     xtime_m = xtime_s / 60.    ! second -> minute
     radt    = dt_radtsw / 60.  ! second -> minute
     xt24 = mod(xtime_m + radt*0.5, 1440.)
     do i = its,ite
        xlat(i)  = latCell(i) / degrad
        xlong(i) = lonCell(i) / degrad
        hrang(i) = 15. *((gmt + xt24/60. + xlong(i)/15.) - 12.) * degrad
     enddo

     !--- forcing compute done

     ! pass in other direct mpas variables
     ! generic mpas variables
     call mpas_pool_get_config(configs,'config_distributed_ahe',config_distributed_ahe)
     !call mpas_get_timeInterval(mpas_get_clock_timestep(clock, ierr), dt=dt)
     call mpas_pool_get_config(configs,'config_dt',dt)
     call mpas_pool_get_dimension(mesh ,'nSoilLevels',nSoilLevels    )
     call mpas_pool_get_dimension(mesh ,'nRoofLevels',num_roof_layers)
     call mpas_pool_get_dimension(mesh ,'nWallLevels',num_wall_layers)
     call mpas_pool_get_dimension(mesh ,'nRoadLevels',num_road_layers)
     call mpas_pool_get_array(sfc_input,'ivgtyp'     ,ivgtyp         )
     call mpas_pool_get_array(sfc_input,'skintemp'   ,skintemp       )
     call mpas_pool_get_array(sfc_input,'frc_urb2d'  ,frc_urb2d      )
     call mpas_pool_get_array(sfc_input,'utype_urb2d',utype_urb2d    )
     call mpas_pool_get_array(sfc_input,'lp_urb2d'   ,lp_urb2d       )
     call mpas_pool_get_array(sfc_input,'lb_urb2d'   ,lb_urb2d       )
     call mpas_pool_get_array(sfc_input,'hgt_urb2d'  ,hgt_urb2d      )
     call mpas_pool_get_array(sfc_input,'mh_urb2d'   ,mh_urb2d       )
     call mpas_pool_get_array(sfc_input,'stdh_urb2d' ,stdh_urb2d     )
     call mpas_pool_get_array(sfc_input,'lf_urb2d'   ,lf_urb2d       )
     call mpas_pool_get_array(sfc_input,'lf_urb2d_s' ,lf_urb2d_s     )
     call mpas_pool_get_array(sfc_input,'z0_urb2d'   ,z0_urb2d       )
     call mpas_pool_get_array(sfc_input,'vegfra'     ,vegfra         )
     call mpas_pool_get_array(sfc_input,'ahe_urb2d'  ,ahe_urb2d      )
     call mpas_pool_get_array(diag     ,'surface_pressure',psfc      )
     ! mpas physics variable
     call mpas_pool_get_array(diag_physics,'coszr'    ,coszr    )
     call mpas_pool_get_array(diag_physics,'swddir'   ,swddir   )
     call mpas_pool_get_array(diag_physics,'swddif'   ,swddif   )
     call mpas_pool_get_array(diag_physics,'glw'      ,glw      )
     call mpas_pool_get_array(diag_physics,'znt'      ,znt      )
     call mpas_pool_get_array(diag_physics,'hfx'      ,hfx      )
     call mpas_pool_get_array(diag_physics,'qfx'      ,qfx      )
     call mpas_pool_get_array(diag_physics,'lh '      ,lh       )
     call mpas_pool_get_array(diag_physics,'grdflx'   ,grdflx   )
     call mpas_pool_get_array(diag_physics,'sfc_emiss',sfc_emiss)
     call mpas_pool_get_array(diag_physics,'qsfc'     ,qsfc     )
     call mpas_pool_get_array(diag_physics,'th2m'     ,th2m     )
     call mpas_pool_get_array(diag_physics,'psim'     ,psim     )
     call mpas_pool_get_array(diag_physics,'psih'     ,psih     )
     call mpas_pool_get_array(diag_physics,'u10'      ,u10      )
     call mpas_pool_get_array(diag_physics,'v10'      ,v10      )
     call mpas_pool_get_array(diag_physics,'gz1oz0'   ,gz1oz0   )
     call mpas_pool_get_array(diag_physics,'chs'      ,chs      )
     call mpas_pool_get_array(diag_physics,'chs2'     ,chs2     )
     call mpas_pool_get_array(diag_physics,'cqs2'     ,cqs2     )
     call mpas_pool_get_array(diag_physics,'ust'      ,ust      )
     ! noahmp diagnostic variable
     call mpas_pool_get_array(output_noahmp,'t2mxy'   ,t2mxy    )
     call mpas_pool_get_array(output_noahmp,'q2mxy'   ,q2mxy    )
     ! generic urban variables
     call mpas_pool_get_array(diag_physics_urban,'sh_urb2d',sh_urb2d)
     call mpas_pool_get_array(diag_physics_urban,'lh_urb2d',lh_urb2d)
     call mpas_pool_get_array(diag_physics_urban,'g_urb2d' ,g_urb2d )
     call mpas_pool_get_array(diag_physics_urban,'rn_urb2d',rn_urb2d)
     call mpas_pool_get_array(diag_physics_urban,'ts_urb2d',ts_urb2d)
     ! Single-layer UCM diagnostic I/O variables
     call mpas_pool_get_array(diag_physics_slucm,'cmr_sfcdif'    ,cmr_sfcdif    )
     call mpas_pool_get_array(diag_physics_slucm,'chr_sfcdif'    ,chr_sfcdif    )
     call mpas_pool_get_array(diag_physics_slucm,'cmc_sfcdif'    ,cmc_sfcdif    )
     call mpas_pool_get_array(diag_physics_slucm,'chc_sfcdif'    ,chc_sfcdif    )
     call mpas_pool_get_array(diag_physics_slucm,'cmgr_sfcdif'   ,cmgr_sfcdif   )
     call mpas_pool_get_array(diag_physics_slucm,'chgr_sfcdif'   ,chgr_sfcdif   )
     call mpas_pool_get_array(diag_physics_slucm,'tr_urb2d'      ,tr_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'tb_urb2d'      ,tb_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'tg_urb2d'      ,tg_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'tc_urb2d'      ,tc_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'qc_urb2d'      ,qc_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'uc_urb2d'      ,uc_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'xxxr_urb2d'    ,xxxr_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'xxxb_urb2d'    ,xxxb_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'xxxg_urb2d'    ,xxxg_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'xxxc_urb2d'    ,xxxc_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'trl_urb3d'     ,trl_urb3d     )
     call mpas_pool_get_array(diag_physics_slucm,'tbl_urb3d'     ,tbl_urb3d     )
     call mpas_pool_get_array(diag_physics_slucm,'tgl_urb3d'     ,tgl_urb3d     )
     call mpas_pool_get_array(diag_physics_slucm,'dzr'           ,dzr           )
     call mpas_pool_get_array(diag_physics_slucm,'dzb'           ,dzb           )
     call mpas_pool_get_array(diag_physics_slucm,'dzg'           ,dzg           )
     call mpas_pool_get_array(diag_physics_slucm,'cmcr_urb2d'    ,cmcr_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'tgr_urb2d'     ,tgr_urb2d     )
     call mpas_pool_get_array(diag_physics_slucm,'tgrl_urb3d'    ,tgrl_urb3d    )
     call mpas_pool_get_array(diag_physics_slucm,'smr_urb3d'     ,smr_urb3d     )
     call mpas_pool_get_array(diag_physics_slucm,'drelr_urb2d'   ,drelr_urb2d   )
     call mpas_pool_get_array(diag_physics_slucm,'drelb_urb2d'   ,drelb_urb2d   )
     call mpas_pool_get_array(diag_physics_slucm,'drelg_urb2d'   ,drelg_urb2d   )
     call mpas_pool_get_array(diag_physics_slucm,'flxhumr_urb2d' ,flxhumr_urb2d )
     call mpas_pool_get_array(diag_physics_slucm,'flxhumb_urb2d' ,flxhumb_urb2d )
     call mpas_pool_get_array(diag_physics_slucm,'flxhumg_urb2d' ,flxhumg_urb2d )
     call mpas_pool_get_array(diag_physics_slucm,'tvg_urb2d'     ,tvg_urb2d     )
     call mpas_pool_get_array(diag_physics_slucm,'xxxvg_urb2d'   ,xxxvg_urb2d   )
     call mpas_pool_get_array(diag_physics_slucm,'tvgl_urb3d'    ,tvgl_urb3d    )
     call mpas_pool_get_array(diag_physics_slucm,'smg_urb3d'     ,smg_urb3d     )
     call mpas_pool_get_array(diag_physics_slucm,'tt_urb2d'      ,tt_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'cmcg_urb2d'    ,cmcg_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'flxhumvg_urb2d',flxhumvg_urb2d)
     call mpas_pool_get_array(diag_physics_slucm,'flxhumt_urb2d' ,flxhumt_urb2d )
     !  Single-layer UCM output-only variable
     call mpas_pool_get_array(output_slucm,'psim_urb2d'  ,psim_urb2d  )
     call mpas_pool_get_array(output_slucm,'psih_urb2d'  ,psih_urb2d  )
     call mpas_pool_get_array(output_slucm,'u10_urb2d'   ,u10_urb2d   )
     call mpas_pool_get_array(output_slucm,'v10_urb2d'   ,v10_urb2d   )
     call mpas_pool_get_array(output_slucm,'gz1oz0_urb2d',gz1oz0_urb2d)
     call mpas_pool_get_array(output_slucm,'akms_urb2d'  ,akms_urb2d  )
     call mpas_pool_get_array(output_slucm,'th2_urb2d'   ,th2_urb2d   )
     call mpas_pool_get_array(output_slucm,'q2_urb2d'    ,q2_urb2d    )
     call mpas_pool_get_array(output_slucm,'ust_urb2d'   ,ust_urb2d   )

     ! run urban model
     call noahmp_urban(NoahmpIO = mpas_noahmp, sf_urban_physics = sf_urban_physics,            &
                       NSOIL = nSoilLevels, IVGTYP = ivgtyp, ITIMESTEP = itimestep,            &
                       DT = dt, COSZ_URB2D = coszr, XLAT_URB2D = xlat,                         &
                       T3D = t_phy, QV3D = qv, U_PHY = u, V_PHY = v, SWDOWN = swdown,          &
                       SWDDIR = swddir, SWDDIF = swddif, GLW = glw, P8W3D = p8w,               &
                       RAINBL = rainbl, DZ8W = dz8w, ZNT = znt, TSK = skintemp,                &
                       HFX = hfx, QFX = qfx, LH = lh, GRDFLX = grdflx, ALBEDO = sfc_albedo,    &
                       EMISS = sfc_emiss, QSFC = qsfc, its = its, ite = ite, kts = kts,        &
                       kte = kte, cmr_sfcdif = cmr_sfcdif, chr_sfcdif = chr_sfcdif,            &
                       cmc_sfcdif = cmc_sfcdif, chc_sfcdif = chc_sfcdif,                       &
                       cmgr_sfcdif = cmgr_sfcdif, chgr_sfcdif = chgr_sfcdif,                   &
                       tr_urb2d = tr_urb2d, tb_urb2d = tb_urb2d,                               &
                       tg_urb2d = tg_urb2d, tc_urb2d = tc_urb2d, qc_urb2d = qc_urb2d,          &
                       uc_urb2d = uc_urb2d, xxxr_urb2d = xxxr_urb2d, xxxb_urb2d = xxxb_urb2d,  &
                       xxxg_urb2d = xxxg_urb2d, xxxc_urb2d = xxxc_urb2d,trl_urb3d = trl_urb3d, &
                       tbl_urb3d = tbl_urb3d, tgl_urb3d = tgl_urb3d, sh_urb2d = sh_urb2d,      &
                       lh_urb2d = lh_urb2d, g_urb2d = g_urb2d, rn_urb2d = rn_urb2d,            &
                       ts_urb2d = ts_urb2d, psim_urb2d = psim_urb2d, psih_urb2d = psih_urb2d,  &
                       u10_urb2d = u10_urb2d, v10_urb2d = v10_urb2d,                           &
                       GZ1OZ0_urb2d = gz1oz0_urb2d, AKMS_URB2D = akms_urb2d,                   &
                       th2_urb2d = th2_urb2d, q2_urb2d = q2_urb2d, ust_urb2d = ust_urb2d,      &
                       declin_urb = declin, omg_urb2d = hrang,                                 &
                       num_roof_layers = num_roof_layers, num_wall_layers = num_wall_layers,   &
                       num_road_layers = num_road_layers, dzr = dzr, dzb = dzb, dzg = dzg,     &
                       cmcr_urb2d = cmcr_urb2d,tgr_urb2d = tgr_urb2d,tgrl_urb3d = tgrl_urb3d,  &
                       smr_urb3d = smr_urb3d, drelr_urb2d = drelr_urb2d,                       &
                       drelb_urb2d = drelb_urb2d, drelg_urb2d = drelg_urb2d,                   &
                       flxhumr_urb2d = flxhumr_urb2d, flxhumb_urb2d = flxhumb_urb2d,           &
                       flxhumg_urb2d = flxhumg_urb2d, tvg_urb2d = tvg_urb2d,                   &
                       xxxvg_urb2d = xxxvg_urb2d,tvgl_urb3d = tvgl_urb3d,smg_urb3d = smg_urb3d,&
                       tt_urb2d = tt_urb2d, cmcg_urb2d = cmcg_urb2d,                           &
                       flxhumvg_urb2d = flxhumvg_urb2d, flxhumt_urb2d = flxhumt_urb2d,         &
                       julian = int(curr_julday), julyr = year, frc_urb2d = frc_urb2d,         &
                       utype_urb2d = utype_urb2d, chs = chs, chs2 = chs2, cqs2 = cqs2,         &
                       lp_urb2d = lp_urb2d, lb_urb2d = lb_urb2d, hgt_urb2d = hgt_urb2d,        &
                       mh_urb2d = mh_urb2d, stdh_urb2d = stdh_urb2d, lf_urb2d = lf_urb2d,      &
                       lf_urb2d_s = lf_urb2d_s, z0_urb2d = z0_urb2d, vegfra = vegfra,          &
                       T2 = t2mxy, Q2 = q2mxy, PSFC = psfc, TH2 = th2m, U10 = u10, V10 = v10,  &
                       PSIM = psim, PSIH = psih, GZ1OZ0 = gz1oz0,                              &
                       distributed_ahe_opt = config_distributed_ahe, gmt = gmt,                &
                       xtime =xtime_m, ahe_urb2d = ahe_urb2d, ust = ust,                       &
                       !--- multi-layer urban vars below (make optional; no BEP,BEM coupling now)
                       !num_urban_ndm, urban_map_zrd, urban_map_zwd, urban_map_gd,             &
                       !urban_map_zd, urban_map_zdf,  urban_map_bd, urban_map_wd,              &
                       !urban_map_gbd, urban_map_fbd, urban_map_zgrd, num_urban_hi,            &
                       !trb_urb4d, tw1_urb4d, tw2_urb4d, tgb_urb4d,                            &
                       !tlev_urb3d, qlev_urb3d, tw1lev_urb3d,   tw2lev_urb3d,                  &
                       !tglev_urb3d, tflev_urb3d, sf_ac_urb3d, lf_ac_urb3d, cm_ac_urb3d,       &
                       !sfvent_urb3d, lfvent_urb3d, sfwin1_urb3d, sfwin2_urb3d,                &
                       !sfw1_urb3d, sfw2_urb3d, sfr_urb3d, sfg_urb3d, ep_pv_urb3d, t_pv_urb3d, & 
                       !trv_urb4d, qr_urb4d, qgr_urb3d, tgr_urb3d, drain_urb4d, draingr_urb3d, &
                       !sfrv_urb3d, lfrv_urb3d, dgr_urb3d, dg_urb3d, lfr_urb3d, lfg_urb3d,     &
                       !hi_urb2d, th_phy, rho, p_phy, julday, xlong, xlat,                     &
                       !a_u_bep, a_v_bep, a_t_bep, a_q_bep, a_e_bep, b_u_bep, b_v_bep,         &
                       !b_t_bep, b_q_bep, b_e_bep, dlg_bep, dl_u_bep, sf_bep, vl_bep           &
                      )

    call mpas_timer_stop('sf_urban')

 endif ! iopt_urban_physics > 0

!call mpas_log_write('--- end subroutine driver_sf_urban.')

 end subroutine driver_sf_urban

!=================================================================================================================
 end module mpas_atmphys_driver_sf_urban
!=================================================================================================================
