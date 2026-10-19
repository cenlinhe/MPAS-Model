! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!=================================================================================================================
 module mpas_atmphys_sf_urbaninit

 use mpas_log
 use mpas_dmpar
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_io_units

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars
 use mpas_atmphys_landuse, only: isurban
!wrf physics
 use module_sf_urban

 implicit none
 private
 public:: init_sf_urban


!Initialization of the urban physics scheme.
!Cenlin He (send comments to cenlinhe@ucar.edu).
!2026-10-15.
!
! subroutines in mpas_atmphys_sf_urbaninit:
! -----------------------------------------
! init_sf_urban    : call urban_param_init and urban_var_init (defined in physics_wrf/module_sf_urban.F).
! urban_param_init : main subroutine to read in urban lookup table parameters.
! urban_var_init   : main initialization subroutine for urban state variables
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! add call to intialization subroutines defined in physics_wrf/module_sf_urban.F
! Cenlin He (cenlinhe@ucar.edu) / 2026-10-15.


 contains

!=================================================================================================================
 subroutine init_sf_urban(dminfo,mesh,configs,diag_physics,sfc_input,diag_physics_urban,diag_physics_slucm)
!=================================================================================================================

! input arguments:
 type(dm_info)       ,intent(in) :: dminfo
 type(mpas_pool_type),intent(in) :: mesh
 type(mpas_pool_type),intent(in) :: configs

! inout arguments:
 type(mpas_pool_type),intent(inout) :: diag_physics
 type(mpas_pool_type),intent(inout) :: sfc_input
 type(mpas_pool_type),intent(inout) :: diag_physics_urban
 type(mpas_pool_type),intent(inout) :: diag_physics_slucm

! local pointer
 logical,pointer :: slucm_distributed_drag,do_restart
 integer,pointer :: nSoilLevels,sf_urban_physics,use_wudapt_lcz
 integer,pointer :: num_roof_layers,num_wall_layers,num_road_layers
 integer,pointer :: nCellsSolve,nVertLevels
 integer,dimension(:),pointer :: ivgtyp,utype_urb2d

 real(kind=RKIND),dimension(:),pointer   :: dzr,dzb,dzg
 real(kind=RKIND),dimension(:),pointer   :: skintemp,tmn,frc_urb2d
 real(kind=RKIND),dimension(:),pointer   :: xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d
 real(kind=RKIND),dimension(:),pointer   :: tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d
 real(kind=RKIND),dimension(:),pointer   :: sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d
 real(kind=RKIND),dimension(:),pointer   :: hgt_urb2d,lb_urb2d,lp_urb2d
 real(kind=RKIND),dimension(:),pointer   :: mh_urb2d,stdh_urb2d,cmcr_urb2d,tgr_urb2d
 real(kind=RKIND),dimension(:),pointer   :: drelr_urb2d,drelb_urb2d,drelg_urb2d
 real(kind=RKIND),dimension(:),pointer   :: flxhumr_urb2d,flxhumb_urb2d,flxhumg_urb2d
 real(kind=RKIND),dimension(:),pointer   :: tvg_urb2d,xxxvg_urb2d,tt_urb2d
 real(kind=RKIND),dimension(:),pointer   :: cmcg_urb2d,flxhumvg_urb2d,flxhumt_urb2d
 real(kind=RKIND),dimension(:,:),pointer :: tslb,trl_urb3d,tbl_urb3d,tgl_urb3d
 real(kind=RKIND),dimension(:,:),pointer :: lf_urb2d,tgrl_urb3d,smr_urb3d
 real(kind=RKIND),dimension(:,:),pointer :: tvgl_urb3d,smg_urb3d


!-----------------------------------------------------------------------------------------------------------------

!-------- read in urban lookup table parameters --------
! generic MPAS variables
 call mpas_pool_get_config(configs,'config_urban_physics',sf_urban_physics)
 call mpas_pool_get_config(configs,'config_urban_use_lcz',use_wudapt_lcz  )
 call mpas_pool_get_config(configs,'config_slucm_distributed_drag',slucm_distributed_drag)
 call mpas_pool_get_dimension(mesh,'nSoilLevels' ,nSoilLevels    )  
 call mpas_pool_get_dimension(mesh,'nRoofLevels' ,num_roof_layers)
 call mpas_pool_get_dimension(mesh,'nWallLevels' ,num_wall_layers)
 call mpas_pool_get_dimension(mesh,'nRoadLevels' ,num_road_layers)
 call mpas_pool_get_array(sfc_input,'frc_urb2d'  ,frc_urb2d      )
 call mpas_pool_get_array(sfc_input,'utype_urb2d',utype_urb2d    )
! single-layer UCM diagnostic I/O variables
 call mpas_pool_get_array(diag_physics_slucm,'dzr', dzr )
 call mpas_pool_get_array(diag_physics_slucm,'dzb', dzb )
 call mpas_pool_get_array(diag_physics_slucm,'dzg', dzg )

! read urban lookup table parameters
 call urban_param_init(dzr,dzb,dzg,nSoilLevels,sf_urban_physics,use_wudapt_lcz, &
                       slucm_distributed_drag,num_roof_layers,num_wall_layers,num_road_layers)


!--------- call main urban initialization subroutine -----------
 if ( associated(frc_urb2d) .and. associated(utype_urb2d) ) then

     ! generic MPAS variables
     call mpas_pool_get_config(configs,'config_do_restart',do_restart )
     call mpas_pool_get_dimension(mesh,'nCellsSolve'      ,nCellsSolve)
     call mpas_pool_get_dimension(mesh,'nVertLevels'      ,nVertLevels)
     call mpas_pool_get_array(sfc_input,'skintemp'        ,skintemp   )
     call mpas_pool_get_array(sfc_input,'tslb'            ,tslb       )
     call mpas_pool_get_array(sfc_input,'tmn'             ,tmn        )
     call mpas_pool_get_array(sfc_input,'ivgtyp'          ,ivgtyp     )
     call mpas_pool_get_array(sfc_input,'lp_urb2d'        ,lp_urb2d   )
     call mpas_pool_get_array(sfc_input,'lb_urb2d'        ,lb_urb2d   )
     call mpas_pool_get_array(sfc_input,'hgt_urb2d'       ,hgt_urb2d  )
     call mpas_pool_get_array(sfc_input,'mh_urb2d'        ,mh_urb2d   )
     call mpas_pool_get_array(sfc_input,'stdh_urb2d'      ,stdh_urb2d )
     call mpas_pool_get_array(sfc_input,'lf_urb2d'        ,lf_urb2d   )
     ! generic urban variables
     call mpas_pool_get_array(diag_physics_urban,'sh_urb2d',sh_urb2d  )
     call mpas_pool_get_array(diag_physics_urban,'lh_urb2d',lh_urb2d  )
     call mpas_pool_get_array(diag_physics_urban,'g_urb2d' ,g_urb2d   )
     call mpas_pool_get_array(diag_physics_urban,'rn_urb2d',rn_urb2d  )
     call mpas_pool_get_array(diag_physics_urban,'ts_urb2d',ts_urb2d  )
     ! Single-layer UCM diagnostic I/O variables
     call mpas_pool_get_array(diag_physics_slucm,'xxxr_urb2d'    ,xxxr_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'xxxb_urb2d'    ,xxxb_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'xxxg_urb2d'    ,xxxg_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'xxxc_urb2d'    ,xxxc_urb2d    )
     call mpas_pool_get_array(diag_physics_slucm,'tr_urb2d'      ,tr_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'tb_urb2d'      ,tb_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'tg_urb2d'      ,tg_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'tc_urb2d'      ,tc_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'qc_urb2d'      ,qc_urb2d      )
     call mpas_pool_get_array(diag_physics_slucm,'trl_urb3d'     ,trl_urb3d     )
     call mpas_pool_get_array(diag_physics_slucm,'tbl_urb3d'     ,tbl_urb3d     )
     call mpas_pool_get_array(diag_physics_slucm,'tgl_urb3d'     ,tgl_urb3d     )
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

     ! main initialization call
     call urban_var_init(ISURBAN = isurban, TSURFACE0_URB = skintemp, TLAYER0_URB = tslb,    &
                         TDEEP0_URB = tmn, IVGTYP = ivgtyp, its = 1, ite = nCellsSolve,      &
                         kts = 1, kte = nVertLevels, num_soil_layers = nSoilLevels,          &
                         num_roof_layers = num_roof_layers,num_wall_layers = num_wall_layers,&
                         num_road_layers = num_road_layers,                                  &
                         LCZ_1 = mpas_noahmp%LCZ_1_TABLE, LCZ_2 = mpas_noahmp%LCZ_2_TABLE,   &
                         LCZ_3 = mpas_noahmp%LCZ_3_TABLE, LCZ_4 = mpas_noahmp%LCZ_4_TABLE,   &
                         LCZ_5 = mpas_noahmp%LCZ_5_TABLE, LCZ_6 = mpas_noahmp%LCZ_6_TABLE,   &
                         LCZ_7 = mpas_noahmp%LCZ_7_TABLE, LCZ_8 = mpas_noahmp%LCZ_8_TABLE,   &
                         LCZ_9 = mpas_noahmp%LCZ_9_TABLE, LCZ_10= mpas_noahmp%LCZ_10_TABLE,  &
                         LCZ_11= mpas_noahmp%LCZ_11_TABLE, restart = do_restart,             &
                         sf_urban_physics = sf_urban_physics, XXXR_URB2D = xxxr_urb2d,       &
                         XXXB_URB2D = xxxb_urb2d, XXXG_URB2D = xxxg_urb2d,                   &
                         XXXC_URB2D = xxxc_urb2d, TR_URB2D = tr_urb2d, TB_URB2D = tb_urb2d,  &
                         TG_URB2D = tg_urb2d, TC_URB2D = tc_urb2d, QC_URB2D = qc_urb2d,      &
                         TRL_URB3D = trl_urb3d, TBL_URB3D = tbl_urb3d,TGL_URB3D = tgl_urb3d, &
                         SH_URB2D = sh_urb2d, LH_URB2D = lh_urb2d, G_URB2D = g_urb2d,        &
                         RN_URB2D = rn_urb2d, TS_URB2D = ts_urb2d, LP_URB2D = lp_urb2d,      &
                         LB_URB2D = lb_urb2d, HGT_URB2D = hgt_urb2d, MH_URB2D = mh_urb2d,    &
                         STDH_URB2D = stdh_urb2d,LF_URB2D = lf_urb2d,CMCR_URB2D = cmcr_urb2d,&
                         TGR_URB2D = tgr_urb2d,TGRL_URB3D = tgrl_urb3d,SMR_URB3D = smr_urb3d,&
                         DRELR_URB2D = drelr_urb2d, DRELB_URB2D = drelb_urb2d,               &
                         DRELG_URB2D = drelg_urb2d, FLXHUMR_URB2D = flxhumr_urb2d,           &
                         FLXHUMB_URB2D = flxhumb_urb2d, FLXHUMG_URB2D = flxhumg_urb2d,       &
                         TVG_URB2D = tvg_urb2d, XXXVG_URB2D = xxxvg_urb2d,                   &
                         TVGL_URB3D = tvgl_urb3d, SMG_URB3D = smg_urb3d,TT_URB2D = tt_urb2d, &
                         CMCG_URB2D = cmcg_urb2d, FLXHUMVG_URB2D = flxhumvg_urb2d,           &
                         FLXHUMT_URB2D = flxhumt_urb2d, FRC_URB2D = frc_urb2d,               &
                         UTYPE_URB2D = utype_urb2d, USE_WUDAPT_LCZ = use_wudapt_lcz          &
                         !--- multi-layer urban vars below (make optional; no BEP,BEM coupling now)
                         !num_urban_ndm, urban_map_zrd, urban_map_zwd, urban_map_gd,         &
                         !urban_map_zd, urban_map_zdf, urban_map_bd, urban_map_wd,           &
                         !urban_map_gbd, urban_map_fbd, urban_map_zgrd, num_urban_hi,        &
                         !TRB_URB4D, TW1_URB4D, TW2_URB4D, TGB_URB4D, TLEV_URB3D,QLEV_URB3D, &
                         !TW1LEV_URB3D, TW2LEV_URB3D, TGLEV_URB3D, TFLEV_URB3D, SF_AC_URB3D, &
                         !LF_AC_URB3D, CM_AC_URB3D, SFVENT_URB3D, LFVENT_URB3D,SFWIN1_URB3D, &
                         !SFWIN2_URB3D, SFW1_URB3D, SFW2_URB3D, SFR_URB3D, SFG_URB3D,        &
                         !EP_PV_URB3D, T_PV_URB3D, TRV_URB4D, QR_URB4D, QGR_URB3D,TGR_URB3D, &
                         !DRAIN_URB4D, DRAINGR_URB3D, SFRV_URB3D, LFRV_URB3D, DGR_URB3D,     &
                         !DG_URB3D, LFR_URB3D, LFG_URB3D, SMOIS_URB, HI_URB2D, A_U_BEP,      &
                         !A_V_BEP, A_T_BEP, A_Q_BEP, A_E_BEP, B_U_BEP, B_V_BEP, B_T_BEP,     &
                         !B_Q_BEP, B_E_BEP, DLG_BEP, DL_U_BEP, SF_BEP, VL_BEP,               &
                        )

 else

     call physics_error_fatal('arguments not present for calling urban model')

 endif ! associated(frc_urb2d)&(utype_urb2d)

 end subroutine init_sf_urban

!=================================================================================================================
 end module mpas_atmphys_sf_urbaninit
!=================================================================================================================
