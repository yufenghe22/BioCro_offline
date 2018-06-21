!CWRF:MEDIATION_LAYER:PHYSICS
!

#define LaxSFC

!+xum add options
#undef DynFVC
!=xum

MODULE module_surface_driver

!----------------------------------------------------------------------
   USE module_wrf_error
   USE module_state_description , ONLY: SFCLAYSCHEME                &
                                       ,MYNNSFCSCHEME               &
                                       ,MYJSFCSCHEME                &
                                       ,QNSESFCSCHEME               &
                                       ,PXSFCSCHEME                 &
                                       ,GFSSFCSCHEME                &
                                       ,GFDLSFCSCHEME               &
                                       ,SLABSCHEME                  &
                                       ,GFDLSLAB                    &
                                       ,SOMSCHEME                   &
                                       ,CSSPSCHEME                  &
                                       ,UCMSCHEME                   &
                                       ,BEPSCHEME                   &
                                       ,NMMLSMSCHEME                &
                                       ,LSMSCHEME                   &
                                       ,RUCLSMSCHEME                &
                                       ,PXLSMSCHEME                 &
                                       ,MYJPBLSCHEME,MYJURBSCHEME   &
                                       ,MYNNPBLSCHEME               &
                                       ,QNSEPBLSCHEME

   USE module_sf_sbc            , ONLY: EXT_SBC
   USE module_sf_som            , ONLY: SOM, SSTSKIN , sominit
   USE module_sf_cssp           , ONLY: CSSP         , csspinit
#ifdef DoVAST
   USE module_sf_cssp           , ONLY: horizontalflow_vast
#endif
   USE module_sf_noahlsm        , ONLY: LSM          , lsminit
   USE module_sf_ruclsm         , ONLY: LSMRUC       , lsmrucinit
   USE module_sf_pxlsm          , ONLY: PXLSM
   USE module_sf_lsm_nmm        , ONLY: NMMLSM       , nmmlsminit
   USE module_sf_gfdl           , ONLY: SF_GFDL
   USE module_sf_mynnsfc        , ONLY: MYNNSFC      , mynnsfcinit
   USE module_sf_sfclay         , ONLY: SFCLAY       , sfclayinit
   USE module_sf_pxsfclay       , ONLY: PXSFCLAY     , pxsfclayinit
   USE module_sf_myjsfc         , ONLY: MYJSFC       , myjsfcinit
   USE module_sf_qnsesfc        , ONLY: QNSESFC      , qnsesfcinit
   USE module_sf_gfs            , ONLY: SF_GFS       , sf_gfsinit
   USE module_sf_slab           , ONLY: SLAB         , slabinit
   USE module_sf_sfcdiags       , ONLY: SFCDIAGS, SFCSTRES, SFCVISCOUS, SFCETP
#ifdef CalUCM
   USE module_sf_urban          , ONLY: UCM          , ucminit
#endif

   USE module_saturation        , ONLY: SATQM

   USE module_ra_zenith         , ONLY: ISNEWDAY

   USE module_model_constants   , ONLY: cp, rovcp => rcp, stbolt, karman, xlf, rhowater

   USE module_sf_som_constant   , ONLY: tsice, thsice, snsice

#ifdef XOML
!^^^^^^^^^^^ FOR XOML ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   USE module_xoml_driver       , ONLY: XOML_driver,XOML_init
!$$$$$$$$$$$ FOR XOML $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#endif
!----------------------------------------------------------------------

!------------Yufeng----------
use module_date_time,         ONLY: get_julgmt
use BioCroX
use biocro_parameters
use module_domain
USE module_sf_cssp_constant, ONLY: denliq,denice
!---------------------------


   IMPLICIT NONE

   PUBLIC surface_init
   PUBLIC surface_driver
#ifdef DoVAST
   PUBLIC surface_finish
#endif

CONTAINS

!--------------------------------------------------------------------
   SUBROUTINE surface_driver(                                       &
              itimestep,DT     ,cdate   ,STEPSF  ,fdtcza  ,         &
              do_sfpr_avg      ,dt_sfpr_avg      ,scheme_soilwater ,&
              biocro_switch,  &
              sf_sfclay_physics,sf_land_physics  ,sf_ocean_physics ,&
              bl_pbl_physics   ,ra_lw_physics    ,DX      ,DY      ,&
              num_soil_layers  ,IFSNOW  ,ISFFLX  ,XICE    ,XLAND   ,&
              XLONG   ,XLAT    ,HT      ,VEGFRA  ,IVGTYP  ,ISLTYP  ,&
              GLWd    ,GSWd    ,GSW     ,                           &
              p_phy   ,u_phy   ,v_phy   ,t_phy   ,th_phy  ,         &
              rho     ,p8w     ,dz8w    ,z       ,pi_phy  ,         &
              qv_curr ,qc_curr ,qi_curr ,                           &
              PRECCV  ,PRECNCV ,SNOWCV  ,SNOWNCV ,EMISS   ,RMOL    ,&
              MAVAIL  ,THC     ,ZS      ,DZS     ,CAPG    ,REGIME  ,&
              HFX     ,QFX     ,TSK     ,ZNT     ,UST     ,QST     ,&
              SMOIS   ,TSLB    ,GRDFLX  ,SNOW    ,SNOWH   ,SNOWF   ,&
              POTEVP  ,SFCEVP  ,SFRUNOFF,UDRUNOFF,ACSNOW  ,ACSNOM  ,&
              SNOPCX  ,ALBEDO  ,CANWAT  ,Q2M     ,T2M     ,TMN     ,&
              BR      ,PSIM    ,PSIH    ,GZ1OZ0  ,WSPD    ,TH2     ,&
              LOWLYR  ,PBLH    ,AKHS    ,AKMS    ,CHKLOWQ ,QSFC    ,&
              QSG     ,QVG     ,QCG     ,                           &
              THZ0    ,QZ0     ,UZ0     ,VZ0     ,CT      ,TKE_PBL ,&
              U10     ,V10     ,LH      ,SH2O    ,                  &
              SMSTAV  ,SMSTOT  ,SFCEXC  ,PORSL   ,                  &
              ALBBCK  ,EMBCK   ,SNOALB  ,SHDMIN  ,SHDMAX  ,Z0      ,&
              SOILT1  ,TSNAV   ,SMFR3D  ,KEEPFR3DFLAG     ,         &
              ISFTCFLX,USTM    ,F       ,scwater ,                  &
              ISWATER ,ISURBAN ,ISICE   ,MMINLU  ,                  &
              RDLAI2D ,USEMONALB        ,                           &
              FRACTIONAL_SEAICE,SFENTH  ,                           &
              !--------for LSM---------------------------------------
              SNOTIME ,                                             &
              !--------for PXLSM-------------------------------------
              vegf_px ,                                             &
              landusef,soilctop,soilcbot,nlcat   ,nscat   ,         &
              pxlsm_smois_init ,                                    &
#ifdef PX_NDG
              pxlsm_soil_nudge ,anal_interval    ,                  &
              t2_ndg_old, t2_ndg_new,    q2_ndg_old, q2_ndg_new,    &
              sn_ndg_old, sn_ndg_new,    t2obs   ,q2obs   ,         &
#endif
              !--------for CSSP--------------------------------------
              charn,                                                &
              nl_soil ,nl_snow ,sw_band ,n_swsfc ,                  &
              SC_WATER,SC_LANDU,MLU,MLU_PER,BEDROCK ,               & !Yufeng MLU
              fr_sand,fr_clay,                                      &
              xbfdir  ,xwtslope,xuptkw  ,                           &
              xfrost_table,xzwt_perched ,xqflx_sub_snow , &
#ifdef MHydro
              xfdshc  ,                                             &
#endif
#ifdef GenSWC
              xgswcA  ,xgswcB  ,xgswcC  ,xgswcD  ,xvwc10  ,xvwchs  ,&
#endif
              xbsw    ,                                             &
              xcsol   ,xporsl  ,xphi0   ,xdkmg   ,xdksatu ,         &
              xdkdry  ,xhksati ,xhksm   ,xrootfr ,xwdry   ,         &
#ifdef DoVAST
              xswdepth,xrtws   ,                                    &
#endif
#ifdef OroRad
              xfcor   ,xfsky   ,rad_oro ,                           &
#endif
#ifdef LXZALB
              xsalf   ,                                             &
#else
              xalbsol ,                                             &
#endif

#ifdef NewLAKE
              xfrlak, xdplak, xcoszn, xzlake, xdzlak,               &
              xstke1, xdpice, xt_ice, xlkice, xlktmp,               &
              lakemodeltype,&
#endif


              xswi    ,xrainc  ,xrainnc ,xsnowc  ,xsnownc ,         &
              xtss    ,xwliq   ,xwice   ,xz      ,xdz     ,         &
              xtg     ,xtlsun  ,xtlsha  ,                           &
              xfveg   ,xlai    ,xsai    ,xgreen  ,                  &
              xsag    ,xsigf   ,                                    &
              xalb    ,xssun   ,xssha   ,                           &
              xfsun   ,xthermk ,xextkb  ,xextkd  ,                  &
              coszan  ,eccfn   ,xalbg   ,xalbv   ,                  &
              xolrg   ,xsabg   ,xsabvg  ,                           &
              xtaux   ,xtauy   ,xfsenl  ,xfevpl  ,                  &
              xetr    ,xfseng  ,xfevpg  ,                           &
              xshfdt  ,xzwt    ,xtsw    ,xwt     ,xqrchrg ,         &
              xrsur   ,xrbas   ,xrdrn   ,xrsat   ,xrnof   ,xfcov   ,&
              xtstar  ,ZOL     ,xram    ,xrah    ,xraq    ,         &
              xfm     ,xfh     ,xfq     ,                           &
              xfh2    ,xfq2    ,xf10    ,xfq10   ,                  &
              xparsun ,xparsha ,xsabvsun,xsabvsha,                  &
              xrst    ,xassim  ,xrespc  ,xxerr   ,xzerr   ,         &
              !--------for surface boundary conditions---------------
              LAIB_1  ,LAIB_2  ,SAIB_1  ,SAIB_2  ,GRNB_1  ,GRNB_2  ,&
              VEGB_1  ,VEGB_2  ,dtveg   ,setveg  ,vegfrq  ,schveg  ,&
              ALBB_1  ,ALBB_2  ,setalb  ,                           &
              SSTB_1  ,SSTB_2  ,dtsst   ,setsst  ,sstfrq  ,schsst  ,&
              !--------optional SSTSKIN------------------------------
              skin_sst,dtwsk   ,                                    &
              !--------optional MYNNSFC------------------------------
              tsq,qsq,cov,cqc,                                      &
#ifdef CalUCM
              !--------optional URBAN -------------------------------
              sf_urban_physics,                                     &
              frc_urb2d,utype_urb2d,                                &
              num_roof_layers,num_wall_layers,num_road_layers,      &
              dzr,dzb,dzg,                                          &
              tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,         &
              trl_urb3d,tbl_urb3d,tgl_urb3d,                        &
              xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,          &
              num_urban_layers,                                     &
              trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,              &
              sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,            &
              a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,              &
              b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,              &
              dlg_bep,dl_u_bep,sf_bep,vl_bep,                       &
              ts_urb2d,qs_urb2d,t2_urb2d,q2_urb2d,                  &
              u10_urb2d,v10_urb2d,uc_urb2d,ust_urb2d,               &
              sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,                   &
              alb_urb2d,ems_urb2d,                                  &
              psim_urb2d,psih_urb2d,gz1oz0_urb2d,akms_urb2d,        &
#endif
              !--------optional OML ---------------------------------
              omlcall,oml_gamma,tml,t0ml,hml,h0ml,huml,hvml,        &
              !--------optional observation nudging------------------
              uratx   ,vratx   ,tratx   ,                           &
              !--------optional adaptive time step-------------------
              sfdt    ,curr_secs, adapt_step_flag,                  &
              !--------optional wind perturbation--------------------
              u_frame ,v_frame ,                                    &
              !--------optional diagnostic output--------------------
              ck,cka,cd,cda, t2min,t2max,                           &
              aswd,aswu,alwd,alwu,ahfx,alfx,agfx,                   &
              !------------------------------------------------------
#ifdef XOML
!^^^^^^^^^^^ FOR XOML ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              oml_immpene, oml_jrvtype, oml_rlaxopt, &
              oml_kbkspar, oml_rlaxtau, oml_isrpene,    &
              oml_dt, oml_tskin,couple_oml, &
              nl_oml         ,              &
              omlzly         ,omlzlv    ,   &
              omldz          ,sst_nstep ,   &
              t_oml        ,s_oml   ,   &
              u_oml        ,v_oml   ,   &
              e_oml                 ,   &
              qcl_oml        ,buo_oml   ,   &
              gbo_oml        ,eps_oml   ,   &
              vtz_oml        ,ksz_oml   ,   &
              kqz_oml        ,rit_oml   ,   &
              eln_oml                   ,   &
              hobl           ,              &
              dtw1           , sstsk    ,   &
              prevdsst       ,tempsst   ,   &
              obsprevdsst    ,obstempsst,   &
!$$$$$$$$$$$ FOR XOML $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#endif

!Yufeng------------
              TTc,Grain,Leaf, Stem,Root,Rhizome,&
              LeafN,Sp,vmax,alpha,stomataWS,cws,biolai,&
              bioev,biotr,xcoszn,&
!-------------
              ids,ide ,jds,jde ,kds,kde ,&
              ims,ime ,jms,jme ,kms,kme ,&
              i_start,i_end,j_start,j_end,kts,kte,num_tiles)
!
!======================================================================
!
!  This driver calls subroutines for surface parameterizations
!
!             External Surface Boundary Conditions (EXT_SBC)
!
!             Common Land Model+3D Hydrology (CSSP)
!
!             NOAH Land Surface Model (LSM) ---------+
!                                                    |
!             RUC  Land Surface Model (LSMRUC) ------+
!                                                    | + Surface Layer
!             NMM  Land Surface Model (NMMLSM) ------+     (SFCLAY)
!                                                    |
!             PX   Land Surface Model (PXLSM) -------+
!                                                    |
!             Slab Model (SLAB) ---------------------+
!
!             Urban Canopy Model (URBAN) coupled with -----------+
!                                                                |
!                               (CSSP LSM LSMRUC NMMLSM PXLSM) --+
!
!             Simple Ocean Model (SOM)
!
!
! Authors: Xin-Zhong Liang @ ISWS/UIUC, xliang@uiuc.edu, October 2004
! Revised: Xin-Zhong Liang, August 2008
!
!======================================================================
!
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level (at the surface)
!
!======================================================================
!
   INTEGER, INTENT(IN) ::  &
         ids,ide,          &! start/end index for i in domain
         jds,jde,          &! start/end index for j in domain
         kds,kde,          &! start/end index for k in domain
         ims,ime,          &! start/end index for i in memory
         jms,jme,          &! start/end index for j in memory
         kms,kme,          &! start/end index for k in memory
         kts,kte            ! start/end index for k in tile

   INTEGER, INTENT(IN) ::  &
         num_tiles          ! number of tiles
   INTEGER, INTENT(IN), DIMENSION(num_tiles) :: &
         i_start,i_end,    &! start/end index for i in tile
         j_start,j_end      ! start/end index for j in tile

   INTEGER, INTENT(IN) ::  &
         scheme_soilwater ,&! scheme id for soilwater calculation
         biocro_switch,    &! Yufeng
         sf_sfclay_physics,&! scheme id for surface layer
         sf_land_physics,  &! scheme id for land process model
         sf_ocean_physics, &! scheme id for ocean process model
         bl_pbl_physics,   &! scheme id for planetary boundary layer
         ra_lw_physics      ! scheme id for longwave radiation

   REAL, INTENT(IN) ::  &
         charn        ! charn coef
   LOGICAL, INTENT(IN) ::  &
         do_sfpr_avg        ! true to do precip averaging for sfc
   REAL, INTENT(INOUT) ::  &
         dt_sfpr_avg        ! time period length in last precip averaging for sfc [s]

   INTEGER, INTENT(IN) ::  &
         nl_soil,          &! number of soil layers (>0) (CSSP)
         nl_snow,          &! max number of snow layers (<0) (CSSP)
         sw_band,          &! number of solar bands=4(visible,near-infrared;direct,diffuse)
         n_swsfc,          &! number of sfc sw spectral bands
         num_soil_layers,  &! number of soil layer (LSM)
         STEPSF             ! number of dynamic steps between SF calls

   LOGICAL, INTENT(IN), OPTIONAL :: &
         adapt_step_flag    ! flag to use adaptive time step
   REAL, INTENT(IN), OPTIONAL :: &
         sfdt,             &! time interval [minute] calling this surface_driver
         curr_secs          ! current time  [second]

   INTEGER, INTENT(IN) ::  &
         itimestep,        &! model time step index
         IFSNOW,           &! 1 to accounte for snow-cover effects (SLAB)
         ISFFLX,           &! 0 to turn off surface flux calculation (Not CSSP)
         ISFTCFLX,         &! 1 to use (AHW) surface roughness over water
         FRACTIONAL_SEAICE,&! 1 to use fractional sea-ice
         ISWATER,          &! land use id for water category
         ISURBAN,          &! land use id for major urban category
         ISICE              ! land use id for land-ice category

   LOGICAL, INTENT(IN) ::  &
         RDLAI2D,          &! true to read in data for xlai
         USEMONALB          ! true to use monthly albedo

   CHARACTER(LEN=*), INTENT(IN) :: &
         MMINLU             ! land use classification type
 
   REAL, INTENT(IN)    ::  &
         DX,               &! x-dir grid spacing [m]
         DY,               &! y-dir grid spacing [m]
         DT,               &! dynamic model time step interval [seconds]
         dtsst,            &! time lag [seconds] since SSTB_1 is read
         dtveg,            &! time lag [seconds] since VEGB_1 is read
         fdtcza             ! fraction of dtimesf shifted for cosza
                            ! 0-1: dtcza=dtimesf*fdtcza (0.5=mid of dtimesf)
                            !  >1: cosza averaged during dtimesf

   REAL, INTENT(IN)    ::  &
         SFENTH             ! enthalpy flux factor for zot over ocean
                            ! 0: zot via charnock; >0: zot enhanced if wind>15m/s

   REAL, INTENT(IN), DIMENSION(1:num_soil_layers) :: &
         ZS,               &! depths of centers of soil layers [m] (LSM,SLAB)
         DZS                ! thickneses of soil layers [m] (LSM,SLAB)

   REAL, INTENT(IN), DIMENSION(ims:ime,kms:kme,jms:jme) :: &
         p_phy,            &! pressure [Pa]
         u_phy,            &! u-velocity on theta points [m/s]
         v_phy,            &! v-velocity on theta points [m/s]
         t_phy,            &! temperature [K]
         th_phy,           &! potential temperature [K]
         pi_phy,           &! exner function = t_phy/th_phy (dimensionless)
         rho,              &! dry air density [kg/m^3]
         z,                &! height above sea level [m]
         p8w,              &! pressure at full levels [Pa]
         dz8w               ! dz between full levels [m]

   REAL, INTENT(IN), DIMENSION(ims:ime,kms:kme,jms:jme) :: &
         qv_curr,          &! water vapor mixing ratio [kg/kg]
         qc_curr,          &! cloud water mixing ratio [kg/kg]
         qi_curr            ! cloud ice mixing ratio [kg/kg]

   INTEGER, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         IVGTYP,           &! vegetation type (LSM)
         ISLTYP             ! soil type (LSM)

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         MAVAIL,           &! surface moisture availability [0-1]   ! LSMRUC (IO) else (I)
         XLAND,            &! land(1)/water(2) flag (Not CSSP)
         XICE,             &! sea ice (1) flag (LSM)
         TMN                ! soil temperature at lower boundary [K] (LSM,SLAB)
                            ! bulk temperature with no diurnal variation (SSTSKIN)

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         THC                ! thermal inertia [Cal/cm/K/s^0.5] (SLAB)

   LOGICAL, INTENT(IN) ::  &
         setsst,           &! true for updating sea surface temperature
         setveg,           &! true for updating vegetation cover
         setalb             ! true for updating surface albedos

   INTEGER, INTENT(IN) ::  &
         schsst,           &! scheme id for sea surface variation:
                            ! 1=linear-in-time interpolation from wrfsst input
                            ! 2=dynamic-relaxation of model prediction with wrfsst
                            ! 3=dynamic ocean model
         schveg             ! scheme id for vegetation variation:
                            ! 1=linear-in-time interpolation from wrfveg input
                            ! 2=empirical vegetation changes with temperature
                            ! 3=Ecodynamic model

   REAL, INTENT(IN)    ::  &
         sstfrq,           &! sst bc input frequency [seconds]
         vegfrq             ! veg bc input frequency [seconds]

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         VEGB_1,           &! fraction of vegetation cover [-]            at time level 1 (LSM)
         VEGB_2,           &! fraction of vegetation cover [-]            at time level 2 (LSM)
         LAIB_1,           &! leaf area index                             at time level 1 (CSSP)
         LAIB_2,           &! leaf area index                             at time level 2 (CSSP)
         SAIB_1,           &! stem area index                             at time level 1 (CSSP)
         SAIB_2,           &! stem area index                             at time level 2 (CSSP)
         GRNB_1,           &! greenness (fraction of green not dead leaf) at time level 1 (CSSP)
         GRNB_2,           &! greenness (fraction of green not dead leaf) at time level 2 (CSSP)
         ALBB_1,           &! surface albedo (background or snow-free)    at time level 1 (LSM)
         ALBB_2             ! surface albedo (background or snow-free)    at time level 2 (LSM)

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         SSTB_1,           &! sea surface temperature [K]                 at time level 1
         SSTB_2             ! sea surface temperature [K]                 at time level 2

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         PBLH,             &! PBL height [m]
         GLWd,             &! downwelling longwave  flux at the ground surface [W/m^2]
         GSWd               ! downwelling shortwave flux at the ground surface [W/m^2]

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         GSW,              &! net shortwave flux (absorbed) at the ground surface [W/m^2]
         TSK,              &! surface skin temperature [K]
         HFX,              &! upward sensible heat flux at the surface [W/m^2]
         LH,               &! upward latent heat flux at the surface [W/m^2]
         QFX,              &! upward moisture flux at the surface [kg/m^2/s]
         ZNT,              &! roughness length [m]
         UST,              &! u* in similarity theory [m/s]
         USTM,             &! u* in similarity theory [m/s] without vconv correction
         QST,              &! q* (similarity theory) [kg/kg]
         xtstar,           &! t* in similarity theory [K]
         Q2M,              &! water vapor mixing ratio at 2m [kg/kg]
         T2M,              &! temperature at 2m [K]
         TH2,              &! potential temperature at 2m (not 1000hPa) [K]
         CANWAT,           &! canopy water [mm] (lxz: m=>mm)
         VEGFRA,           &! vegetation fraction [0-1]
         SNOWF,            &! snow cover fraction [0-1]
         SNOW,             &! snow water equivalent [mm]
         SNOWH,            &! snow depth [m]
         EMISS,            &! surface emissivity [0-1] (SLAB)
         ALBEDO,           &! surface albedo [0-1]
         GRDFLX,           &! ground heat flux [W/m2]
         CAPG,             &! heat capacity for soil [J/K/m^3]
         REGIME,           &! flag indicating PBL regime (stable, unstable, etc.)
         PSIM,             &! similarity stability function for momentum
         PSIH,             &! similarity stability function for heat
         GZ1OZ0,           &! log(z/z0) where z0 is roughness length
         WSPD,             &! wind speed at lowest model level [m/s]
         BR,               &! bulk Richardson number in surface layer
         ZOL,              &! (lowest model layer height above surface)/(Monin-Obukhov length)
         RMOL,             &! reciprocal of Monin-Obukhov length
         SNOTIME,          &! number of timesteps since last snowfall [-]
         SFRUNOFF,         &! accumulated surface runoff [mm]
         UDRUNOFF,         &! accumulated underground runoff [mm]
         SFCEVP,           &! accumulated surface evaporation [mm]
         ACSNOW,           &! accumulated snowfall [mm]
         ACSNOM,           &! accumulated melted snow [mm]
         SNOPCX,           &! accumulated snow phase change heat flux [W/m^2]
         POTEVP             ! accumulated potential evaporation [mm]

   INTEGER, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         scwater            ! surface category of water characteristics by EXT_SBC
                            ! 1=urban and built-up, 2=soil, 3=wetland, 4=land ice
                            ! 5=shallow lake, 6=deep lake,  7=sea ice, 8=ocean
!----------------------------------------------------------------------
!--MYJ VARIABLEs
!----------------------------------------------------------------------
   INTEGER, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         LOWLYR             ! index of lowest model layer above ground
   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         THZ0,             &! potential temperature at roughness length [K]
         UZ0,              &! u wind component at roughness length [m/s]
         VZ0,              &! v wind component at roughness length [m/s]
         QZ0,              &! specific humidity at roughness length [kg/kg]
         QSFC,             &! specific humidity at lower boundary [kg/kg]
         QSG,              &! saturated mixing ratio at the ground surface [kg/kg]      !-xyuan
         QVG,              &! water vapor mixing ratio at the ground surface [kg/kg]    !-xyuan
         QCG,              &! cloud water mixing ratio at the ground surface [kg/kg]    !-xyuan
         AKHS,             &! sfc exchange coefficient of heat/moisture
         AKMS,             &! sfc exchange coefficient of momentum
         CHKLOWQ            ! check for low moisture availability
   REAL, INTENT(INOUT), DIMENSION(ims:ime,kms:kme,jms:jme) :: &
         TKE_PBL            ! turbulence kinetic energy from PBL [m^2/s^2]
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         CT,               &! countergradient term [K]
       ! TSHLTR,           &! diagnostic 2-m theta [K] = TH2
       ! QSHLTR,           &! diagnostic 2-m specific humidity [kg/kg] = Q2M
         U10,              &! diagnostic 10-m u component [m/s]
         V10                ! diagnostic 10-m v component [m/s]
!----------------------------------------------------------------------
!--MYNN VARIABLEs
!----------------------------------------------------------------------
   REAL, INTENT(IN), DIMENSION(ims:ime,kms:kme,jms:jme), OPTIONAL :: &
         tsq,              &! variance of Theta_l [K^2]
         qsq,              &! variance of Q_w [-]
         cov                ! covariance of Theta_l and Q_w [K]
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         cqc                ! surface exchange coefficient for cloud water [m/s]
!----------------------------------------------------------------------
!--LSM/SLAB VARIABLEs
!----------------------------------------------------------------------
   REAL, INTENT(INOUT), DIMENSION(ims:ime,1:num_soil_layers,jms:jme) :: &
         TSLB,             &! soil temperature [K] (LSM,SLAB)
         SMOIS,            &! volumetric soil moisture [m3/m3] (LSM)
         SH2O               ! unfrozen volumetric soil moisture [m3/m3] (LSM)
   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         ALBBCK,           &! background surface albedo [-] (LSM)
         Z0                 ! background vroughness length [m] (LSM)
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         EMBCK,            &! background surface emissivity [-] (LSM)
         SNOALB,           &! upper bound on maximum albedo over deep snow [-] (LSM)
         SHDMIN,           &! minimum areal coverage of annual green vegetation [-] (LSM)
         SHDMAX             ! maximum areal coverage of annual green vegetation [-] (LSM)
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         PORSL,            &! soil porosity or maximum volumetric moisture content [m3/m3] (LSM)
         SFCEXC,           &! surface exchange coefficient [m/s]
         SMSTOT,           &! column total (liquid+ice) soil moisture content [mm]
         SMSTAV             ! soil moisture availability for evapotranspiration [-]
                            !=ratio of volumetric soil moisture at wilting vs saturation
!----------------------------------------------------------------------
!--RUC VARIABLEs
!----------------------------------------------------------------------
   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         SOILT1,           &! snow temperature at 7.5 cm depth [K]
         TSNAV              ! average temperature of snow pack [C]
   REAL, INTENT(INOUT), DIMENSION(ims:ime,1:num_soil_layers,jms:jme) :: &
         SMFR3D,           &! soil ice
         KEEPFR3DFLAG       ! frozen soil flag [0,1]
!----------------------------------------------------------------------
!--PX VARIABLEs
!----------------------------------------------------------------------
   INTEGER, INTENT(IN) ::  &
         nlcat,            &! number of landuse categories
         nscat,            &! number of soil categories
         pxlsm_smois_init   ! flag to intialize deep soil moisture to a value derived from moisture availiability
                            ! this is typically done for the first simulation
   REAL, INTENT(IN), DIMENSION(ims:ime,1:nlcat,jms:jme) :: &
         landusef           ! landuse fraction [-]
   REAL, INTENT(IN), DIMENSION(ims:ime,1:nscat,jms:jme) :: &
         soilctop,         &!top soil fraction
         soilcbot           ! bottom soil fraction
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
       ! ra,               &! aerodynamic resistance [s/m] = xram
       ! rs,               &! stomatal resistance [s/m] = xrst
         vegf_px            ! vegetation fraction [-]
#ifdef PX_NDG
   INTEGER, INTENT(IN) ::  &
         pxlsm_soil_nudge, &! flag to use soil moisture and temperature nudging
         anal_interval      ! interval of analyses for nudging soil moisture and temperature [s]
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         t2_ndg_old,       &! analysis temperature  [K]     prior to current time
         t2_ndg_new,       &! analysis temperature  [K]     ahead of current time
         q2_ndg_old,       &! analysis mixing ratio [kg/kg] prior to current time
         q2_ndg_new,       &! analysis mixing ratio [kg/kg] ahead of current time
         sn_ndg_old,       &! analysis snow water   [mm]    prior to current time
         sn_ndg_new         ! analysis snow water   [mm]    ahead of current time
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         t2obs,            &! analysis temperature interpolated from prior and next in time analysese
         q2obs              ! analysis moisture    interpolated from prior and next in time analysese
#endif
!----------------------------------------------------------------------
!--CSSP VARIABLEs
!----------------------------------------------------------------------
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         XLONG,            &! longitude [deg]
         XLAT,             &! latitude [deg]
         F,                &! coriorlis parameter [1/s]
         HT,               &! surface elevation [m]
         PRECCV,           &! convective total precipitation rate in last dynamic step [mm/s]
         PRECNCV,          &! resolved   total precipitation rate in last dynamic step [mm/s]
         SNOWCV,           &! convective snowfall rate in last dynamic step [mm/s]
         SNOWNCV            ! resolved   snowfall rate in last dynamic step [mm/s]

   CHARACTER(LEN=*), INTENT(IN) :: &
         cdate              ! current model date ('YYYY-MM-DD HH:MM:%s.ffff')

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         SC_WATER           ! surface category of water characteristics
                            ! 1=urban and built-up, 2=soil, 3=wetland, 4=land ice
                            ! 5=shallow lake, 6=deep lake,  7=sea ice, 8=ocean
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         SC_LANDU,     &           ! surface category of land use (vegetation)
         MLU   !Yufeng Mariginal Land Use (ha)
        ! MLU_PER
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
        BEDROCK       ! soil depth to bedrock [m]
   REAL, INTENT(INOUT), DIMENSION(ims:ime,sw_band,jms:jme) :: &
         xswi               ! solar insolation onto the surface [W/m2]
                            ! (1->4): (visible<0.7um, near-infrared) {direct, diffuse}

   REAL, INTENT(IN), DIMENSION(ims:ime,1:nl_soil,jms:jme) :: &
         xcsol,            &! soil heat capacity of solids [J/m3/K]
         xporsl,           &! soil fraction that is voids [-]
         xwdry,            &! soil dry or minimum residual water content [-]
         xphi0,            &! soil minimum suction [mm]
#ifdef GenSWC
         xgswcA,           &! GenSWC coef A [mm]
         xgswcB,           &! GenSWC coef B [-]
         xgswcC,           &! GenSWC coef C [-]
         xgswcD,           &! GenSWC coef D [-]
         xvwc10,           &! volumetric water content (vwc) at 10 kPa [-]
         xvwchs,           &! volumetric water content (vwc) at hygroscopic [-]
#endif
         xbsw,             &! soil clapp and hornbereger (b) parameter [-]
         xdkmg,            &! soil thermal conductivity of minerals   [W/m/K]
         xdksatu,          &! soil thermal conductivity at saturation [W/m/K]
         xdkdry,           &! soil thermal conductivity when dry      [W/m/K]
         xhksati,          &! soil hydraulic conductivity at saturation [mm h2o/s]
         xhksm,            &! soil hydraulic conductivity at saturation averaged in a layer [mm h2o/s]
         xrootfr            ! soil fraction of vegetation root [-]

   INTEGER, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         xbfdir             ! baseflow direction
   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme):: xqflx_sub_snow  !sublimation rate from snow pack (mm H2O /s) [+]
   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme):: xfrost_table    ! frost table depth (m)
   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme):: xzwt_perched    ! perched water table depth (m)

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         xcoszn,           &! Yufeng: NewLAKE turned off, so we need to define xcoszn here
         xfveg,            &! fraction of vegetation cover [-]
#ifdef MHydro
         xfdshc,           &! decay constant of saturated hydraulic conductivity [1/m]
#endif
         xwtslope           ! water table slope [m/m]

#ifdef OroRad
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         xfsky,            &! sky view factor for orographic radiative effect [-]
         xfcor              ! correction factor for orographic solar slope-aspect-shadowing effects [-]
                            ! ** make sure that xfsky=1 & xfcor=1 if no orographic effect
   INTEGER, INTENT(IN) :: &
         rad_oro            ! =0 if no orography radiative effect
#endif

#ifdef LXZALB
   REAL, INTENT(IN), DIMENSION(ims:ime,sw_band,jms:jme) :: &
         xsalf              ! soil albedo localization factor [-]
                            ! (visible,near-infrared;direct,diffuse)
#else
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         xalbsol            ! soil albedo (visible beam) at saturation [-]
#endif

#ifdef NewLAKE
    integer, intent(in) :: &!
        lakemodeltype
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        xfrlak,             &! lake area fraction [%]
        xdplak,             &! depth of lake [m]
        xcoszn               ! cosine of solar zenith angle [-]

    !real, intent(in), dimension(ims:ime,1:nl_lake,jms:jme) :: &
    real, intent(in), dimension(ims:ime,1:nl_soil,jms:jme) :: &
        xzlake,             &! depth of lake layer [m]
        xdzlak               ! thickness of lake layer [m]

    real, intent(inout), dimension(ims:ime,jms:jme) :: &
        xstke1,             &! first lake layer thermal conductivity [W/m/K]
        xdpice,             &! lake ice thicknes [m]
        xt_ice               ! lake ice temperature [K], only for GEREL

    !real, intent(inout), dimension(ims:ime,1:nl_lake,jms:jme) :: &
    real, intent(inout), dimension(ims:ime,1:nl_soil,jms:jme) :: &
        xlkice,             &! ice mass fraction [-]
        xlktmp               ! lake temperature [kelvin]
#endif

   REAL, INTENT(INOUT), DIMENSION(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &
         xtss,             &! land layer temperature [K]
         xwliq,            &! land layer liquid water [kg/m2]
         xwice,            &! land layer ice lens [kg/m2]
         xz,               &! land layer node depth [m]
         xdz                ! land layer thickness [m]

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         xtg,              &! ground surface temperature [K]
         xtlsun,           &! sunlit leaf temperature [K]
         xtlsha,           &! shaded leaf temperature [K]
       ! xldew,            &! depth of water on foliage [mm] = CANWAT
         xsag,             &! dimensionless snow age [-]
       ! xscv,             &! snow cover, water equivalent [mm] = SNOW
       ! xsnowdp,          &! snow depth [m] = SNOWH
       ! xfsno,            &! fraction of snow cover on ground = SNOWF
         xsigf,            &! fraction of veg. cover free of snow [-]
         xgreen,           &! greenness (fraction of green not dead leaf)
         xlai,             &! leaf area index
         xsai               ! stem area index

   REAL, INTENT(INOUT), DIMENSION(ims:ime,sw_band,jms:jme) :: &
         xalb,             &! total averaged albedo [-]
         xssun,            &! sunlit canopy absorption of solar radiation [-]
         xssha              ! shaded canopy absorption of solar radiation [-]
                            ! (1->4): (visible<0.7um, near-infrared) {direct, diffuse}

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         xfsun,            &! fraction of sunlit canopy
         xthermk,          &! canopy gap fraction for infrared radiation [-]
                            ! fraction of LW pass through canopy versus directly to ground
         xextkb,           &! (k, g(mu)/mu) direct solar extinction coefficient [-]
         xextkd,           &! diffuse and scattered diffuse PAR extinction coefficient [-]
         xwt,              &! total water storage (unsaturated soil water + groundwater) (mm)
         xtsw               ! total soil water below the water table [mm] !-choi2

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         coszan             ! cosine of solar zenith angle at next CSSP step
   REAL, INTENT(OUT) ::    &
         eccfn              ! earth-sun distance factor <(1/r)^2> at next CSSP step

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         xrainc,           &! averaged convective rainfall rate within this CSSP step [mm/s]
         xrainnc,          &! averaged resolved   rainfall rate within this CSSP step [mm/s]
         xsnowc,           &! averaged convective snowfall rate within this CSSP step [mm/s]
         xsnownc,          &! averaged resolved   snowfall rate within this CSSP step [mm/s]
         xxerr,            &! accumulated water  balance error [mm]
         xzwt,             &! water talble depth [m]
         xzerr              ! accumulated energy balance error [W/m2]

   REAL, INTENT(OUT), DIMENSION(ims:ime,sw_band,jms:jme) :: &
         xalbg,            &! ground albedo [-]
         xalbv              ! vegetation albedo [-]
                            ! (1->4): (visible<0.7um, near-infrared) {direct, diffuse}

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         xtaux,            &! x-wind stress [kg/m/s2]
         xtauy,            &! y-wind stress [kg/m/s2]
       ! xfsena,           &! sensible heat from canopy height to atmosphere [W/m2] = HFX
       ! xfevpa,           &! evapotranspiration from canopy height to atmosphere [mm/s] = QFX
       ! xlfevpa,          &! latent heat flux from canopy height to atmosphere [W/m2] = LH
         xfsenl,           &! sensible heat from leaves [W/m2]
         xfevpl,           &! evaporation+transpiration from leaves [mm/s]
         xetr,             &! transpiration rate [mm/s]
         xfseng,           &! sensible heat flux from ground [W/m2]
         xfevpg,           &! evaporation heat flux from ground [mm/s]
         xolrg,            &! outgoing long-wave radiation from ground+canopy [W/m2]
         xshfdt,           &! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]
       ! xemis,            &! averaged bulk surface emissivity = EMISS
       ! xfgrnd,           &! ground heat flux [W/m2] = GRDFLX
       ! xtavg,            &! averaged bulk surface temperature [K] = xtg or xtrad
       ! xqref,            &! 2 m height specific humidity [kg/kg] = Q2M
         xrsur,            &! surface runoff [mm h2o/s]
         xrbas,            &! baseflow runoff [mm h2o/s]
         xrdrn,            &! bottom drainage runoff [mm h2o/s]
         xrsat,            &! excess saturation runoff [mm h2o/s]
         xrnof,            &! total runoff [mm h2o/s]
         xqrchrg,          &! grounfwater recharge rate [mm h2o/s]
         xfcov,            &! fractional area with the water table at surface
       ! xsm,              &! snowmelt [mm h2o] = SNOMLT
         xsabg,            &! solar absorbed by ground  [W/m2]
         xsabvg,           &! solar absorbed by ground + vegetation [W/m2] = GSW w/o orographic effect
       ! xz0ma,            &! averaged aerodynamic roughness length [m] = ZNT
       ! xustar,           &! u* in similarity theory [m/s] = UST
       ! xqstar,           &! q* in similarity theory [kg/kg] = QST
       ! xzol,             &! (lowest model layer height above surface)/(Monin-Obukhov length) = ZOL
       ! xrib,             &! bulk Richardson number in surface layer = BR
       ! xu10,             &! u-velocity at 10m [m/s] = U10
       ! xv10,             &! v-velocity at 10m [m/s] = V10
         xf10,             &! integral of profile function for momentum at 10m
         xfq10,            &! integral of profile function for moisture at 10m
         xfh2,             &! integral of profile function for heat at 2m
         xfq2,             &! integral of profile function for moisture at 2m
         xfm,              &! integral of profile function for momentum
         xfh,              &! integral of profile function for heat
         xfq,              &! integral of profile function for moisture
         xram,             &! aerodynamic momentum resistance [s/m]
         xrah,             &! aerodynamic thermal resistance [s/m]
         xraq,             &! aerodynamic moisture resistance [s/m]
         xrst,             &! canopy stomatal resistance [s/m]
         xassim,           &! canopy assimilation [mol/m2/s]
         xrespc,           &! canopy respiration [mol/m2/s]
         xparsun,          &! PAR by sunlit leaves [W/m2]
         xparsha,          &! PAR by shaded leaves [W/m2]
         xsabvsun,         &! solar absorbed by sunlit vegetation [W/m2]
         xsabvsha           ! solar absorbed by shaded vegetation [W/m2]
   REAL, INTENT(OUT), DIMENSION(ims:ime,1:nl_soil,jms:jme) :: &
         xuptkw             ! root uptake of soil water [mm/s]    !-xyuan

#ifdef DoVAST
!----------------------------------------------------------------------
!--CSS VARIABLEs
!----------------------------------------------------------------------
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         xswdepth           ! surface water depth [mm]

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         xrtws              ! surface water supply rate [mm/s]
#endif

!----------------------------------------------------------------------
!--SSTSKIN VARIABLEs
!----------------------------------------------------------------------
   LOGICAL, OPTIONAL, INTENT(IN) ::  &
         skin_sst           ! true to predict skin sst

   REAL, OPTIONAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         dtwsk              ! warm sea surface temperature departure [C]

#ifdef CalUCM
!----------------------------------------------------------------------
!--UCM VARIABLEs
!----------------------------------------------------------------------
   INTEGER, INTENT(IN) ::  &
         sf_urban_physics   ! 1 or 2 to call UCM for urban effect

   INTEGER, OPTIONAL, INTENT(IN) :: &
         num_roof_layers,  &! number of urban roof layers
         num_wall_layers,  &! number of urban wall layers
         num_road_layers    ! number of urban road layers

   REAL, OPTIONAL, INTENT(IN), DIMENSION(1:num_soil_layers) :: &
         dzr,              &! urban roof layer thickness [cm]
         dzb,              &! urban building wall layer thickness [cm]
         dzg                ! urban ground layer thickness [cm]

   INTEGER, OPTIONAL, INTENT(IN), DIMENSION(ims:ime, jms:jme) :: &
         utype_urb2d        ! urban type [nonurban=0, urban=1, suburban=2, rural=3]

   REAL, OPTIONAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         frc_urb2d          ! urban fraction of the grid [-]

   REAL, OPTIONAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         tr_urb2d,         &! urban roof temperature [K]
         tb_urb2d,         &! urban building wall temperature [K]
         tg_urb2d,         &! urban road temperature [K]
         tc_urb2d,         &! urban-canopy air temperature [K]
         qc_urb2d,         &! urban-canopy air mixing ratio [kg/kg]
         xxxr_urb2d,       &! urban roof   z/Monin-Obukhov-length [-]
         xxxb_urb2d,       &! urban wall   z/Monin-Obukhov-length [-]
         xxxg_urb2d,       &! urban road   z/Monin-Obukhov-length [-]
         xxxc_urb2d         ! urban-canopy z/Monin-Obukhov-length [-]

   REAL, OPTIONAL, INTENT(INOUT), DIMENSION(ims:ime,1:num_soil_layers,jms:jme) :: &
         trl_urb3d,        &! urban roof layer temperature [K]
         tbl_urb3d,        &! urban wall layer temperature [K]
         tgl_urb3d          ! urban road layer temperature [K]

   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         ts_urb2d,         &! urban surface radiative temperature [K]
         qs_urb2d,         &! urban surface humidity [kg/kg]
         t2_urb2d,         &! urban air temperature at 2m [K]
         q2_urb2d,         &! urban air humidity at 2m [kg/kg]
         u10_urb2d,        &! urban u-velocity at 10m [m/s]
         v10_urb2d,        &! urban v-velocity at 10m [m/s]
         uc_urb2d,         &! urban-canopy wind speed [m/s]
         sh_urb2d,         &! urban sensible heat flux [W/m2]
         lh_urb2d,         &! urban latent heat flux [W/m2]
         g_urb2d,          &! urban heat flux into the ground [W/m2]
         rn_urb2d,         &! urban net radiation at the surface [W/m2]
         alb_urb2d,        &! urban effective albedo [-]
         ems_urb2d,        &! urban effective emissivity [-]
         ust_urb2d,        &! urban u* friction velocity [m/s]
         psim_urb2d,       &! urban PSIM
         psih_urb2d,       &! urban PSIH
         gz1oz0_urb2d,     &! urban GZ1OZ0
         akms_urb2d         ! urban AKMS
!
   INTEGER, INTENT(IN) :: num_urban_layers

   REAL, OPTIONAL, INTENT(INOUT), DIMENSION(ims:ime,1:num_urban_layers,jms:jme) :: &
         trb_urb4d        ,&!
         tw1_urb4d        ,&!
         tw2_urb4d        ,&!
         tgb_urb4d        ,&!
         sfw1_urb3d       ,&!
         sfw2_urb3d       ,&!
         sfr_urb3d        ,&!
         sfg_urb3d          !

       ! implicit and explicit components of the source and sink terms at each levels
       ! fluxes can be computed as: fx = a*x + b   example: t_fluxes = a_t * pt + b_t
   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,kms:kme,jms:jme) :: &
         a_u_bep          ,&! implicit component for momemtum in x-dir (center)
         a_v_bep          ,&! implicit component for momemtum in y-dir (center)
         a_t_bep          ,&! implicit component for temperature
         a_q_bep          ,&! implicit component for humidity
         a_e_bep          ,&! implicit component for TKE
         b_u_bep          ,&! explicit component for momemtum in x-dir (center)
         b_v_bep          ,&! explicit component for momemtum in y-dir (center)
         b_t_bep          ,&! explicit component for temperature
         b_q_bep          ,&! explicit component for humidity
         b_e_bep          ,&! explicit component for TKE
         vl_bep           ,&! volume of the urban grid cells
         dlg_bep          ,&! height above ground (l_ground in formula (24) of the BLM paper)
         sf_bep           ,&! surface of the urban grid cells
         dl_u_bep           ! length scale (l_b in formula (22) of the BLM paper)
#endif
!----------------------------------------------------------------------
!--OML VARIABLEs
!----------------------------------------------------------------------
   INTEGER, OPTIONAL, INTENT(IN) :: &
         omlcall            ! 1 to call OML
   REAL, OPTIONAL, INTENT(IN) :: &
         oml_gamma          ! lapse rate below mixed layer in ocean (default 0.14) [K/m]
   REAL, OPTIONAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         t0ml,             &! initial ocean mixed layer temperature [K]
         h0ml               ! initial ocean mixed layer depth [m]
   REAL, OPTIONAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         tml,              &! ocean mixed layer temperature [K]
         hml,              &! ocean mixed layer depth [m]
         huml,             &! ocean mixed layer u component of wind [m/s]
         hvml               ! ocean mixed layer v component of wind [m/s]
!----------------------------------------------------------------------
!-Observation nudging
!----------------------------------------------------------------------
   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         uratx,            &! ratio of u-wind at the lowest layer versus 10m 
         vratx,            &! ratio of v-wind at the lowest layer versus 10m 
         tratx              ! ratio of potential temperature at the lowest layer versus 2m 
!----------------------------------------------------------------------
!-Wind perturbation
!----------------------------------------------------------------------
   REAL, OPTIONAL, INTENT(IN) :: &
         u_frame,          &! idealized u-wind perturbation [m/s]
         v_frame            ! idealized v-wind perturbation [m/s]
!----------------------------------------------------------------------
!-Diagnostic output
!----------------------------------------------------------------------
   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         ck,               &! enthalpy exchange coeff at 10m
         cd,               &! momentum exchange coeff at 10 meters
         cka,              &! enthalpy exchange coeff at the lowest layer
         cda                ! momentum exchange coeff at the lowest layer
   REAL, OPTIONAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         t2min,            &! daily minimum temperature at 2m [K]
         t2max,            &! daily maximum temperature at 2m [K]
         aswd,             &! daily accumulated downwelling sw flux [W/m2]
         aswu,             &! daily accumulated   upwelling sw flux [W/m2]
         alwd,             &! daily accumulated downwelling lw flux [W/m2]
         alwu,             &! daily accumulated   upwelling lw flux [W/m2]
         ahfx,             &! daily accumulated sensible  heat flux [W/m2]
         alfx,             &! daily accumulated latent    heat flux [W/m2]
         agfx               ! daily accumulated ground    heat flux [W/m2]
!^^^^^^^^^^^ FOR XOML ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifdef XOML
  integer  ,intent(in):: oml_immpene, oml_jrvtype , oml_rlaxopt , oml_isrpene
  real     ,intent(in):: oml_kbkspar, oml_rlaxtau 
  real     ,intent(in):: oml_dt
  integer  ,intent(in):: oml_tskin
  integer  ,intent(in):: couple_oml
  integer,intent(in):: nl_oml
  real,dimension(nl_oml),intent(inout):: omlzly ,omlzlv ,omldz
  integer,intent(inout)::sst_nstep
  real,dimension(ims:ime,nl_oml,jms:jme),intent(inout)::  &
        t_oml,s_oml,u_oml,v_oml,e_oml ,                     &
        qcl_oml,buo_oml,gbo_oml,eps_oml,vtz_oml,ksz_oml,              &
        kqz_oml,rit_oml,eln_oml
  real,dimension(ims:ime,       jms:jme),intent(inout)::  &
        hobl   , dtw1, sstsk 
  real,dimension(ims:ime,       jms:jme),intent(inout)::  &
        prevdsst, tempsst, obsprevdsst, obstempsst
#endif
!$$$$$$$$$$$ FOR XOML $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!----------------------------------------------------------------------
!
!  LOCAL  VARIABLEs
!
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: u_phytmp,v_phytmp

   REAL, DIMENSION(ims:ime,jms:jme) :: &
         RAINCV,           &! convective rainfall rate in last dynamic step [mm/s]
         RAINNCV,          &! resolved   rainfall rate in last dynamic step [mm/s]
         PRECSFC,          &! averaged   total pr rate in this surface step [mm/s]
         SNOWSFC,          &! averaged   snowfall rate in this surface step [mm/s]
         FPRFRZN,          &! fraction of frozen from total precipiation [-]
         GSWu               ! upwelling shortwave flux at the ground surface [W/m^2]

   REAL, DIMENSION(ims:ime,jms:jme) :: &
         NOAHRES,          &! residual of surface energy balance equation terms [W/m2]
         SNOMLT,           &! snow melt during this surface step [mm]
         SOILTB,           &! bottom soil temperature [K]
         sicthk             ! sea-ice thickness [m]
!
   REAL, DIMENSION(ims:ime,jms:jme) :: &
         ZASF,             &! lowest model layer height above surface [m]
         PSFC,             &! surface pressure [Pa]
         FLHC,             &! exchange coefficient for heat = RHO*CHS*CP [W/m2/K]
         FLQC,             &! exchange coefficient for moisture = RHO*CQS [kg/m2/s]
         QGH,              &! saturated mixing ratio at the lowest layer [kg/kg]
         CPM,              &! CP*(1+virt*Qv) [J/kg/K]
         CHS,              &! surface exchange coefficient for heat [m/s]
         CHS2,             &! surface exchange coefficient for heat at 2m [m/s]
         CQS2               ! surface exchange coefficient for moisture at 2m [m/s]

   INTEGER, DIMENSION(ims:ime,jms:jme) :: &
         sclandu            ! vegetation type for CSSP

   INTEGER :: &
         icwater            ! starting scwater id to update TSK by EXT_SBC
!
   REAL ::                 &
         dtimesf            ! time interval [seconds] calling this surface_driver

   LOGICAL ::              &
         run_param,        &! true for doing land surface calculation
         dosfclay,         &! true for sfclay calculation
         dosfcflx,         &! true for surface flux calculation
         doocndyn,         &! true for ocean dynamic calculation
         getveg,           &! true for updating vegetation cover
         radiation,        &! true for radiation calculation in SLAB
         dosfcvisc,        &! true for viscous sublayer needed by MYJ PBL scheme
         usemyjsfc,        &! true for using MYJ surface scheme
#ifdef CalUCM
         dourban,          &! true for calling urban canopy model
#endif
         obsnudge,         &! true for saving observation nudging variables
         convtran,         &! true for converting suface transfer functions
         sfcqwcld,         &! true for setting surface exchange coefficient for cloud water
         doregime,         &! true for setting (regime) from bulk Richard number (br)
         psimpsih           ! true for setting (psim,psih,gz1oz0) from (xfm,xfh,znt)

   INTEGER ::              &
#ifdef CalUCM
         ucmdoit,          &! 1 for actually calling UCM
#endif
         ntsflg,           &! 1=do GFDL SLAB
         schscw,           &! scheme id for surface water category
                            ! 1=SC_WATER
                            ! 2=IVGTYP
                            ! 3=XLAND
         lisveg,           &! list of veg variables:
                            ! 1 for out=fveg
                            ! 2 for inp=fveg out=lai,sai=0,green=1
                            ! 3 for inp=fveg out=lai,sai,green=1
                            ! 4 for inp=fveg,out=lai,sai,green
                            ! 5 for out=fveg,lai,sai,green
         lisalb,           &! list of alb variables:
                            ! 1 for out=albbck
         on_lake,          &! 1=update TSK over in-land lake as well
         Ic_water,         &! in-land water category in IVGTYP (USGS=16)
         ic_ocean,         &! ocean category in IVGTYP (USGS=0)
         ic_ice             ! ice category in IVGTYP (USGS=24)

   REAL ::                 &
         fc_land,          &! threshold for XLAND to be a land grid (defult=1.5)
         fc_ice,           &! threshold for fractional cover to be an ice grid (defult=0.5)
         xice_threshold,   &! threshold for sea-ice grid
         f_past, f_this,   &!
         kofm,kof10,ustk,xlfrdt
!
   INTEGER :: i,j,k,l,n,its,ite,jts,jte,itf,jtf,ktf,ij
   INTEGER :: initial
!

!---------------------Yufengs
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   CHARACTER*256  :: outstring1, outstring2
   integer:: juldy,curyr,doy,doy0,doy1,s
   real::gmtim,solar_dir,solar_diff,temp,rh,windspeed,precip,hour
   real,dimension(soilLayers)::water_status,root_distribution,iWatCont,vol_ice,&
                               eff_porosity,fieldc,wiltp 
   type(nitroPstructure)::nitroP
   type(biocro_outputs)::bio_out
   real, DIMENSION(ims:ime,jms:jme) ::TTc,Grain,Leaf,Stem,Root,Rhizome,LeafN,&
                          Sp,vmax,alpha,stomataWS,biolai,MLU_PER,bioev,biotr
   real,dimension(ims:ime,1:soilLayers,jms:jme)::cws,fr_sand,fr_clay
   real :: sand,clay,om_frac,om_tkm,om_watsat,om_b,bsw,sucsat,watsat,tkm
   logical, save :: bio_onoff
!TYPE ( domain ) :: grid
!---------------------

!----------------------------------------------------------------------
!  Turn on all switches
!----------------------------------------------------------------------
!
   IF (.NOT. (sf_land_physics/=0 .OR. sf_ocean_physics/=0) ) RETURN
!
!  ....Modified for adaptive time step.............................
   IF (PRESENT(sfdt) .AND. PRESENT(curr_secs) .AND.         &
       PRESENT(adapt_step_flag) .AND. (adapt_step_flag)) THEN
      run_param = (itimestep == 1) .or. (sfdt == 0) .or.    &
                  (curr_secs+DT >= (INT(curr_secs/(sfdt*60))+1)*sfdt*60)
   ELSE
      run_param = (itimestep == 1) .or. (STEPSF  == 1) .or. &
                  (MOD(itimestep,STEPSF) == 0)
   ENDIF

!  ===================================
   surface_prep: IF ( run_param ) THEN
!  ===================================

!  ...Surface schemes need SF time step for updates
   dtimesf = DT * max(STEPSF,1)
   IF (PRESENT(sfdt)) THEN
      IF (sfdt == 0) THEN
         dtimesf= DT
      ELSEIF (PRESENT(adapt_step_flag) .and. (adapt_step_flag)) THEN
         call wrf_message("WARNING: When using an adaptive time-step, the surface module " // &
                          "time-step should be 0 (i.e., equivalent to model time-step). "  // &
                          "In order to proceed, the surface module time-step will be "     // &
                          "rounded to the nearest minute, possibly inaccuracy.")
         dtimesf= DT * NINT(sfdt*60/DT)
      ENDIF
   ENDIF
!  ................................................................

   dosfclay = .true.
   IF (sf_land_physics==   CSSPSCHEME)  dosfclay = .false.
   IF ((bl_pbl_physics== MYJPBLSCHEME) .OR. &
       (bl_pbl_physics== MYJURBSCHEME) .OR. &
       (bl_pbl_physics==QNSEPBLSCHEME)) dosfclay = .true.
       !-enable {MYJSFC,QNSESFC} to set (THZ0,QZ0,UZ0,VZ0)

   IF (dosfclay) THEN
#ifndef LaxSFC
      IF (bl_pbl_physics   == MYJPBLSCHEME .AND. &
          sf_sfclay_physics/= MYJSFCSCHEME) THEN  ! it is relaxed by dosfcvisc
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/pbl = ', &
                 sf_sfclay_physics,bl_pbl_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (bl_pbl_physics   == MYJURBSCHEME .AND. &
          sf_sfclay_physics/= MYJSFCSCHEME) THEN  ! it is relaxed by dosfcvisc
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/pbl = ', &
                 sf_sfclay_physics,bl_pbl_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (sf_land_physics  == RUCLSMSCHEME .AND. &
          sf_sfclay_physics/= MYJSFCSCHEME) THEN  ! it is relaxed
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/land = ', &
                 sf_sfclay_physics,sf_land_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (sf_land_physics  ==  PXLSMSCHEME .AND. &
          sf_sfclay_physics/=  PXSFCSCHEME) THEN  ! it is relaxed
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/land = ', &
                 sf_sfclay_physics,sf_land_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
#endif
      IF (sf_sfclay_physics==MYNNSFCSCHEME .AND. &
          bl_pbl_physics   /=MYNNPBLSCHEME) THEN  ! it is required to set (tsq,qsq,cov)
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/pbl = ', &
                 sf_sfclay_physics,bl_pbl_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (sf_sfclay_physics==0 .AND. &
         (sf_land_physics  ==   SLABSCHEME .OR. &
          sf_land_physics  ==     GFDLSLAB .OR. &
          sf_land_physics  == NMMLSMSCHEME .OR. &
          sf_land_physics  ==    LSMSCHEME)) THEN
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/land = ', &
                 sf_sfclay_physics,sf_land_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
   ENDIF

   dosfcflx = ISFFLX/=0
   IF (sf_land_physics==CSSPSCHEME) dosfcflx = .true.

   doocndyn = .false.
   IF (sf_land_physics==CSSPSCHEME .AND. &
       sf_ocean_physics/=0) doocndyn = .true.

   getveg = setveg
   IF (sf_land_physics==SLABSCHEME .OR. &
       sf_land_physics==GFDLSLAB) getveg = .false.

   lisveg = 1            ! need VEGFRA only
   lisalb = 1            ! need ALBBCK
   IF (sf_land_physics==CSSPSCHEME) THEN

#ifdef DynFVC
      lisveg = 5         ! +xum to out fvegm lsai, green not used now
#else
      lisveg = 3         ! no need for VEGFRA but inp=fveg & out=lai,sai,green=1
#endif

      lisalb = 0         ! no need for ALBBCK
   ENDIF
   IF (.NOT.setalb) lisalb = 0

   xice_threshold = 0.5
   IF (FRACTIONAL_SEAICE==1) xice_threshold = 0.02

   fc_ice   = xice_threshold
   fc_land  = 1.5
   ic_ice   = ISICE
   ic_water = ISWATER
   ic_ocean = ic_water   ! IVGTYP uses same ISWATER to define all water bodies
   on_lake  = 1  ! update TSK over in-land lake as well as ocean
   schscw   = 3  ! SLAB or PXLSM uses XLAND  to define surface water category
   IF (dosfcflx) THEN
      IF (sf_land_physics==CSSPSCHEME) THEN

         !+xum 
#if (defined NewLAKE || defined OldLAKE)
         on_lake = 0      ! CSSP has a built-in predictive lake model
#else
         on_lake = 1      ! update TSK over in-land lake as well as ocean
#endif
         !=xum

         schscw   = 1    ! CSSP uses SC_WATER to define surface water category
         ic_ocean = 0    ! SC_WATER uses 0   to define ocean
      ELSEIF (sf_land_physics==LSMSCHEME .OR. &
              sf_land_physics==NMMLSMSCHEME .OR. &
              sf_land_physics==RUCLSMSCHEME) THEN
         schscw   = 2    ! LSM uses IVGTYP   to define surface water category
      ENDIF
   ENDIF

   usemyjsfc = sf_sfclay_physics==MYJSFCSCHEME

   dosfcvisc = (ISFFLX/=0) .AND. ( &
               (bl_pbl_physics==MYJPBLSCHEME)  .OR. &
               (bl_pbl_physics==MYJURBSCHEME)  .OR. &
               (bl_pbl_physics==QNSEPBLSCHEME) ) .AND. ( &
               (sf_sfclay_physics/=MYJSFCSCHEME) .AND. &
               (sf_sfclay_physics/=QNSESFCSCHEME) )

#ifdef CalUCM
   ucmdoit = 0
   dourban = ( &
   !--optional io for urban canopy model
      (sf_urban_physics==UCMSCHEME     .AND. num_soil_layers==4               .AND. &
      PRESENT(num_roof_layers)         .AND. num_roof_layers==num_soil_layers .AND. &
      PRESENT(num_wall_layers)         .AND. num_wall_layers==num_soil_layers .AND. &
      PRESENT(num_road_layers)         .AND. num_road_layers==num_soil_layers .AND. &
      PRESENT(dzr)       .AND. PRESENT(dzb)       .AND. PRESENT(dzg)          .AND. &
      PRESENT(tr_urb2d)                .AND. PRESENT(tb_urb2d)                .AND. &
      PRESENT(tg_urb2d)                .AND. PRESENT(tc_urb2d)                .AND. &
      PRESENT(qc_urb2d)                .AND.                                        &
      PRESENT(xxxr_urb2d)              .AND. PRESENT(xxxb_urb2d)              .AND. &
      PRESENT(xxxg_urb2d)              .AND. PRESENT(xxxc_urb2d)              .AND. &
      PRESENT(trl_urb3d)               .AND. PRESENT(tbl_urb3d)               .AND. &
      PRESENT(tgl_urb3d)                                                   )  .OR.  &
   !--optional io for urban multi-layer model
      (sf_urban_physics==BEPSCHEME     .AND. num_urban_layers>=4              .AND. &
      PRESENT(trb_urb4d)               .AND. PRESENT(tw1_urb4d)               .AND. &
      PRESENT(tw2_urb4d)               .AND. PRESENT(tgb_urb4d)               .AND. &
      PRESENT(sfw1_urb3d)              .AND. PRESENT(sfw2_urb3d)              .AND. &
      PRESENT(sfr_urb3d)               .AND. PRESENT(sfg_urb3d)               .AND. &
      PRESENT(a_u_bep)                 .AND. PRESENT(b_u_bep)                 .AND. &
      PRESENT(a_v_bep)                 .AND. PRESENT(b_v_bep)                 .AND. &
      PRESENT(a_t_bep)                 .AND. PRESENT(b_t_bep)                 .AND. &
      PRESENT(a_q_bep)                 .AND. PRESENT(b_q_bep)                 .AND. &
      PRESENT(a_e_bep)                 .AND. PRESENT(b_e_bep)                 .AND. &
      PRESENT(dlg_bep)                 .AND. PRESENT(dl_u_bep)                .AND. &
      PRESENT(sf_bep)                  .AND. PRESENT(vl_bep)               ) ).AND. &
   !--optional io for all urban models
      PRESENT(utype_urb2d)             .AND. PRESENT(frc_urb2d)               .AND. &
      PRESENT(ts_urb2d)                .AND. PRESENT(qs_urb2d)                .AND. &
      PRESENT(t2_urb2d)                .AND. PRESENT(q2_urb2d)                .AND. &
      PRESENT(u10_urb2d)               .AND. PRESENT(v10_urb2d)               .AND. &
      PRESENT(uc_urb2d)                .AND. PRESENT(ust_urb2d)               .AND. &
      PRESENT(sh_urb2d)                .AND. PRESENT(lh_urb2d)                .AND. &
      PRESENT(g_urb2d)                 .AND. PRESENT(rn_urb2d)                .AND. &
      PRESENT(alb_urb2d)               .AND. PRESENT(ems_urb2d)               .AND. &
      PRESENT(psim_urb2d)              .AND. PRESENT(psih_urb2d)              .AND. &
      PRESENT(gz1oz0_urb2d)            .AND. PRESENT(akms_urb2d)
   IF (sf_urban_physics==UCMSCHEME .OR. sf_urban_physics==BEPSCHEME) THEN
      IF (dourban) THEN
         ucmdoit = sf_urban_physics
      ELSE
         print*,'**** Warning: Inappropriate IO => Abort calling UCM'
      ENDIF
   ENDIF
#endif
!  ==================
   ENDIF surface_prep
!  ==================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  DO ALL PROCESSES THAT ARE LOCAL, WITHOUT INFLUENCE OF ADJACENT GRIDS
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!=)Set start and end points for tiles

   
   !+xum initilize local variables
   RAINCV = 0.
   RAINNCV = 0.
!
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij ,its,ite,jts,jte,itf,jtf,ktf, i,j,k,l,n )

   DO ij  = 1 , num_tiles
      its = i_start(ij)
      ite = i_end(ij)
      jts = j_start(ij)
      jte = j_end(ij)
      itf = min(ite,ide-1)
      jtf = min(jte,jde-1)
      ktf = min(kte,kde-1)
!
!----------------------------------------------------------------------
!  Average (rainfall,snowfall) rate within dtimesf at each dynamic step
!----------------------------------------------------------------------
!
      do i = its,itf
      do j = jts,jtf
          RAINCV(i,j) =  PRECCV(i,j) -  SNOWCV(i,j) 
         RAINNCV(i,j) = PRECNCV(i,j) - SNOWNCV(i,j)
      enddo
      enddo

   IF (do_sfpr_avg .and. dt_sfpr_avg>0.) THEN
      f_past = dt_sfpr_avg
      f_this = f_past+dt
      f_past = f_past/f_this
      f_this = 1.-f_past

      !--the result is always the average of the last(dt_sfpr_avg+dt) period 

      do i = its,itf
      do j = jts,jtf
          xrainc(i,j) =  RAINCV(i,j)*f_this +  xrainc(i,j)*f_past
          xsnowc(i,j) =  SNOWCV(i,j)*f_this +  xsnowc(i,j)*f_past
         xrainnc(i,j) = RAINNCV(i,j)*f_this + xrainnc(i,j)*f_past
         xsnownc(i,j) = SNOWNCV(i,j)*f_this + xsnownc(i,j)*f_past
      enddo
      enddo
   ELSE
      do i = its,itf
      do j = jts,jtf
          xrainc(i,j) =  RAINCV(i,j)
          xsnowc(i,j) =  SNOWCV(i,j)
         xrainnc(i,j) = RAINNCV(i,j)
         xsnownc(i,j) = SNOWNCV(i,j)
      enddo
      enddo
   ENDIF

!  ===================
   IF (.NOT.run_param) CYCLE ! next ij
!  ===================

   IF ((lisveg==2) .or. (lisveg==3) .or. (lisveg==4)) THEN
      !----- xfveg is on input as a static field from wrfinput or wrfrst
      !----- it is used only in CSSP with a different definition from VEGFRA
      DO j = jts,jtf
      DO i = its,itf
         VEGFRA(i,j) = xfveg(i,j)
      ENDDO
      ENDDO
   ENDIF
!
!----------------------------------------------------------------------
!  Update external surface boundary conditions (wrfveg,wrfsst data)
!----------------------------------------------------------------------
!
   CALL EXT_SBC(cdate   ,setsst  ,getveg  ,schscw  ,         &
                SC_WATER,IVGTYP  ,XLAND   ,XICE    ,on_lake ,&
                ic_water,ic_ocean,fc_land ,ic_ice  ,fc_ice  ,&
                TSK     ,TMN     ,scwater ,icwater ,         &
                VEGFRA  ,xlai    ,xsai    ,xgreen  ,         &
                tsice   ,thsice  ,snsice  ,SNOWH   ,sicthk  ,&
                SSTB_1  ,SSTB_2  ,dtsst   ,sstfrq  ,schsst  ,&
                LAIB_1  ,LAIB_2  ,SAIB_1  ,SAIB_2  ,         &
                GRNB_1  ,GRNB_2  ,VEGB_1  ,VEGB_2  ,         &
                dtveg   ,vegfrq  ,schveg  ,lisveg  ,         &
                ALBB_1  ,ALBB_2  ,lisalb  ,ALBBCK  ,         &
                MLU_PER, biolai,           & !Yufeng
                ids,ide ,jds,jde ,kds,kde ,&
                ims,ime ,jms,jme ,kms,kme ,&
                its,itf ,jts,jtf ,kts,ktf  )

   IF (setsst .AND. schsst==1 .AND. &
      (sf_land_physics==    LSMSCHEME .OR. &
       sf_land_physics== NMMLSMSCHEME .OR. &
       sf_land_physics== RUCLSMSCHEME .OR. &
       sf_land_physics==  PXLSMSCHEME .OR. &
       sf_land_physics==   SLABSCHEME .OR. &
       sf_land_physics==     GFDLSLAB)) THEN
      DO j = jts,jtf
      DO i = its,itf
         if (scwater(i,j)<icwater) cycle
         TSLB(i,:,j) = TSK(i,j)
            xtg(i,j) = TSK(i,j)
      ENDDO
      ENDDO
   ENDIF
!
!----------------------------------------------------------------------
!  Prepare data for this tile
!----------------------------------------------------------------------
!
   k = kts                        ! at the surface
   DO j = jts,jtf
   DO i = its,itf
      ZASF(i,j) =   z(i,k,j) -  HT(i,j)    ! = 0.5*dz8w(i,k,j)
      PSFC(i,j) = p8w(i,k,j)
      GSWu(i,j) =  GSWd(i,j) - GSW(i,j)
   ENDDO
   ENDDO

#ifdef OroRad
!--Make orographic correction only on net SW as original WRF

   IF (sf_land_physics/=CSSPSCHEME .AND. rad_oro/=0) THEN
      DO j = jts,jtf
      DO i = its,itf
         GSW(i,j) = GSW(i,j) * xfcor(i,j)
      ENDDO
      ENDDO
   ENDIF
#endif

   IF (sw_band==4) THEN
      IF (n_swsfc/=sw_band) THEN                    ! no xswi but having GSWd
         DO j = jts,jtf
         DO i = its,itf
            xswi(i,1,j) = GSWd(i,j)*0.28            ! vis dir
            xswi(i,2,j) = GSWd(i,j)*0.31            ! nir dir
            xswi(i,3,j) = GSWd(i,j)*0.24            ! vis dif
            xswi(i,4,j) = GSWd(i,j)*0.17            ! nir dif
          ! vis/nir solor insolation = 0.47074/0.52926 at the top of atmosphere
          ! but the ratio reversed at the surface due to greater h2o absorption in nir than vis
          ! the above partition factors (sum=1.0) are subjective and can be better estimated
         ENDDO
         ENDDO
      ENDIF
   ENDIF
!
   IF (PRESENT(u_frame) .AND. PRESENT(v_frame)) THEN
      DO j = jts,jtf
      DO k = kts,ktf
      DO i = its,itf
         u_phytmp(i,k,j) = u_phy(i,k,j)+u_frame
         v_phytmp(i,k,j) = v_phy(i,k,j)+v_frame
      ENDDO
      ENDDO
      ENDDO
   ELSE
      DO j = jts,jtf
      DO k = kts,ktf
      DO i = its,itf
         u_phytmp(i,k,j) = u_phy(i,k,j)
         v_phytmp(i,k,j) = v_phy(i,k,j)
      ENDDO
      ENDDO
      ENDDO
   ENDIF
!
   obsnudge = .false.
   convtran = .false.
   sfcqwcld = .false.
   doregime = .false.
   psimpsih = .false.
!
!  ==================
   IF (dosfclay) THEN
!  ==================
!
!----------------------------------------------------------------------
!  Invoke SFCLAY: surface layer (between surface and pbl)
!----------------------------------------------------------------------
!
   sfclay_select: SELECT CASE(sf_sfclay_physics)

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (SFCLAYSCHEME)
         CALL wrf_debug(100,'in SFCLAY')
         CALL SFCLAY(u_phytmp,v_phytmp,t_phy,qv_curr,               &
                  p_phy,dz8w,PSFC,CHS,CHS2,CQS2,CPM,                &
                  ZNT,UST,PBLH,MAVAIL,ZOL,xtstar,REGIME,PSIM,PSIH,  &
                  XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL,     &
                  U10,V10,TH2,T2M,Q2M,USTM,                         &
                  GZ1OZ0,WSPD,BR,ISFFLX,ISFTCFLX,DX,                &
                  xfm,xfh,xfq,xf10,xfq10,                           &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCSTRES(U10,V10,T2M,PSFC,UST,xtaux,xtauy,            &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,itf, jts,jtf, kts,ktf  )
         obsnudge = .true.
         convtran = .true.
         sfcqwcld = .true.

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (MYNNSFCSCHEME)
         CALL wrf_debug(100,'in MYNNSFC')
         CALL MYNNSFC(u_phytmp,v_phytmp,t_phy,qv_curr,              &
                  p_phy,dz8w,PSFC,CHS,CHS2,CQS2,CPM,                &
                  ZNT,UST,PBLH,MAVAIL,ZOL,xtstar,REGIME,PSIM,PSIH,  &
                  XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL,     &
                  U10,V10,TH2,T2M,Q2M,USTM,                         &
                  GZ1OZ0,WSPD,BR,ISFFLX,ISFTCFLX,DX,                &
                  xfm,xfh,xfq,xf10,xfq10,                           &
                  itimestep,th_phy,pi_phy,qc_curr,                  &
                  tsq,qsq,cov,cqc,qcg,                              &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCSTRES(U10,V10,T2M,PSFC,UST,xtaux,xtauy,            &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,itf, jts,jtf, kts,ktf  )
         obsnudge = .true.
         convtran = .true.

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (PXSFCSCHEME)
         CALL wrf_debug(100,'in PXSFCLAY')
         CALL PXSFCLAY(u_phytmp,v_phytmp,t_phy,th_phy,qv_curr,      &
                  p_phy,dz8w,PSFC,CHS,CHS2,CQS2,CPM,                &
                  ZNT,UST,PBLH,MAVAIL,ZOL,xtstar,REGIME,PSIM,PSIH,  &
                  XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL,     &
                  U10,V10,                                          &
                  GZ1OZ0,WSPD,BR,ISFFLX,DX,                         &
                  TH2,T2M,Q2M,USTM,                                 &
                  xfm,xfh,xfq,xf10,xfq10,                           &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCSTRES(U10,V10,T2M,PSFC,UST,xtaux,xtauy,            &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,itf, jts,jtf, kts,ktf  )
         obsnudge = .true.
         convtran = .true.
         sfcqwcld = .true.

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (MYJSFCSCHEME)
         CALL wrf_debug(100,'in MYJSFC')
         CALL MYJSFC(itimestep,HT,dz8w,p_phy,p8w,th_phy,t_phy,      &
                  qv_curr,qc_curr,u_phytmp,v_phytmp,                &
                  TSK,QSFC,THZ0,QZ0,UZ0,VZ0,                        &
                  LOWLYR,XLAND,UST,ZNT,Z0,PBLH,                     &
                  MAVAIL,RMOL,AKHS,AKMS,BR,                         &
                  CHS,CHS2,CQS2,HFX,QFX,                            &
                  LH,FLHC,FLQC,QGH,CPM,CT,                          &
                  U10,V10,T2M,TH2,Q2M,USTM,xtstar,ZOL,              &
                  xfm,xfh,xfq,xf10,xfq10,CHKLOWQ,                   &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCSTRES(U10,V10,T2M,PSFC,UST,xtaux,xtauy,            &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,itf, jts,jtf, kts,ktf  )
         DO j = jts,jtf
         DO i = its,itf
            WSPD(i,j) = max(0.1,sqrt(u_phytmp(i,kts,j)**2+v_phytmp(i,kts,j)**2))
         ENDDO
         ENDDO
         obsnudge = .true.
         convtran = .true.
         sfcqwcld = .true.
         doregime = .true.
         psimpsih = .true.

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (QNSESFCSCHEME)
         CALL wrf_debug(100,'in QNSESFC')
         CALL QNSESFC(itimestep,HT,dz8w,p_phy,p8w,th_phy,t_phy,     &
                  qv_curr,qc_curr,u_phytmp,v_phytmp,                &
                  TSK,QSFC,THZ0,QZ0,UZ0,VZ0,                        &
                  LOWLYR,XLAND,UST,ZNT,Z0,PBLH,                     &
                  MAVAIL,RMOL,AKHS,AKMS,BR,                         &
                  CHS,CHS2,CQS2,HFX,QFX,                            &
                  LH,FLHC,FLQC,QGH,CPM,CT,                          &
                  U10,V10,T2M,TH2,Q2M,USTM,xtstar,ZOL,              &
                  xfm,xfh,xfq,xf10,xfq10,CHKLOWQ,                   &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCSTRES(U10,V10,T2M,PSFC,UST,xtaux,xtauy,            &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,itf, jts,jtf, kts,ktf  )
         DO j = jts,jtf
         DO i = its,itf
            WSPD(i,j) = max(0.1,sqrt(u_phytmp(i,kts,j)**2+v_phytmp(i,kts,j)**2))
         ENDDO
         ENDDO
         obsnudge = .true.
         convtran = .true.
         sfcqwcld = .true.
         doregime = .true.
         psimpsih = .true.

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (GFSSFCSCHEME)
         CALL wrf_debug(100, 'in GFSSFC')
         CALL SF_GFS(u_phytmp,v_phytmp,t_phy,qv_curr,p_phy,         &
                  PSFC,CHS,CHS2,CQS2,CPM,                           &
                  ZNT,UST,PSIM,PSIH,                                &
                  XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,                   &
                  QGH,QSFC,U10,V10,TH2,T2M,Q2M,                     &
                  RMOL,ZOL,xtstar,USTM,                             &
                  xfm,xfh,xfq,xf10,xfq10,                           &
                  GZ1OZ0,WSPD,BR,ISFFLX,                            &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCSTRES(U10,V10,T2M,PSFC,UST,xtaux,xtauy,            &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,itf, jts,jtf, kts,ktf  )
         obsnudge = .true.
         convtran = .true.
         sfcqwcld = .true.
         doregime = .true.

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (GFDLSFCSCHEME)
         IF (sf_land_physics == GFDLSLAB) THEN
            ntsflg = 1
            CALL wrf_debug(100, 'in GFDLSLAB')
         ELSE
            ntsflg = 0
            CALL wrf_debug(100, 'in GFDLSFC')
         ENDIF
         CALL SF_GFDL(u_phytmp,v_phytmp,t_phy,qv_curr,p_phy,        &
                  MAVAIL,                                           &
                  xfm,xfh,xfq,xf10,xfq10,xtstar,ZOL,RMOL,           &
                  xtaux,xtauy,GSW,GLWd,                             &
                  PSFC,CHS,CHS2,CQS2,CPM,dtimesf,                   &
                  ZNT,UST,PSIM,PSIH,                                &
                  XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,                   &
                  QGH,QSFC,U10,V10,                                 &
                  GZ1OZ0,WSPD,BR,ISFFLX,ntsflg,SFENTH,              &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         DO j = jts,jtf
         DO i = its,itf
            USTM(i,j) = UST(i,j)
         ENDDO
         ENDDO
         obsnudge = .true.
         convtran = .true.
         sfcqwcld = .true.
         doregime = .true.

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE DEFAULT
         CALL wrf_message('xxxx The surface layer scheme does not exist')

   END SELECT sfclay_select
!
      IF (convtran) THEN
         DO j = jts,jtf
         DO i = its,itf
                 ustk = UST(i,j)*karman
            xram(i,j) = xfm(i,j)/ustk
            xrah(i,j) = xfh(i,j)/ustk
            xraq(i,j) = xfq(i,j)/ustk
            xfh2(i,j) = ustk/CHS2(i,j)
            xfq2(i,j) = ustk/CQS2(i,j)
         ENDDO
         ENDDO
      ENDIF
!
      IF (psimpsih) THEN
         DO j = jts,jtf
         DO i = its,itf
          GZ1OZ0(i,j) = ALOG(ZASF(i,j)/ZNT(i,j))
            PSIM(i,j) = GZ1OZ0(i,j)-xfm(i,j)
            PSIH(i,j) = GZ1OZ0(i,j)-xfh(i,j)
         ENDDO
         ENDDO
      ENDIF

!  =====
   ENDIF
!  =====
!
!  ==================
   IF (dosfcflx) THEN
!  ==================
!
   xlfrdt = xlf*rhowater/dtimesf       ! conversion for SNOMELT=>SNOPCX
!
!----------------------------------------------------------------------
!  Invoke land surface model
!----------------------------------------------------------------------
!
   land_select: SELECT CASE(sf_land_physics)

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (CSSPSCHEME)
         CALL wrf_debug(100,'in CSSP')
#ifdef CalUCM
!--CSSP uses as (IO) no WRF-interfacing urban variable
#endif
!Yufeng : a way to initilize biolai!? :D
         DO j = jts,jtf
         DO i = its,itf
         if(MLU_PER(i,j)>0.01) then
         biolai(i,j)=max(biolai(i,j),0.0007 * 1.7)
         endif
         ENDDO
         ENDDO
         CALL CSSP (biocro_switch, scheme_soilwater ,                  &
                  cdate   ,dtimesf ,nl_soil ,nl_snow ,sw_band ,fdtcza  ,&
                  p_phy   ,t_phy   ,qv_curr ,u_phy   ,v_phy   ,z       ,&
                  p8w     ,xrainc  ,xrainnc ,xsnowc  ,xsnownc ,PBLH    ,&
                  GLWd    ,xswi    ,HT      ,XLONG   ,XLAT    ,xfveg   ,&
                  SC_WATER,SC_LANDU,BEDROCK ,                           &
                  xbfdir  ,xwtslope,xuptkw  ,                           &
                  xfrost_table,xzwt_perched , xqflx_sub_snow ,          &
#ifdef MHydro
                  xfdshc  ,                                             &
#endif
#ifdef GenSWC
                  xgswcA  ,xgswcB  ,xgswcC  ,xgswcD  ,xvwc10  ,xvwchs  ,&
#endif
                  xbsw    ,                                             &
                  xcsol   ,xporsl  ,xphi0   ,xdkmg   ,xdksatu ,         &
                  xdkdry  ,xhksati ,xhksm   ,xrootfr ,xwdry   ,         &
#ifdef DoVAST
                  xswdepth,xrtws   ,                                    &
#endif
#ifdef OroRad
                  xfcor   ,xfsky   ,rad_oro ,GSW     ,                  &
#endif
#ifdef LXZALB
                  xsalf   ,                                             &
#else
                  xalbsol ,                                             &
#endif

#ifdef NewLAKE
                  xfrlak, xdplak, xcoszn, xzlake, xdzlak,               &
                  xstke1, xdpice, xt_ice, xlkice, xlktmp,               &
                  lakemodeltype,&
#endif


                  xtss    ,xwliq   ,xwice   ,xz      ,xdz     ,         &
                  TSK     ,xtlsun  ,xtlsha  ,CANWAT  ,                  &!-lxz xtg=>TSK
                  SNOWH   ,SNOW    ,xsag    ,SNOWF   ,xsigf   ,         &
                  xlai    ,xsai    ,xgreen  ,schveg  ,getveg  ,setalb  ,&
                  xalb    ,xssun   ,xssha   ,biolai  ,MLU_PER ,         &
                  xfsun   ,xthermk ,xextkb  ,xextkd  ,                  &
                  EMISS   ,coszan  ,eccfn   ,xalbg   ,xalbv   ,         &
                  xtg     ,T2M     ,Q2M     ,xsabg   ,xsabvg  ,xolrg   ,&!-lxz xtavg=>xtg for rad with emiss
                  xtaux   ,xtauy   ,HFX     ,QFX     ,LH      ,         &
                  xfsenl  ,xfevpl  ,xetr    ,xfseng  ,xfevpg  ,GRDFLX  ,&
                  xshfdt  ,SNOMLT  ,xzwt    ,xtsw    ,xwt     ,xqrchrg ,&
                  xrsur   ,xrbas   ,xrdrn   ,xrsat   ,xrnof   ,xfcov   ,&
                  UST     ,QST     ,xtstar  ,ZNT     ,ZOL     ,BR      ,&
                  xram    ,xrah    ,xraq    ,xfm     ,xfh     ,xfq     ,&
                  xfh2    ,xfq2    ,xf10    ,xfq10   ,U10     ,V10     ,&
                  xparsun ,xparsha ,xsabvsun,xsabvsha,WSPD    ,         &
                  xrst    ,xassim  ,xrespc  ,xxerr   ,xzerr   ,         &
                  VEGFRA  ,                                             &
#ifdef CalUCM
                  ucmdoit ,                  &
#endif
                  ! control parameter for different xemis for different rad
                  ! should be changed later
                  ra_lw_physics,                                        &
                  charn,                                                &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         obsnudge = .true.
         sfcqwcld = .true.
         doregime = .true.
         psimpsih = .true.

         !---note that xsabvg is updated by CSSP due to changing albedo
         !   while GSW accounts for both changing albedo & orographic effects

         k = kts                   ! at the lowest atmospheric layer
         DO j = jts,jtf
         DO i = its,itf
            PRECSFC(i,j) = xrainnc(i,j) + xsnownc(i,j) &
                         +  xrainc(i,j) +  xsnowc(i,j)
            FPRFRZN(i,j) =(xsnownc(i,j) +  xsnowc(i,j)) &
                         /(PRECSFC(i,j) + 1.E-8)

            if (GSWd(i,j) >= 0.01) then
                 GSWu(i,j) =  GSWd(i,j) -  xsabvg(i,j)
             if (GSWu(i,j)>= 0.01) &    ! otherwise no change but using input of previous state
               ALBEDO(i,j) = min(.999,GSWu(i,j)/GSWd(i,j))
            endif

             GZ1OZ0(i,j) = ALOG(ZASF(i,j)/ZNT(i,j))
               PSIM(i,j) = GZ1OZ0(i,j)-xfm(i,j)
               PSIH(i,j) = GZ1OZ0(i,j)-xfh(i,j)
                QGH(i,j) = SATQM(t_phy(i,k,j),p_phy(i,k,j))
               RMOL(i,j) = ZOL(i,j)/ZASF(i,j)
                TH2(i,j) = T2M(i,j)*(1.E5/PSFC(i,j))**ROVCP
               USTM(i,j) = UST(i,j)
                    ustk = UST(i,j)*karman
                CHS(i,j) = ustk/xfh(i,j)      ! = 1.0/xrah(i,j)
               CHS2(i,j) = ustk/xfh2(i,j)
               CQS2(i,j) = ustk/xfq2(i,j)
            sclandu(i,j) = nint(SC_LANDU(i,j))

            if (HFX(i,j) > 0.01) then
               QSFC(i,j) = Q2M(i,j) + CP*(TSK(i,j)-T2M(i,j))* &
                                        (CHS2(i,j)*QFX(i,j))/(CQS2(i,j)*HFX(i,j))
               QSFC(i,j) = max(0.,min(QSFC(i,j),SATQM(TSK(i,j),PSFC(i,j))))
            else
               QSFC(i,j) = Q2M(i,j)
            endif

         ENDDO
         ENDDO
         k = kts                   ! at the lowest atmospheric layer
         CALL SFCETP(TSK,CHS,PSFC,                                  &
                  p_phy(:,k,:),t_phy(:,k,:),qv_curr(:,k,:),QGH,     &
                  PRECSFC,SNOWF,SNOW,FPRFRZN,                       &
                  sclandu,XICE,ic_ice,fc_ice,                       &
                  EMISS,GLWd,GSW,GRDFLX,POTEVP,dtimesf,QFX,         &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf )
#ifdef CalUCM
         IF (dourban) THEN
         CALL UCM(dtimesf,cdate,XLONG,XLAT,                         &
                  u_phy,v_phy,t_phy,qv_curr,                        &
                  th_phy,rho,p_phy,p8w,dz8w,                        &
                  PRECSFC,ZNT,CHS,CHS2,EMISS,ALBEDO,                & 
                  TSK,QSFC,HFX,QFX,LH,GRDFLX,GSWd,GLWd,             &
                  TH2,T2M,Q2M,U10,V10,UST,                          &
                  PSIM,PSIH,GZ1OZ0,AKHS,AKMS,                       &
                  ucmdoit,                                          &
                  frc_urb2d,utype_urb2d,                            &
                  num_roof_layers,num_wall_layers,num_road_layers,  &
                  dzr,dzb,dzg,                                      &
                  tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,     &
                  trl_urb3d,tbl_urb3d,tgl_urb3d,                    &
                  xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,      &
                  num_urban_layers,                                 &
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,          &
                  sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,        &
                  a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,          &
                  b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,          &
                  dlg_bep,dl_u_bep,sf_bep,vl_bep,                   &
                  ts_urb2d,t2_urb2d,q2_urb2d,u10_urb2d,v10_urb2d,   &
                  uc_urb2d,ust_urb2d,sh_urb2d,lh_urb2d,g_urb2d,     &
                  rn_urb2d,qs_urb2d,alb_urb2d,ems_urb2d,            &
                  psim_urb2d,psih_urb2d,gz1oz0_urb2d,akms_urb2d,    &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf )
         DO j = jts,jtf
         DO i = its,itf
            IF (utype_urb2d(i,j) > 0 .AND. frc_urb2d(i,j) < 0.99) THEN
               xtg(i,j) = (ts_urb2d(i,j)-TSK(i,j))*frc_urb2d(i,j)/(1.-frc_urb2d(i,j)) &
                              + xtg(i,j)
            ENDIF
         ENDDO
         ENDDO
         ENDIF
#endif
         l = 1                          ! at the top soil layer
         n = l - nl_snow                ! at the top soil layer below snow
         DO j = jts,jtf
         DO i = its,itf
            if (nint(SC_WATER(i,j))>3) then   ! water or ice
             PORSL(i,j) = 1.
            MAVAIL(i,j) = 1.
            else
             PORSL(i,j) = xporsl(i,l,j)       ! skin soil
            MAVAIL(i,j) = (xwliq(i,n,j)+xwice(i,n,j))/(xporsl(i,l,j)*xdz(i,n,j))
            endif

            SFCEVP(i,j) = SFCEVP(i,j) + QFX(i,j)*dtimesf

            TSLB(i,l,j) = xtss(i,n,j)

            if (nint(SC_WATER(i,j))>6) cycle     ! below=initial over ocean
            SFRUNOFF(i,j) = SFRUNOFF(i,j) + xrsur(i,j)*dtimesf
            UDRUNOFF(i,j) = UDRUNOFF(i,j) + (xrnof(i,j)-xrsur(i,j))*dtimesf
            ACSNOW(i,j) = ACSNOW(i,j) + (xsnowc(i,j)+xsnownc(i,j))*dtimesf
            ACSNOM(i,j) = ACSNOM(i,j) + SNOMLT(i,j)
            SNOPCX(i,j) = SNOPCX(i,j) - SNOMLT(i,j)*xlfrdt
         ENDDO
         ENDDO

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (NMMLSMSCHEME)
         CALL wrf_debug(100,'in NMM LSM')
         DO j = jts,jtf
         DO i = its,itf
            PRECSFC(i,j) = xrainnc(i,j) + xsnownc(i,j) &
                         +  xrainc(i,j) +  xsnowc(i,j)
            FPRFRZN(i,j) =(xsnownc(i,j) +  xsnowc(i,j)) &
                         /(PRECSFC(i,j) + 1.E-8)
         ENDDO
         ENDDO
#ifdef CalUCM
!--NMM uses as (IO) the WRF-interfacing urban variables
!  which need to be adjusted to remove the urban effect
!  for the LSM solution of the residual(non-urban) area
         IF (dourban) THEN
         DO j = jts,jtf
         DO i = its,itf
            IF (utype_urb2d(i,j) > 0 .AND. frc_urb2d(i,j) < 0.99) THEN
                  TSK(i,j) = (   TSK(i,j)- ts_urb2d(i,j)*frc_urb2d(i,j))/(1.-frc_urb2d(i,j))
               ALBEDO(i,j) = (ALBEDO(i,j)-alb_urb2d(i,j)*frc_urb2d(i,j))/(1.-frc_urb2d(i,j))
            ENDIF
         ENDDO
         ENDDO
         ENDIF
#endif
         CALL NMMLSM(dz8w,qv_curr,p8w,rho,                          &
                  t_phy,th_phy,TSK,CHS,                             &
                  HFX,QFX,QGH,GSW,GLWd,LH,                          &
                  SMSTAV,SMSTOT,SFRUNOFF,                           &
                  UDRUNOFF,IVGTYP,ISLTYP,VEGFRA,POTEVP,             &
                  GRDFLX,SFCEXC,ACSNOW,ACSNOM,SNOPCX,               &
                  ALBBCK,TMN,XLAND,XICE,                            &
                  SNOWF,QSFC,SOILTB,PRECSFC,xrst,PORSL,             &
                  num_soil_layers,dtimesf,DZS,                      &
                  SMOIS,TSLB,SNOW,CANWAT,CPM,FPRFRZN,               &
                  ALBEDO,SNOALB,SH2O,SNOWH,                         &
#ifdef CalUCM
                  ucmdoit ,                  &
#endif
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2M,TH2,Q2M,PSFC, &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
#ifdef CalUCM
         IF (dourban) THEN
         CALL UCM(dtimesf,cdate,XLONG,XLAT,                         &
                  u_phy,v_phy,t_phy,qv_curr,                        &
                  th_phy,rho,p_phy,p8w,dz8w,                        &
                  PRECSFC,ZNT,CHS,CHS2,EMISS,ALBEDO,                & 
                  TSK,QSFC,HFX,QFX,LH,GRDFLX,GSWd,GLWd,             &
                  TH2,T2M,Q2M,U10,V10,UST,                          &
                  PSIM,PSIH,GZ1OZ0,AKHS,AKMS,                       &
                  ucmdoit,                                          &
                  frc_urb2d,utype_urb2d,                            &
                  num_roof_layers,num_wall_layers,num_road_layers,  &
                  dzr,dzb,dzg,                                      &
                  tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,     &
                  trl_urb3d,tbl_urb3d,tgl_urb3d,                    &
                  xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,      &
                  num_urban_layers,                                 &
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,          &
                  sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,        &
                  a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,          &
                  b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,          &
                  dlg_bep,dl_u_bep,sf_bep,vl_bep,                   &
                  ts_urb2d,t2_urb2d,q2_urb2d,u10_urb2d,v10_urb2d,   &
                  uc_urb2d,ust_urb2d,sh_urb2d,lh_urb2d,g_urb2d,     &
                  rn_urb2d,qs_urb2d,alb_urb2d,ems_urb2d,            &
                  psim_urb2d,psih_urb2d,gz1oz0_urb2d,akms_urb2d,    &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         ENDIF
#endif
         l = 1                          ! at the top soil layer
         DO j = jts,jtf
         DO i = its,itf
            if (XLAND(i,j)>1.5) then    ! water
            MAVAIL(i,j) = 1.
            else                        ! land
            MAVAIL(i,j) = SMOIS(i,l,j)/PORSL(i,j)
            xtg(i,j)    = TSK(i,j)
            endif

            SFCEVP(i,j) = SFCEVP(i,j) + QFX(i,j)*dtimesf
            xolrg(i,j)  = STBOLT*EMISS(i,j)* TSK(i,j)**4
            xalb(i,:,j) = ALBEDO(i,j)
         ENDDO
         ENDDO

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (LSMSCHEME)
         CALL wrf_debug(100,'in NOAH LSM')
         DO j = jts,jtf
         DO i = its,itf
            PRECSFC(i,j) = xrainnc(i,j) + xsnownc(i,j) &
                         +  xrainc(i,j) +  xsnowc(i,j)
            FPRFRZN(i,j) =(xsnownc(i,j) +  xsnowc(i,j)) &
                         /(PRECSFC(i,j) + 1.E-8)
         ENDDO
         ENDDO
#ifdef CalUCM
         !--LSM uses as (IO) the WRF-interfacing urban variables
         !  which need to be adjusted to remove the urban effect
         !  for the LSM solution of the residual(non-urban) area
         IF (dourban) THEN
         DO j = jts,jtf
         DO i = its,itf
            IF (utype_urb2d(i,j) > 0 .AND. frc_urb2d(i,j) < 0.99) THEN
                  TSK(i,j) = (   TSK(i,j)- ts_urb2d(i,j)*frc_urb2d(i,j))/(1.-frc_urb2d(i,j))
                 QSFC(i,j) = (  QSFC(i,j)- qs_urb2d(i,j)*frc_urb2d(i,j))/(1.-frc_urb2d(i,j))
                EMISS(i,j) = ( EMISS(i,j)-ems_urb2d(i,j)*frc_urb2d(i,j))/(1.-frc_urb2d(i,j))
            ENDIF
         ENDDO
         ENDDO
         ENDIF
#endif
         CALL LSM(dz8w,p8w,p_phy,rho,qv_curr,                       &
                  u_phy,v_phy,th_phy,t_phy,                         &
                  TSK,HFX,QFX,LH,GRDFLX,QGH,                        &
                  GSW,GSWd,GLWd,SMSTAV,SMSTOT,                      &
                  SFRUNOFF,UDRUNOFF,IVGTYP,ISLTYP,                  &
                  ISURBAN,ISICE,MMINLU,VEGFRA,                      &
                  ALBEDO,ALBBCK,ZNT,Z0,TMN,XLAND,XICE,              &
                  EMISS,EMBCK,SNOWF,QSFC,PRECSFC,FPRFRZN,           &
                  RDLAI2D,xlai,USEMONALB,                           &
                  num_soil_layers,dtimesf,DZS,xrst,PORSL,           &
                  SMOIS,TSLB,SNOW,CANWAT,                           &
                  CHS,CHS2,CQS2,CPM,BR,UST,                         &
                  SH2O,SNOWH,SNOALB,SHDMIN,SHDMAX,                  &
                  ACSNOW,ACSNOM,SNOPCX,POTEVP,                      &
                  SNOTIME,xice_threshold,NOAHRES,                   &
#ifdef CalUCM
                  ucmdoit ,                  &
#endif
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2M,TH2,Q2M,PSFC, &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
#ifdef CalUCM
         IF (dourban) THEN
         CALL UCM(dtimesf,cdate,XLONG,XLAT,                         &
                  u_phy,v_phy,t_phy,qv_curr,                        &
                  th_phy,rho,p_phy,p8w,dz8w,                        &
                  PRECSFC,ZNT,CHS,CHS2,EMISS,ALBEDO,                & 
                  TSK,QSFC,HFX,QFX,LH,GRDFLX,GSWd,GLWd,             &
                  TH2,T2M,Q2M,U10,V10,UST,                          &
                  PSIM,PSIH,GZ1OZ0,AKHS,AKMS,                       &
                  ucmdoit,                                          &
                  frc_urb2d,utype_urb2d,                            &
                  num_roof_layers,num_wall_layers,num_road_layers,  &
                  dzr,dzb,dzg,                                      &
                  tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,     &
                  trl_urb3d,tbl_urb3d,tgl_urb3d,                    &
                  xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,      &
                  num_urban_layers,                                 &
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,          &
                  sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,        &
                  a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,          &
                  b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,          &
                  dlg_bep,dl_u_bep,sf_bep,vl_bep,                   &
                  ts_urb2d,t2_urb2d,q2_urb2d,u10_urb2d,v10_urb2d,   &
                  uc_urb2d,ust_urb2d,sh_urb2d,lh_urb2d,g_urb2d,     &
                  rn_urb2d,qs_urb2d,alb_urb2d,ems_urb2d,            &
                  psim_urb2d,psih_urb2d,gz1oz0_urb2d,akms_urb2d,    &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         ENDIF
#endif
         l = 1                          ! at the top soil layer
         DO j = jts,jtf
         DO i = its,itf
            if (XLAND(i,j)>1.5) then    ! water
            MAVAIL(i,j) = 1.
            else                        ! land
            MAVAIL(i,j) = SMOIS(i,l,j)/PORSL(i,j)
            xtg(i,j)    = TSK(i,j)
            endif

            SFCEVP(i,j) = SFCEVP(i,j) + QFX(i,j)*dtimesf
            xolrg(i,j)  = STBOLT*EMISS(i,j)* TSK(i,j)**4
            xalb(i,:,j) = ALBEDO(i,j)
            xzerr(i,j)  = xzerr(i,j) + NOAHRES(i,j)
         ENDDO
         ENDDO

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (RUCLSMSCHEME)
         CALL wrf_debug(100,'in RUC LSM')
         DO j = jts,jtf
         DO i = its,itf
            PRECSFC(i,j) = xrainnc(i,j) + xsnownc(i,j) &
                         +  xrainc(i,j) +  xsnowc(i,j)
            FPRFRZN(i,j) =(xsnownc(i,j) +  xsnowc(i,j)) &
                         /(PRECSFC(i,j) + 1.E-8)
         ENDDO
         ENDDO
#ifdef CalUCM
         !--RUC uses as (IO) the WRF-interfacing urban variables
         !  which need to be adjusted to remove the urban effect
         !  for the RUC solution of the residual(non-urban) area
         IF (dourban) THEN
         DO j = jts,jtf
         DO i = its,itf
            IF (utype_urb2d(i,j) > 0 .AND. frc_urb2d(i,j) < 0.99) THEN
                  TSK(i,j) = (   TSK(i,j)- ts_urb2d(i,j)*frc_urb2d(i,j))/(1.-frc_urb2d(i,j))
               ALBEDO(i,j) = (ALBEDO(i,j)-alb_urb2d(i,j)*frc_urb2d(i,j))/(1.-frc_urb2d(i,j))
            ENDIF
         ENDDO
         ENDDO
         ENDIF
#endif
         CALL LSMRUC(dtimesf,itimestep,num_soil_layers,             &
                  ZS,PRECSFC,SNOW,SNOWH,SNOWF,FPRFRZN,              &
                  z,p8w,t_phy,qv_curr,qc_curr,rho,                  &
                  GLWd,GSW,EMISS,PORSL,                             &
                  CHS,FLQC,FLHC,MAVAIL,CANWAT,VEGFRA,               &
                  ALBEDO,ZNT,SNOALB,ALBBCK,                         &
                  QSFC,QSG,QVG,QCG,SOILT1,TSNAV,                    &
                  TMN,IVGTYP,ISLTYP,XLAND,XICE,                     &
                  SMOIS,SH2O,SMSTAV,SMSTOT,                         &
                  TSLB,TSK,HFX,QFX,LH,                              &
                  SFRUNOFF,UDRUNOFF,SFCEXC,                         &
                  GRDFLX,ACSNOW,ACSNOM,                             &
                  SMFR3D,KEEPFR3DFLAG,usemyjsfc,                    &
#ifdef CalUCM
                  ucmdoit ,                  &
#endif
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         CALL SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2M,TH2,Q2M,PSFC, &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         k = kts                   ! at the lowest atmospheric layer
         CALL SFCETP(TSK,CHS,PSFC,                                  &
                  p_phy(:,k,:),t_phy(:,k,:),qv_curr(:,k,:),QGH,     &
                  PRECSFC,SNOWF,SNOW,FPRFRZN,                       &
                  IVGTYP,XICE,ic_ice,fc_ice,                        &
                  EMISS,GLWd,GSW,GRDFLX,POTEVP,dtimesf,QFX,         &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf )
#ifdef CalUCM
         IF (dourban) THEN
         CALL UCM(dtimesf,cdate,XLONG,XLAT,                         &
                  u_phy,v_phy,t_phy,qv_curr,                        &
                  th_phy,rho,p_phy,p8w,dz8w,                        &
                  PRECSFC,ZNT,CHS,CHS2,EMISS,ALBEDO,                & 
                  TSK,QSFC,HFX,QFX,LH,GRDFLX,GSWd,GLWd,             &
                  TH2,T2M,Q2M,U10,V10,UST,                          &
                  PSIM,PSIH,GZ1OZ0,AKHS,AKMS,                       &
                  ucmdoit,                                          &
                  frc_urb2d,utype_urb2d,                            &
                  num_roof_layers,num_wall_layers,num_road_layers,  &
                  dzr,dzb,dzg,                                      &
                  tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,     &
                  trl_urb3d,tbl_urb3d,tgl_urb3d,                    &
                  xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,      &
                  num_urban_layers,                                 &
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,          &
                  sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,        &
                  a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,          &
                  b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,          &
                  dlg_bep,dl_u_bep,sf_bep,vl_bep,                   &
                  ts_urb2d,t2_urb2d,q2_urb2d,u10_urb2d,v10_urb2d,   &
                  uc_urb2d,ust_urb2d,sh_urb2d,lh_urb2d,g_urb2d,     &
                  rn_urb2d,qs_urb2d,alb_urb2d,ems_urb2d,            &
                  psim_urb2d,psih_urb2d,gz1oz0_urb2d,akms_urb2d,    &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         ENDIF
#endif
         l = 1                          ! at the top soil layer
         DO j = jts,jtf
         DO i = its,itf
            if (XLAND(i,j)>1.5) then    ! water
            MAVAIL(i,j) = 1.
            else                        ! land
!-yx        MAVAIL(i,j) = SMOIS(i,l,j)/PORSL(i,j)
            xtg(i,j)    = TSK(i,j)
            endif

            SFCEVP(i,j) = SFCEVP(i,j) + QFX(i,j)*dtimesf
            xolrg(i,j)  = STBOLT*EMISS(i,j)* TSK(i,j)**4
            xalb(i,:,j) = ALBEDO(i,j)
         !- xrst is not available
         ENDDO
         ENDDO

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (PXLSMSCHEME)
         CALL wrf_debug(100,'in P-X LSM')
         DO j = jts,jtf
         DO i = its,itf
            PRECSFC(i,j) = xrainnc(i,j) + xsnownc(i,j) &
                         +  xrainc(i,j) +  xsnowc(i,j)
            SNOWSFC(i,j) = xsnownc(i,j) +  xsnowc(i,j)
            FPRFRZN(i,j) = SNOWSFC(i,j) /(PRECSFC(i,j) + 1.E-8)
             SNOMLT(i,j) =    SNOW(i,j) + SNOWSFC(i,j)
         ENDDO
         ENDDO
         !--- make sure that [frc_urb2d=sum(landusef=1,31-33)]
         CALL PXLSM(u_phy,v_phy,dz8w,qv_curr,t_phy,th_phy,rho,      &     
                  PSFC,GSW,GLWd,PRECSFC,EMISS,                      &
                  itimestep,num_soil_layers,DT,                     & 
                  XLAND,ALBBCK,ALBEDO,SNOALB,                       &
                  SMOIS,TSLB,MAVAIL,T2M,Q2M,                        &
                  ZS,DZS,PSIH,                                      &
                  landusef,soilctop,soilcbot,VEGFRA,vegf_px,        &
                  PORSL,xram,xrst,xlai,nlcat,nscat,                 &
                  HFX,QFX,LH,TSK,ZNT,CANWAT,                        &
                  GRDFLX,xfevpg,SHDMIN,SHDMAX,                      &
                  SNOWF,PBLH,RMOL,UST,CAPG,dtimesf,                 &
                  SNOW,SNOWH,SNOWSFC,                               &
                  pxlsm_smois_init,                                 &
#ifdef PX_NDG
                  pxlsm_soil_nudge,anal_interval,                   &
                  t2_ndg_old,t2_ndg_new,q2_ndg_old,q2_ndg_new,      &
                  sn_ndg_old,sn_ndg_new,t2obs,q2obs,                &
#endif
#ifdef CalUCM
                  ucmdoit ,                  &
#endif
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         k = kts                   ! at the lowest atmospheric layer
         CALL SFCETP(TSK,CHS,PSFC,                                  &
                  p_phy(:,k,:),t_phy(:,k,:),qv_curr(:,k,:),QGH,     &
                  PRECSFC,SNOWF,SNOW,FPRFRZN,                       &
                  IVGTYP,XICE,ic_ice,fc_ice,                        &
                  EMISS,GLWd,GSW,GRDFLX,POTEVP,dtimesf,QFX,         &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf )
         DO j = jts,jtf
         DO i = its,itf
               TH2(i,j) = T2M(i,j)*(1.E5/PSFC(i,j))**ROVCP
            ACSNOW(i,j) = ACSNOW(i,j) + (xsnowc(i,j)+xsnownc(i,j))*dtimesf
            SNOMLT(i,j) = SNOMLT(i,j) -   SNOW(i,j)
            ACSNOM(i,j) = ACSNOM(i,j) + SNOMLT(i,j)
            SNOPCX(i,j) = SNOPCX(i,j) - SNOMLT(i,j)*xlfrdt
         ENDDO
         ENDDO
#ifdef CalUCM
         !--PX uses as (IO) no WRF-interfacing urban variable
         IF (dourban) THEN
         CALL UCM(dtimesf,cdate,XLONG,XLAT,                         &
                  u_phy,v_phy,t_phy,qv_curr,                        &
                  th_phy,rho,p_phy,p8w,dz8w,                        &
                  PRECSFC,ZNT,CHS,CHS2,EMISS,ALBEDO,                & 
                  TSK,QSFC,HFX,QFX,LH,GRDFLX,GSWd,GLWd,             &
                  TH2,T2M,Q2M,U10,V10,UST,                          &
                  PSIM,PSIH,GZ1OZ0,AKHS,AKMS,                       &
                  ucmdoit,                                          &
                  frc_urb2d,utype_urb2d,                            &
                  num_roof_layers,num_wall_layers,num_road_layers,  &
                  dzr,dzb,dzg,                                      &
                  tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,     &
                  trl_urb3d,tbl_urb3d,tgl_urb3d,                    &
                  xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,      &
                  num_urban_layers,                                 &
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,          &
                  sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,        &
                  a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,          &
                  b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,          &
                  dlg_bep,dl_u_bep,sf_bep,vl_bep,                   &
                  ts_urb2d,t2_urb2d,q2_urb2d,u10_urb2d,v10_urb2d,   &
                  uc_urb2d,ust_urb2d,sh_urb2d,lh_urb2d,g_urb2d,     &
                  rn_urb2d,qs_urb2d,alb_urb2d,ems_urb2d,            &
                  psim_urb2d,psih_urb2d,gz1oz0_urb2d,akms_urb2d,    &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )
         ENDIF
#endif
         l = 1                          ! at the top soil layer
         DO j = jts,jtf
         DO i = its,itf
            if (XLAND(i,j)>1.5) then    ! water
            MAVAIL(i,j) = 1.
            else                        ! land
            MAVAIL(i,j) = SMOIS(i,l,j)/PORSL(i,j)
            xtg(i,j)    = TSK(i,j)
            endif

            SFCEVP(i,j) = SFCEVP(i,j) + QFX(i,j)*dtimesf
            xolrg(i,j)  = STBOLT*EMISS(i,j)* TSK(i,j)**4
            xalb(i,:,j) = ALBEDO(i,j)
         ENDDO
         ENDDO

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (SLABSCHEME)
         CALL wrf_debug(100,'in SLAB')
         radiation = ra_lw_physics>0
         CALL SLAB(t_phy,qv_curr,p_phy,FLHC,FLQC,                   &
                 PSFC,XLAND,TMN,HFX,QFX,LH,TSK,QSFC,                &
                 GSW,GLWd,CAPG,THC,SNOWF,EMISS,                     &
                 dtimesf,IFSNOW,                                    &
                 TSLB,ZS,DZS,num_soil_layers,radiation,             &
                 ids,ide ,jds,jde ,kds,kde ,&
                 ims,ime ,jms,jme ,kms,kme ,&
                 its,itf ,jts,jtf ,kts,ktf ,&
                 TML,T0ML,HML,H0ML,HUML,HVML,UST,U_PHY,V_PHY,       &
                 F,OMLCALL,OML_GAMMA)
         DO j = jts,jtf
         DO i = its,itf
            SFCEVP(i,j) = SFCEVP(i,j) + QFX(i,j)*dtimesf
            xolrg(i,j)  = STBOLT*EMISS(i,j)* TSK(i,j)**4
            xalb(i,:,j) = ALBEDO(i,j)
         !  CQS2 accounts for MAVAIL for SFCDIAGS Q2M
            CQS2(i,j)   = CQS2(i,j)*MAVAIL(i,j)
            xtg(i,j)    = TSK(i,j)
         ENDDO
         ENDDO
         CALL SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2M,TH2,Q2M,PSFC, &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (GFDLSLAB)
         DO j = jts,jtf
         DO i = its,itf
            SFCEVP(i,j) = SFCEVP(i,j) + QFX(i,j)*dtimesf
            xolrg(i,j)  = STBOLT*EMISS(i,j)* TSK(i,j)**4
            xalb(i,:,j) = ALBEDO(i,j)
            xtg(i,j)    = TSK(i,j)
         !  CQS2 accounts for MAVAIL for SFCDIAGS Q2M
            CQS2(i,j)   = CQS2(i,j)*MAVAIL(i,j)
         ENDDO
         ENDDO
         CALL SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2M,TH2,Q2M,PSFC, &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE DEFAULT
         CALL wrf_debug(100,'xxxx The land surface model does not exist')

   END SELECT land_select

!  =====
   ENDIF
!  =====
!
      sfcqwcld = sfcqwcld .AND. (bl_pbl_physics==MYNNPBLSCHEME)
      IF (sfcqwcld) THEN
         DO j = jts,jtf
         DO i = its,itf
            cqc(i,j) = CHS(i,j)
         ENDDO
         ENDDO
      ENDIF
!
      IF (doregime) THEN
         DO j = jts,jtf
         DO i = its,itf
            IF     (BR(i,j) >= 0.2) THEN
                REGIME(i,j) =  1.0        ! nighttime stable condition
            ELSEIF (BR(i,j) >  0.0) THEN
                REGIME(i,j) =  2.0        ! damped mechanical turbulent condition
            ELSEIF (BR(i,j) == 0.0) THEN
                REGIME(i,j) =  3.0        ! forced convection condition
            ELSE
                REGIME(i,j) =  4.0        ! free convection condition
            ENDIF
         ENDDO
         ENDDO
      ENDIF
!
!  ===============================================================
   IF (PRESENT(skin_sst) .AND. PRESENT(dtwsk) .AND. skin_sst) THEN
!  ===============================================================
!
!----------------------------------------------------------------------
!  Do skin sea surface temperature prediction
!----------------------------------------------------------------------
!
         CALL SSTSKIN(DELT=dtimesf,DTWSK=dtwsk,TSK=tsk,TMN=tmn,     &
                  SCWATER=scwater,ICWATER=icwater,GLW=glwd,GSW=gsw, &
                  HFX=hfx,LFX=lh,UST=ust,EMISS=emiss,               &
                  IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,  &
                  IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme,  &
                  ITS=its,ITE=itf,JTS=jts,JTE=jtf,KTS=kts,KTE=ktf)
!  =====
   ENDIF
!  =====
!
!  ===================
   IF (dosfcvisc) THEN
!  ===================
!
!----------------------------------------------------------------------
!  Do diagnostic surface viscous sublayer over oceans for MYJ PBL
!----------------------------------------------------------------------
!
         initial = itimestep-1
         CALL SFCVISCOUS(XLAND=xland,FM=xfm,FH=xfh,UST=ust,         &
                  TSK=tsk,QSFC=qsfc,PSFC=psfc,QGH=qgh,              &
                  ULOW=u_phytmp(ims:ime,kts:kts,jms:jme),           &
                  VLOW=v_phytmp(ims:ime,kts:kts,jms:jme),           &
                  QLOW=qv_curr(ims:ime,kts:kts,jms:jme),            &
                  THLOW=th_phy(ims:ime,kts:kts,jms:jme),            &
                  UZ0=uz0,VZ0=vz0,QZ0=qz0,THZ0=thz0,                &
                  AKMS=akms,AKHS=akhs,CHKLOWQ=chklowq,CGT=ct,       &
!-                SETCGT=1,                                         &
                  INITIAL=initial,                                  &
!-                NITERAT=5,                                        &
!-                QCLOW=qc_curr(ims:ime,kts:kts,jms:jme),           &
!-                TLOW=t_phy(ims:ime,kts:kts,jms:jme),              &
!-                DZLOW=dz8w(ims:ime,kts:kts,jms:jme),              &
!-                PBLH=pblh,                                        &
                  IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,  &
                  IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme,  &
                  ITS=its,ITE=itf,JTS=jts,JTE=jtf,KTS=kts,KTE=ktf)
!  =====
   ENDIF
!  =====
!
!  ==================
   IF (doocndyn) THEN
!  ==================
!
!----------------------------------------------------------------------
!  Invoke ocean surface model
!----------------------------------------------------------------------
!
#ifdef XOML
!^^^^^^^^^^^ FOR XOML ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! chao to prevent unstable OI in bluewaters
!        do j = jts, jtf
!           do i = its, itf
!              if(SSTB_1(i,j) <= 150. .and. SSTB_1(i,j) >= 1) then
!                 print*, 'xxx', i, j, SSTB_1(i,j), SC_WATER(i,j)
!              end if
!           end do
!        end do

        !print*, 'yyy', gsw
        !print*, 'zzz', tsk
        !print*,  oml_immpene  , oml_jrvtype , oml_tskin  , couple_oml,    &
        !oml_rlaxopt , oml_kbkspar , oml_rlaxtau , oml_isrpene ,      &
        !cdate, itimestep , sfdt,         &
        !oml_dt, omlzly ,omlzlv ,omldz  , sst_nstep ,        &
        !ids,ide, jds,jde, nl_oml,                      &
        !ims,ime, jms,jme,                              &
        !its,ite, jts,jte,                              &
        !sc_water,xlat,                             &
        !SSTB_1,gsw,glwd-emiss*5.67051e-8*(tsk**4), &
!       !5.67051e-8 ! Stefan-Boltzmann constant  [W/m2/K4]
        !raincv+rainncv, &
        !qfx, &
        !hfx,&
        !lh, &
        !xtaux, &
        !xtauy ,     &
        !t_oml,s_oml,u_oml,v_oml,e_oml ,                   &
        !qcl_oml,buo_oml,gbo_oml,eps_oml,vtz_oml,ksz_oml,            &
        !kqz_oml,rit_oml,eln_oml,                    &
        !hobl   ,  dtw1, sstsk  ,                    &
        !prevdsst, tempsst, obsprevdsst, obstempsst, &
        !nl_soil,nl_snow,  TSK     ,xtg     ,xtss 


        CALL XOML_driver(                                           &
        oml_immpene  , oml_jrvtype , oml_tskin  , couple_oml,    &
        oml_rlaxopt , oml_kbkspar , oml_rlaxtau , oml_isrpene ,      &
        cdate, itimestep , sfdt,         &
        oml_dt, omlzly ,omlzlv ,omldz  , sst_nstep ,        &
        ids,ide, jds,jde, nl_oml,                      &
        ims,ime, jms,jme,                              &
        its,ite, jts,jte,                              &
        sc_water,xlat,                             &
        SSTB_1,gsw,glwd-emiss*5.67051e-8*(tsk**4), &
!       5.67051e-8 ! Stefan-Boltzmann constant  [W/m2/K4]
        raincv+rainncv, qfx, hfx,lh, xtaux, xtauy ,     &
        t_oml,s_oml,u_oml,v_oml,e_oml ,                   &  
        qcl_oml,buo_oml,gbo_oml,eps_oml,vtz_oml,ksz_oml,            &  
        kqz_oml,rit_oml,eln_oml,                    &  
        hobl   ,  dtw1, sstsk  ,                    &  
        prevdsst, tempsst, obsprevdsst, obstempsst, &
        nl_soil,nl_snow,  TSK     ,xtg     ,xtss   )

!$$$$$$$$$$$ FOR XOML $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#endif
   ocean_select: SELECT CASE(sf_ocean_physics)

      CASE (SOMSCHEME)
         CALL wrf_debug(100,'in SOM')
         CALL SOM(cdate   ,dtimesf ,nl_soil ,nl_snow ,              &
                  SC_WATER,XLAND   ,XICE    ,                       &
                  TSK     ,xtg     ,xtss    ,SNOW    ,SNOWH   ,     &
                  GRDFLX  ,xshfdt  ,EMISS   ,                       &
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,itf ,jts,jtf ,kts,ktf  )

      CASE DEFAULT
         CALL wrf_debug(100,'xxxx The ocean surface model does not exist')

   END SELECT ocean_select

!  =====
   ENDIF
!  =====
!
!  ==================
   IF (obsnudge) THEN
!  ==================
!
!----------------------------------------------------------------------
!  Save for observation nudging
!----------------------------------------------------------------------
!
      IF ( PRESENT(uratx) .AND. PRESENT(vratx) .AND. PRESENT(tratx) ) THEN
         k = kts                        ! at the lowest layer
         DO j = jts,jtf
         DO i = its,itf
            tratx(i,j) = th_phy(i,k,j)/TH2(i,j)
            uratx(i,j) = 1.2
            vratx(i,j) = 1.2
            if (abs(U10(i,j)) > 1.e-10) uratx(i,j) = u_phytmp(i,k,j)/U10(i,j)
            if (abs(V10(i,j)) > 1.e-10) vratx(i,j) = v_phytmp(i,k,j)/V10(i,j)
         ENDDO
         ENDDO
      ENDIF

!  =====
   ENDIF
!  =====
!
!----------------------------------------------------------------------
!  Save diagnostic output
!----------------------------------------------------------------------
!
      IF ( PRESENT(ck)    .AND. PRESENT(cd)    .AND. &
           PRESENT(cka)   .AND. PRESENT(cda)   ) THEN
         DO j = jts,jtf
         DO i = its,itf
            kofm = karman/xfm(i,j)
            kof10 = karman/xf10(i,j)
             ck(i,j) = kof10*(karman/xfq10(i,j))
             cd(i,j) = kof10*(kof10)
            cka(i,j) = kofm*(karman/xfq(i,j))
            cda(i,j) = kofm*(kofm)
         ENDDO
         ENDDO
      ENDIF
!
      IF ( PRESENT(t2min) .AND. PRESENT(t2max) ) THEN
!-       IF ( (itimestep == 1) .OR. ISNEWDAY(cdate,dtimesf) ) THEN
         IF ( (itimestep == 1) ) THEN   ! do re-initialization in {med_hist_out(stream=0)}
            DO j = jts,jtf
            DO i = its,itf
               t2min(i,j) = T2M(i,j)
               t2max(i,j) = T2M(i,j)
            ENDDO
            ENDDO
         ELSE
            DO j = jts,jtf
            DO i = its,itf
               t2min(i,j) = min(t2min(i,j),T2M(i,j))
               t2max(i,j) = max(t2max(i,j),T2M(i,j))
            ENDDO
            ENDDO
         ENDIF
      ENDIF
!
      IF ( PRESENT(aswd)  .AND. PRESENT(aswu)  .AND. &
           PRESENT(alwd)  .AND. PRESENT(alwu)  .AND. &
           PRESENT(ahfx)  .AND. PRESENT(alfx)  .AND. &
           PRESENT(agfx)  ) THEN
!-       IF ( (itimestep == 1) .OR. ISNEWDAY(cdate,dtimesf) ) THEN
         IF ( (itimestep == 1) ) THEN   ! do re-initialization in {med_history_out}
            DO j = jts,jtf
            DO i = its,itf
               aswd(i,j) =   GSWd(i,j)
               aswu(i,j) =   GSWu(i,j)
               alwd(i,j) =   GLWd(i,j)
               alwu(i,j) =  xolrg(i,j)
               ahfx(i,j) =    HFX(i,j)
               alfx(i,j) =     LH(i,j)
               agfx(i,j) = GRDFLX(i,j)
            ENDDO
            ENDDO
         ELSE
            DO j = jts,jtf
            DO i = its,itf
               aswd(i,j) =   GSWd(i,j) + aswd(i,j)
               aswu(i,j) =   GSWu(i,j) + aswu(i,j)
               alwd(i,j) =   GLWd(i,j) + alwd(i,j)
               alwu(i,j) =  xolrg(i,j) + alwu(i,j)
               ahfx(i,j) =    HFX(i,j) + ahfx(i,j)
               alfx(i,j) =     LH(i,j) + alfx(i,j)
               agfx(i,j) = GRDFLX(i,j) + agfx(i,j)
            ENDDO
            ENDDO
         ENDIF
      ENDIF
!--------------Yufengs------------------------
!variables needed for BioCro
!To calc RH, we need, PSFC, Q2M and tair.
!rh=0.01*0.263*psfc*q2m*exp(17.67*(tair0-tair)/(tair-29.65))
!rh[rh>1]=1
!xswi:solar radiation; T2M:temperature; :RH; WSAVG:windspeed; PRAVG:precipitation
  ! IF ( wrf_dm_on_monitor() ) THEN
   !   WRITE(*,*)'start & end memory at i index are: ',ims, ime
   !   WRITE(*,*)'start & end domain at i index are: ',ids, ide
   !   WRITE(*,*)'start & end tile at i index are: ',its, itf
   !   WRITE(*,*)'# of tiles is: ',num_tiles
   !ENDIF

!ims,ime
!IF ( wrf_dm_on_monitor() ) THEN
        ! get the year, julian day and gmt time
!         CALL get_julgmt(cdate, curyr, juldy, gmtim)

         ! aussume the soil profile is same

 !        print*, 'dates are...',  juldy,gmtim,gmtim==INT(gmtim)
!print *, 'YYYYYYY size of xwliq is', size(xwliq,2)
!ENDIF
if (biocro_switch == 1) then
CALL get_julgmt(cdate, curyr, juldy, gmtim)
doy=juldy
hour=gmtim
print *, 'dates are...',  cdate,curyr,juldy,gmtim,abs(hour-NINT(hour)) < 0.05
!print *, 'bio_onoff is',  bio_onoff
!print *, 'itimestep,sfdt are',itimestep,sfdt
!print *,"hour is",hour
if ((abs(hour-NINT(hour)) < 0.05)) then
!print *,"integer time hour is",hour
!soilLayers=11
do j=jts,jtf
do i=its,itf
if ((MLU_PER(i,j) > 0.01) .and. (nint(SC_WATER(i,j))==2 .or. nint(SC_WATER(i,j))==3)) then
solar_dir  = xswi(i,1,j)*2.1!(xswi(i,1,j) + xswi(i,2,j)) * 2.1 ! W/m2 to umol PAR m-2 s-1
solar_diff = xswi(i,3,j)*2.1!(xswi(i,3,j) + xswi(i,4,j)) * 2.1 ! W/m2 to umol PAR m-2 s-1
temp=T2M(i,j)
rh=0.01*0.263*PSFC(i,j)*Q2M(i,j)*exp(17.67*(273.15-temp)/(temp-29.65))
rh=min(1.0,rh)
temp=temp-273.15  ! K to C degree
windspeed = WSPD(i,j)
precip    = PRECSFC(i,j)*3600 ! mm/s to mm/hour
do s=1,nl_soil
sand=fr_sand(i,s,j)*100.
clay=fr_clay(i,s,j)*100.
om_frac=0.025
om_tkm=0.25
om_watsat=0.9
om_b=2.7
watsat = 0.489 - 0.00126*sand
bsw    = 2.91 + 0.159*clay
sucsat = 10.0 * ( 10.0**(1.88-0.0131*sand) )
!bd     = (1.0 - watsat)*2.7e3
watsat = (1. - om_frac)*watsat + om_watsat*om_frac   
tkm    = (1.-om_frac)*(8.80*sand+2.92*clay)/(sand+clay)+om_tkm*om_frac                               
bsw    = (1.-om_frac)*bsw + om_frac*om_b
wiltp(s)   = (2.57e5/sucsat)**(-1.0/bsw)*watsat
fieldc(s) = watsat * (0.1/(xhksati(i,s,j)*86400.0))**(1.0/(2.0*bsw+3.0))
vol_ice(s)  = min(xporsl(i,s,j), xwice(i,s+5,j)/(xdz(i,s+5,j)*denice))
eff_porosity(s) = xporsl(i,s,j)-vol_ice(s)
!iWatCont(s)=min(eff_porosity(s),xwliq(i,s+5,j)/(denliq*xdz(i,s+5,j)))! This should be soil water volumetric content (%)
!iWatCont(s)=max(iWatCont(s),0.3)
!print *,"vol_ice,xwice,xporsl,eff_porosity,iWatCont,xwliq,s,i,j are",&
!vol_ice(s),xwice(i,s+5,j),xporsl(i,s,j),eff_porosity(s),iWatCont(s),xwliq(i,s+5,j),s,i,j
cws(i,s,j)=min(eff_porosity(s),(xwliq(i,s+5,j)*1.e-3)/xdz(i,s+5,j))
!fieldc(s) = xporsl(i,s,j) * (0.1/(xhksati(i,s,j)*86400.0))**(1.0/(2.0*xbsw(i,s,j)+3.0))
!wiltp(s)  = (1.5e5/xphi0(i,s,j))**(-1.0/xbsw(i,s,j)) * (xporsl(i,s,j)-vol_ice(s))
enddo
doy0 = floor(doy0_coefs(1)*XLAT(i,j)+ doy0_coefs(2)*XLONG(i,j)+doy0_coefs(3)*XLAT(i,j)*XLAT(i,j)+&
doy0_coefs(4)*XLONG(i,j)*XLONG(i,j) + doy0_coefs(5)*XLAT(i,j)*XLONG(i,j))
doy1 = floor(doy1_coefs(1)*XLAT(i,j)+ doy1_coefs(2)*XLONG(i,j)+doy1_coefs(3)*XLAT(i,j)*XLAT(i,j)+&
doy1_coefs(4)*XLONG(i,j)*XLONG(i,j) + doy1_coefs(5)*XLAT(i,j)*XLONG(i,j))
!print *,"i,j,doy0,doy1 are",i,j,doy0,doy1
!call srand(86456)
!print *,"hour,coszan,xcoszn,i,j are",hour,coszan(i,j),xcoszn(i,j),i,j
if (doy >= doy0 .and. doy <= doy1) then
bio_out=GroX(XLAT(i,j),doy,xcoszn(i,j),solar_dir,solar_diff,temp,rh,windspeed,precip,&
    kd,chil,leafwidth,et_equation,heightf,nlayers,initial_biomass,sencoefs,&
    timestep,iSp,SpD,dbpcoefs,thermalp,tbase,vmax1,alpha1,kparm,theta,&
    beta,Rd,Catm,b0,b1,soilcoefs,ileafn,kLN,vmaxb1,alphab1,mresp,wsFun,&
    ws,soilLayers,soilDepths,cws(i,:,j),&
    secs,kpLN,lnb0,lnb1,lnfun,upperT,lowerT,nitroP,stomataWS(i,j),&
    TTc(i,j),Grain(i,j),Leaf(i,j),Stem(i,j),Root(i,j),Rhizome(i,j),&
    LeafN(i,j),Sp(i,j),vmax(i,j),alpha(i,j),&
    doy0,doy1,iWatCont,hour,fieldc,wiltp,BEDROCK(i,j))
!print *,"GroX has ran and TTc,xcoszn(i,j) are", TTc(i,j),xcoszn(i,j)
TTc(i,j)=bio_out%TTc
Grain(i,j)=bio_out%Grain
Leaf(i,j)=bio_out%Leaf
Stem(i,j)=bio_out%Stem
Root(i,j)=bio_out%Root
Rhizome(i,j)=bio_out%Rhizome
LeafN(i,j)=bio_out%LeafN
Sp(i,j)=bio_out%Sp
vmax(i,j)=bio_out%vmax
alpha(i,j)=bio_out%alpha
stomataWS(i,j)=bio_out%stomataWS
cws(i,:,j)=bio_out%cws
biolai(i,j)=bio_out%biolai            !  working!
bioev(i,j)=bio_out%bioev
biotr(i,j)=bio_out%biotr
!xwliq(i,6:16,j)=cws(i,:,j)*1000.0*xdz(i,6:16,j)*MLU_PER(i,j)+&
!                (1.0-MLU_PER(i,j))*xwliq(i,6:16,j)
xfevpl(i,j) = biotr(i,j)/(3.6e4) * MLU_PER(i,j)+(1.0-MLU_PER(i,j))*xfevpl(i,j)
xfevpg(i,j) = bioev(i,j)/(3.6e4) * MLU_PER(i,j)+(1.0-MLU_PER(i,j))*xfevpg(i,j)
else
TTc(i,j)=.0
Grain(i,j)=.0
Leaf(i,j)= 7e-4
Stem(i,j)=0.007
Root(i,j)=0.007
Rhizome(i,j)=7.0
LeafN(i,j)=2.0
Sp(i,j)=.0
vmax(i,j)=.0
alpha(i,j)=.0
stomataWS(i,j)=1.0
!cws(i,:,j)=.0
biolai(i,j)=0.00119
!bioev(i,j)=.0
!biotr(i,j)=.0
endif
end if !if at MLU>0.01

enddo
enddo
end if !if at integer hour 
endif !if biocro is on
   ENDDO  ! end num_tiles Yufeng
   !$OMP END PARALLEL DO
!
!----------------------------------------------------------------------
!  Set control for the use at next step
!----------------------------------------------------------------------
! 
   IF (do_sfpr_avg) THEN
      IF (run_param) THEN
         dt_sfpr_avg = 0.
      ELSE
         dt_sfpr_avg = dt_sfpr_avg + dt
      ENDIF
   ENDIF

!-------------------------------------------------------
   END SUBROUTINE surface_driver


#ifdef DoVAST
!----------------------------------------------------------------------
   SUBROUTINE surface_finish(                                       &
              itimestep,sf_land_physics,                            &
              no_flowdir,                                           &
              nl_soil,nl_snow,STEPSF,DT,DX,DY,                      &
              SC_WATER,xdz,xhksm,xbsw,xphi0,xporsl,xwdry,           &
              xwliq,xwice,xrsat,xrnof,xxerr,                        &
              ZSFCM,xslpxs,xslpys,                                  &
              xdirid,xdirwd,xchleng,xbedslp,                        &
              xrtws,xswdepth,xoutflow,xtsw,                         &
              sfdt    ,curr_secs, adapt_step_flag,                  &
              !------------------------------------------------------
              ids,ide ,jds,jde ,kds,kde ,&
              ims,ime ,jms,jme ,kms,kme ,&
              i_start,i_end,j_start,j_end,kts,kte,num_tiles         &
              )
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  DO ALL PROCESSES THAT ARE NONLOCAL, WITH INFLUENCE OF ADJACENT GRIDS
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
   INTEGER, INTENT(IN) ::  &
         ids,ide,          &! start/end index for i in domain
         jds,jde,          &! start/end index for j in domain
         kds,kde,          &! start/end index for k in domain
         ims,ime,          &! start/end index for i in memory
         jms,jme,          &! start/end index for j in memory
         kms,kme,          &! start/end index for k in memory
         kts,kte            ! start/end index for k in tile

   INTEGER, INTENT(IN) ::  &
         num_tiles          ! number of tiles
   INTEGER, INTENT(IN), DIMENSION(num_tiles) :: &
         i_start,i_end,    &! start/end index for i in tile
         j_start,j_end      ! start/end index for j in tile

   INTEGER, INTENT(IN) ::  &
         itimestep,        &! model time step index
         sf_land_physics    ! scheme id for land process model

   INTEGER, INTENT(IN) ::  &
         no_flowdir,       &! number of flow direction
         nl_soil,          &! number of soil layers (>0) (CSSP)
         nl_snow,          &! max number of snow layers (<0) (CSSP)
         STEPSF             ! number of dynamic steps between SF calls

   LOGICAL, INTENT(IN), OPTIONAL :: &
         adapt_step_flag    ! flag to use adaptive time step
   REAL, INTENT(IN), OPTIONAL :: &
         sfdt,             &! time interval [minute] calling this surface_driver
         curr_secs          ! current time  [second]

   REAL, INTENT(IN)    ::  &
         DX,               &! x-dir grid spacing [m]
         DY,               &! y-dir grid spacing [m]
         DT                 ! dynamic model time step interval [seconds]

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         SC_WATER           ! surface category of water characteristics
                            ! 1=urban and built-up, 2=soil, 3=wetland, 4=land ice
                            ! 5=shallow lake, 6=deep lake,  7=sea ice, 8=ocean
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &                         ! halo 8
         ZSFCM,            &! modified surface elevation after filling sinks [m]
         xslpxs,           &! x-dir sdev slope for each grid [m/m]
         xslpys             ! y-dir sdev slope for each grid [m/m]

   REAL, INTENT(IN), DIMENSION(ims:ime,1:nl_soil,jms:jme) :: &               ! halo 8
         xbsw,             &! Clapp-Hornberger (b)
         xphi0,            &! soil minimum suction [mm]
         xporsl,           &! soil porosity (fraction that is voids) [-]
         xwdry,            &! soil dry or minimum residual water content [-]
         xhksm              ! soil hydraulic conductivity at saturation averaged in a layer [mm h2o/s]

   REAL, INTENT(IN), DIMENSION(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &       ! halo 8
         xdz                ! land layer thickness [m]

   REAL, INTENT(IN), DIMENSION(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &       ! halo 8
         xwice              ! land layer ice lens [kg/m2]

   REAL, INTENT(INOUT), DIMENSION(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &    ! halo 8
         xwliq              ! land layer liquid water [kg/m2]
!
!c CS
!!-xum   INTEGER, INTENT(IN), DIMENSION(ims:ime,8,jms:jme) :: &
!!+xum
!   INTEGER, INTENT(IN), DIMENSION(ims:ime,9,jms:jme) :: &
!!=xum
   INTEGER, INTENT(IN), DIMENSION(ims:ime,no_flowdir ,jms:jme) :: &
!. CS
         xdirid             ! 8-directional flow direction index             ! halo 8
                            !    +---+---+---+
                            !    | 6 | 7 | 8 |
                            !    +---+---+---+
                            !    | 5 | o | 1 |
                            !    +---+---+---+
                            !    | 4 | 3 | 2 |
                            !    +---+---+---+
                            ! note xdirid(i,l,j)>0 means inflow  from
                            !    the neighbor grid(l) to the center grid(i,j)
                            !    & xdirid(i,l,j)<0 means outflow from
                            !    the center grid(i,j) to the neighbor grid(l)

!  REAL, INTENT(IN), DIMENSION(ims:ime,9,jms:jme) :: &                       ! halo 8
!c CS
   REAL, INTENT(IN), DIMENSION(ims:ime,no_flowdir,jms:jme) :: &                       ! halo 8
!.CS
         xdirwd,           &! 8-directional channel width [m]
         xchleng,          &! 8-directional channel length [m]
         xbedslp            ! 8-directional channel bottom slope [m/m]

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &                         ! halo 8
         xrtws              ! surface water supply rate [mm/s]

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &                      ! halo 24
         xswdepth,         &! surface water depth [mm]
         xoutflow           ! surface outflow [m3/s]

!----------------------- Need no Halo ---------------------------------
   real, intent(inout), dimension(ims:ime,jms:jme) :: &
         xtsw,             &! total soil water below the water table [mm]
         xrsat,            &! excess saturation runoff [mm h2o/s]
         xrnof,            &! total runoff [mm h2o/s]
         xxerr              ! accumulated water  balance error [mm]

!----------------------- Local ----------------------------------------
   LOGICAL ::              &
         run_param          ! true for doing land surface calculation

   REAL ::                 &
         dtimesf            ! time interval [seconds] calling this surface_driver

   INTEGER :: its,ite,jts,jte,ij
!
!----------------------------------------------------------------------
!
   IF (sf_land_physics /= CSSPSCHEME) RETURN
!
!  ....Modified for adaptive time step.............................
   IF (PRESENT(sfdt) .AND. PRESENT(curr_secs) .AND.         &
       PRESENT(adapt_step_flag) .AND. (adapt_step_flag)) THEN
      run_param = (itimestep == 1) .or. (sfdt == 0) .or.    &
                  (curr_secs+DT >= (INT(curr_secs/(sfdt*60))+1)*sfdt*60)
   ELSE
      run_param = (itimestep == 1) .or. (STEPSF  == 1) .or. &
                  (MOD(itimestep,STEPSF) == 0)
   ENDIF
   IF (.NOT. (run_param) ) RETURN

!  ...Surface schemes need SF time step for updates
   dtimesf = DT * max(STEPSF,1)
   IF (PRESENT(sfdt)) THEN
      IF (sfdt == 0) THEN
         dtimesf= DT
      ELSEIF (PRESENT(adapt_step_flag) .and. (adapt_step_flag)) THEN
         call wrf_message("WARNING: When using an adaptive time-step, the surface module " // &
                          "time-step should be 0 (i.e., equivalent to model time-step). "  // &
                          "In order to proceed, the surface module time-step will be "     // &
                          "rounded to the nearest minute, possibly inaccuracy.")
         dtimesf= DT * NINT(sfdt*60/DT)
      ENDIF
   ENDIF
!  ................................................................
!
!=)Set start and end points for tiles
!
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij ,its,ite,jts,jte )

   DO ij  = 1 , num_tiles
      its = i_start(ij)
      ite = i_end(ij)
      jts = j_start(ij)
      jte = j_end(ij)
!
!----------------------------------------------------------------------
!  Do horizontal flow (lateralflow + surfaceflow)
!----------------------------------------------------------------------
!
      !--need halo 8  points for (SC_WATER,ZSFCM,xslpxs,xslpys
      !                          ,xdirid,xdirwd,xchleng,xbedslp
      !                          ,xbsw,xphi0,xporsl,xwdry,xhksm
      !                          ,xdz,xwice,xwliq,xrtws)
      !     & halo 24 points for (xswdepth,xoutflow)

         CALL wrf_debug(100,'in CSSP horizontalflow_vast')
         CALL horizontalflow_vast (                                 &
              SC_WATER,nl_soil,nl_snow,no_flowdir,dtimesf,DX,DY,    &
              ZSFCM,xdz,xhksm,xbsw,xphi0,xporsl,xwdry,              &
              xslpxs,xslpys,xwliq,xwice,xrsat,xrnof,xxerr,          &
              xdirid,xdirwd,xchleng,xbedslp,                        &
              xrtws,xswdepth,xoutflow,xtsw,                         &
              ids,ide ,jds,jde ,kds,kde ,&
              ims,ime ,jms,jme ,kms,kme ,&
              its,ite ,jts,jte ,kts,kte )



   ENDDO  ! end num_tiles
   !$OMP END PARALLEL DO

   END SUBROUTINE surface_finish
#endif


   SUBROUTINE surface_init(scheme_soilwater,                        &
              cdate   ,DT      ,SFDT    ,STEPSF  ,fdtcza  ,         &
              do_sfpr_avg      ,dt_sfpr_avg      ,                  &
              sf_sfclay_physics,sf_land_physics  ,sf_ocean_physics ,&
              bl_pbl_physics   ,DX      ,DY      ,                  &
              num_soil_layers  ,ISFFLX  ,XICE    ,XLAND   ,         &
              IVGTYP  ,ISLTYP  ,XLONG   ,XLAT    ,EMISS   ,MAVAIL  ,&
              VEGFRA  ,CANWAT  ,SNOW    ,SNOWH   ,SNOWF   ,         &
              TSK     ,TMN     ,TSLB    ,ZS      ,DZS     ,         &
              SMOIS   ,SMSTAV  ,SMSTOT  ,POTEVP  ,SFCEVP  ,         &
              SFRUNOFF,UDRUNOFF,ACSNOW  ,ACSNOM  ,SNOPCX  ,         &
              LOWLYR  ,UST     ,ZNT     ,PBLH    ,QCG     ,         &
              USTM    ,xtstar  ,QSFC    ,HFX     ,QFX     ,         &
              SH2O    ,SMFR3D  ,FNDSOILW,FNDSNOWH,                  &
              FRACTIONAL_SEAICE,RDMAXALB,ISWATER ,ISICE   ,MMINLU  ,&
              !--------for surface boundary conditions---------------
              setveg  ,schveg  ,setalb  ,                           &
              !--------for LSM---------------------------------------
              SNOALB  ,SNOTIME ,ALBBCK  ,                           &
              !--------for CSSP--------------------------------------
              nl_soil ,nl_snow ,sw_band ,                           &
              SC_WATER,SC_LANDU,                                    &
              BEDROCK ,FR_SAND ,FR_CLAY , OC_3D          ,          &
              xfrost_table,xzwt_perched , xqflx_sub_snow ,          &
#ifdef MHydro
              xwtslope,xfdshc  ,                                    &
#endif
#ifdef GenSWC
              xgswcA  ,xgswcB  ,xgswcC  ,xgswcD  ,xvwc10  ,xvwchs  ,&
#endif
              xbsw    ,                                             &
              xcsol   ,xporsl  ,xphi0   ,xdkmg   ,xdksatu ,         &
              xdkdry  ,xhksati ,xhksm   ,xrootfr ,xwdry   ,         &
#ifdef DoVAST
              xswdepth,xoutflow,                                    &
#endif
#ifdef LXZALB
              xsalf   ,                                             &
#else
              xalbsol ,                                             &
#endif

#ifdef NewLAKE
              xfrlak, xdplak,                                       & !+xum in
              xz0ma,  xstke1, xustar,                               &
              xzlake, xdzlak, xlktmp, xlkice,                       &
#endif

              xtss    ,xwliq   ,xwice   ,xz      ,xdz     ,xzwt    ,&
              xtg     ,xtlsun  ,xtlsha  ,xsag    ,xsigf   ,         &
              xfveg   ,xlai    ,xsai    ,xgreen  ,                  &
              xalb    ,xssun   ,xssha   ,xtsw    ,xwt     ,         &
              xfsun   ,xthermk ,xextkb  ,xextkd  ,                  &
              coszan  ,eccfn   ,xalbg   ,xalbv   ,xolrg   ,         &
              xrainc  ,xrainnc ,xsnowc  ,xsnownc ,xxerr   ,xzerr   ,&
              !--------optional MYNNSFC------------------------------
              mynn_closure_order,                                   &
              !--------optional SSTSKIN------------------------------
              skin_sst,dtwsk   ,                                    &
#ifdef CalUCM
              !--------optional URBAN -------------------------------
              sf_urban_physics,                                     &
              frc_urb2d,utype_urb2d,                                &
              num_roof_layers,num_wall_layers,num_road_layers,      &
              dzr,dzb,dzg,                                          &
              tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,         &
              trl_urb3d,tbl_urb3d,tgl_urb3d,                        &
              xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,          &
              num_urban_layers,                                     &
              trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,              &
              sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,            &
              a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,              &
              b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,              &
              dlg_bep,dl_u_bep,sf_bep,vl_bep,                       &
#endif
              !--------optional OML ---------------------------------
              omlcall,oml_hml0,tml,t0ml,hml,h0ml,huml,hvml,         &
              !------------------------------------------------------
#ifdef XOML
!^^^^^^^^^^^ FOR XOML ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       couple_oml ,oml_tskin  ,oml_vgrdopt,oml_mdldepth,oml_stretch, &
       oml_dt     ,oml_initopt,oml_rlaxopt,oml_rlaxtau ,             &
       oml_immpene,oml_isrpene,oml_jrvtype,oml_kbkspar ,             &
!       restart, cdate,             &
       nl_oml,   &
!       sc_water,xlat,                             &
       SSTB_1,                                                     &
       omlzly ,omlzlv ,omldz  , sst_nstep ,                        &
       t_oml,s_oml,u_oml,v_oml,e_oml,                   &  
       qcl_oml,buo_oml,gbo_oml,eps_oml,vtz_oml,ksz_oml,            &  
       kqz_oml,eln_oml,                    &  
       hobl   ,                                                    &  
       prevdsst, tempsst, obsprevdsst, obstempsst,                 &
!       nl_soil,nl_snow,  TSK     ,xtg     ,xtss  ,                 &
       dtw1,sstsk ,                                                &
        !grid_psfc,grid_swdns,grid_swups,grid_lwups,grid_lwdns,     &
!$$$$$$$$$$$ FOR XOML $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#endif

              restart ,allowed_to_read  ,&
              ids,ide ,jds,jde ,kds,kde ,&
              ims,ime ,jms,jme ,kms,kme ,&
              its,ite ,jts,jte ,kts,kte )

!----------------------------------------------------------------------
!--CONTROL VARIABLEs
!----------------------------------------------------------------------
   INTEGER, INTENT(IN) ::  &
         scheme_soilwater ! scheme id for soilwater calculation
   INTEGER, INTENT(IN) :: &
         ids,ide,     &! start/end index for i in domain
         jds,jde,     &! start/end index for j in domain
         kds,kde,     &! start/end index for k in domain
         ims,ime,     &! start/end index for i in memory
         jms,jme,     &! start/end index for j in memory
         kms,kme,     &! start/end index for k in memory
         its,ite,     &! start/end index for i in tile
         jts,jte,     &! start/end index for j in tile
         kts,kte       ! start/end index for k in tile
!
   LOGICAL, INTENT(IN) :: &
         restart,     &! true for restart
         allowed_to_read

   INTEGER, INTENT(IN) ::  &
         sf_sfclay_physics,&! scheme id for surface layer
         sf_land_physics,  &! scheme id for land process model
         sf_ocean_physics, &! scheme id for ocean process model
         bl_pbl_physics     ! scheme id for planetary boundary layer

   LOGICAL, INTENT(IN) ::  &
         do_sfpr_avg        ! true to do precip averaging for sfc
   REAL, INTENT(INOUT) ::  &
         dt_sfpr_avg        ! time period length in last precip averaging for sfc [s]

   INTEGER, INTENT(IN) ::  &
         ISFFLX,           &! 0 to turn off surface flux calculation (Not CSSP)
         FRACTIONAL_SEAICE,&! 1 to use fractional sea-ice
         ISWATER,          &! land use id for water category
         ISICE              ! land use id for land-ice category
 
   REAL, INTENT(IN)    ::  &
         DX,               &! x-dir grid spacing [m]
         DY,               &! y-dir grid spacing [m]
         DT,               &! dynamic model time step interval [seconds]
         SFDT,             &! time interval [minute] calling this surface_driver
         fdtcza             ! fraction of dtimesf shifted for cosza
                            ! 0-1: dtcza=dtimesf*fdtcza (0.5=mid of dtimesf)
                            !  >1: cosza averaged during dtimesf

   INTEGER, INTENT(INOUT) :: &
         STEPSF             ! number of dynamic steps between SF calls

   INTEGER, INTENT(IN) ::  &
         num_soil_layers    ! number of soil layer (LSM)

   REAL, INTENT(IN), DIMENSION(1:num_soil_layers) :: &
         ZS,               &! depths of centers of soil layers [m] (LSM,SLAB)
         DZS                ! thickneses of soil layers [m] (LSM,SLAB)

   LOGICAL, INTENT(IN) ::  &
         RDMAXALB           ! true to read in data for SNOALB

   CHARACTER(LEN=*), INTENT(IN) :: &
         MMINLU             ! land use classification type
!----------------------------------------------------------------------
!--COMMON VARIABLEs
!----------------------------------------------------------------------
   INTEGER, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         IVGTYP,           &! vegetation type (LSM)
         ISLTYP             ! soil type (LSM)

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         XLAND              ! land(1)/water(2) flag (Not CSSP)

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         XICE               ! sea ice (1) flag (LSM)

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         TSK,              &! surface skin temperature [K]
         TMN,              &! soil temperature at lower boundary [K] (LSM,SLAB)
         ZNT,              &! roughness length [m]
         UST,              &! u* in similarity theory [m/s]
         SNOW,             &! snow water equivalent [mm]
         SNOWF,            &! snow cover fraction [0-1]
         SNOWH,            &! snow depth [m]
         SNOALB,           &! snow albedo [-]
         ALBBCK,           &! background surface albedo [-] (LSM)
         CANWAT,           &! canopy water [mm] (lxz: m=>mm)
         EMISS,            &! surface emissivity [0-1] (SLAB)
         VEGFRA,           &! vegetation fraction [0-1]
         SMSTOT,           &! column total (liquid+ice) soil moisture content [mm]
         SMSTAV,           &! soil moisture availability for evapotranspiration [-]
         SNOTIME,          &! number of timesteps since last snowfall [-]
         SFRUNOFF,         &! accumulated surface runoff [mm]
         UDRUNOFF,         &! accumulated underground runoff [mm]
         ACSNOW,           &! accumulated snowfall [mm]
         ACSNOM,           &! accumulated melted snow [mm]
         SNOPCX,           &! accumulated snow phase change heat flux [W/m^2]
         POTEVP,           &! accumulated potential evaporation [mm]
         SFCEVP,           &! accumulated surface evaporation [mm]
         MAVAIL             ! surface moisture availability [0-1]   ! LSMRUC (IO) else (I)

   REAL, INTENT(INOUT), DIMENSION(ims:ime,1:num_soil_layers,jms:jme) :: &
         TSLB,             &! soil temperature [K] (LSM,SLAB)
         SMOIS,            &! volumetric soil moisture [m3/m3] (LSM)
         SH2O,             &! unfrozen volumetric soil moisture [m3/m3] (LSM)
         SMFR3D             ! soil ice

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         USTM,             &! u* in similarity theory [m/s] without vconv correction
         xtstar,           &! t* in similarity theory [K]
         QCG,              &! cloud water mixing ratio at the ground surface [kg/kg]
         QSFC,             &! specific humidity at lower boundary [kg/kg]
         HFX,              &! upward sensible heat flux at the surface [W/m^2]
         QFX                ! upward moisture flux at the surface [kg/m^2/s]


   INTEGER, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         LOWLYR             ! index of lowest model layer above ground

   LOGICAL, INTENT(IN) :: FNDSOILW, FNDSNOWH
!
   INTEGER, INTENT(IN) ::  &
         mynn_closure_order ! MYNNPBL tke closure order: 3=level 3; else=level 2.5

!----------------------------------------------------------------------
!--SSTSKIN VARIABLEs
!----------------------------------------------------------------------
   LOGICAL, OPTIONAL, INTENT(IN) ::  &
         skin_sst           ! true to predict skin sst

   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         dtwsk              ! warm sea surface temperature departure [C]
!
#ifdef CalUCM
!----------------------------------------------------------------------
!--UCM VARIABLEs
!----------------------------------------------------------------------
   INTEGER, INTENT(IN) ::  &
         sf_urban_physics   ! 1 or 2 to call UCM for urban effect

   INTEGER, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         utype_urb2d        ! urban type [nonurban=0, urban=1, suburban=2, rural=3]

   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         frc_urb2d          ! urban fraction of the grid [-]

   INTEGER, OPTIONAL, INTENT(IN) :: &
         num_roof_layers,  &! number of urban roof layers
         num_wall_layers,  &! number of urban wall layers
         num_road_layers    ! number of urban road layers

   REAL, OPTIONAL, INTENT(OUT), DIMENSION(1:num_soil_layers) :: &
         dzr,              &! urban roof layer thickness [cm]
         dzb,              &! urban building wall layer thickness [cm]
         dzg                ! urban ground layer thickness [cm]

   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         tr_urb2d,         &! urban roof temperature [K]
         tb_urb2d,         &! urban building wall temperature [K]
         tg_urb2d,         &! urban road temperature [K]
         tc_urb2d,         &! urban-canopy air temperature [K]
         qc_urb2d,         &! urban-canopy air mixing ratio [kg/kg]
         xxxr_urb2d,       &! urban roof   z/Monin-Obukhov-length [-]
         xxxb_urb2d,       &! urban wall   z/Monin-Obukhov-length [-]
         xxxg_urb2d,       &! urban road   z/Monin-Obukhov-length [-]
         xxxc_urb2d         ! urban-canopy z/Monin-Obukhov-length [-]

   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,1:num_soil_layers,jms:jme) :: &
         trl_urb3d,        &! urban roof layer temperature [K]
         tbl_urb3d,        &! urban wall layer temperature [K]
         tgl_urb3d          ! urban road layer temperature [K]
!
   INTEGER, INTENT(IN) :: num_urban_layers

   REAL, OPTIONAL, INTENT(INOUT), DIMENSION(ims:ime,1:num_urban_layers,jms:jme) :: &
         trb_urb4d        ,&!
         tw1_urb4d        ,&!
         tw2_urb4d        ,&!
         tgb_urb4d        ,&!
         sfw1_urb3d       ,&!
         sfw2_urb3d       ,&!
         sfr_urb3d        ,&!
         sfg_urb3d          !

       ! implicit and explicit components of the source and sink terms at each levels
       ! fluxes can be computed as: fx = a*x + b   example: t_fluxes = a_t * pt + b_t
   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: &
         a_u_bep          ,&! implicit component for momemtum in x-dir (center)
         a_v_bep          ,&! implicit component for momemtum in y-dir (center)
         a_t_bep          ,&! implicit component for temperature
         a_q_bep          ,&! implicit component for humidity
         a_e_bep          ,&! implicit component for TKE
         b_u_bep          ,&! explicit component for momemtum in x-dir (center)
         b_v_bep          ,&! explicit component for momemtum in y-dir (center)
         b_t_bep          ,&! explicit component for temperature
         b_q_bep          ,&! explicit component for humidity
         b_e_bep          ,&! explicit component for TKE
         vl_bep           ,&! volume of the urban grid cells
         dlg_bep          ,&! height above ground (l_ground in formula (24) of the BLM paper)
         sf_bep           ,&! surface of the urban grid cells
         dl_u_bep           ! length scale (l_b in formula (22) of the BLM paper)
#endif
!----------------------------------------------------------------------
!--OML VARIABLEs
!----------------------------------------------------------------------
   INTEGER, OPTIONAL, INTENT(IN) :: omlcall
   REAL, OPTIONAL, INTENT(IN) :: oml_hml0
   REAL, OPTIONAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         tml,t0ml,hml,h0ml,huml,hvml
!----------------------------------------------------------------------
!--CSSP VARIABLEs
!----------------------------------------------------------------------
   INTEGER, INTENT(IN) :: &
         nl_soil,          &! number of soil layers (>0)
         nl_snow,          &! max number of snow layers (<0)
         sw_band            ! number of solar bands (visible,near-infrared;direct,diffuse)

   CHARACTER(LEN=*), INTENT(IN) :: &
         cdate              ! current model date ('YYYY-MM-DD HH:MM:%s.ffff')

   LOGICAL, INTENT(IN) :: &
         setveg,           &! true to update veg b.c.
         setalb             ! true to update surface albedos

   INTEGER, INTENT(IN) :: &
         schveg             ! scheme id for vegetation variation:
                            ! 1=linear-in-time interpolation from wrfveg input
                            ! 2=empirical vegetation changes with temperature
                            ! 3=Ecodynamic model

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         SC_WATER           ! surface category of water characteristics
                            ! 1=urban and built-up, 2=soil, 3=wetland, 4=land ice
                            ! 5=shallow lake, 6=deep lake,  7=sea ice, 8=ocean

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         SC_LANDU           ! surface category of land use (vegetation)

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         XLAT,             &! latitude [degrees]
         XLONG,            &! longitude [degrees]
#ifdef MHydro
         xwtslope,         &! water table slope [m/m]
#endif
         BEDROCK            ! soil depth to bedrock [m]

   REAL, INTENT(IN), DIMENSION(ims:ime,1:nl_soil,jms:jme) :: &
         OC_3d,            &! organic carbon percent of weight 
         FR_SAND,          &! soil fraction of sand
         FR_CLAY            ! soil fraction of clay

   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         xfveg,            &! fraction of vegetation cover [-]
       ! xscv,             &! snow cover, water equivalent [mm] = SNOW
         xtg                ! ground surface temperature [K]

#ifdef LXZALB
   REAL, INTENT(IN), DIMENSION(ims:ime,sw_band,jms:jme) :: &
         xsalf              ! soil albedo localization factor [-]
                            ! (visible,near-infrared;direct,diffuse)
#else
   REAL, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: &
         xalbsol            ! soil albedo (visible beam) at saturation [-]
#endif

   REAL, INTENT(OUT), DIMENSION(ims:ime,1:nl_soil,jms:jme) :: &
         xcsol,            &! soil heat capacity of solids [J/m3/K]
         xporsl,           &! soil porosity (fraction that is voids) [-]
         xwdry,            &! soil dry or minimum residual water content [-]
         xphi0,            &! soil minimum suction (=- negative water potential) [mm]
#ifdef GenSWC
         xgswcA,           &! GenSWC coef A [mm]
         xgswcB,           &! GenSWC coef B [-]
         xgswcC,           &! GenSWC coef C [-]
         xgswcD,           &! GenSWC coef D [-]
         xvwc10,           &! volumetric water content (vwc) at 10 kPa [-]
         xvwchs,           &! volumetric water content (vwc) at hygroscopic [-]
#endif
         xbsw,             &! soil clapp and hornbereger (b) parameter [-]
         xdkmg,            &! soil thermal conductivity of minerals   [W/m/K]
         xdksatu,          &! soil thermal conductivity at saturation [W/m/K]
         xdkdry,           &! soil thermal conductivity when dry      [W/m/K]
         xhksati,          &! soil hydraulic conductivity at saturation [mm h2o/s]
         xhksm,            &! soil hydraulic conductivity at saturation averaged in a layer [mm h2o/s]
         xrootfr            ! soil fraction of vegetation root [-]

   REAL, INTENT(INOUT), DIMENSION(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &
         xtss,             &! land layer temperature [K]
         xwliq,            &! land layer liquid water [kg/m2]
         xwice,            &! land layer ice lens [kg/m2]
         xz,               &! land layer node depth [m]
         xdz                ! land layer thickness [m]

   REAL, INTENT(INOUT), DIMENSION(ims:ime,jms:jme) :: &
         xgreen,           &! greenness (fraction of green not dead leaf)
         xlai,             &! leaf area index
         xsai               ! stem area index

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         xtlsun,           &! sunlit leaf temperature [K]
         xtlsha,           &! shaded leaf temperature [K]
       ! xldew,            &! depth of water on foliage [mm] = CANWAT
         xsag,             &! dimensionless snow age [-]
       ! xsnowdp,          &! snow depth [m] = SNOWH
       ! xfsno,            &! fraction of snow cover on ground = SNOWF
         xsigf              ! fraction of veg. cover free of snow [-]

   REAL, INTENT(OUT), DIMENSION(ims:ime,sw_band,jms:jme) :: &
         xalb,             &! total averaged albedo [-]
         xalbg,            &! ground albedo [-]
         xalbv,            &! vegetation albedo [-]
         xssun,            &! sunlit canopy absorption of solar radiation
         xssha              ! shaded canopy absorption of solar radiation
                            ! (1->4): (visible, near-infrared) {direct, diffuse}

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         xfsun,            &! fraction of sunlit canopy
         xthermk,          &! canopy gap fraction for infrared radiation
         xextkb,           &! (k, g(mu)/mu) direct solar extinction coefficient
         xextkd             ! diffuse and scattered diffuse PAR extinction coefficient

#ifdef DoVAST
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         xswdepth,         &! surface water depth [mm]
         xoutflow           ! surface outflow [m3/s]
#endif

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
         coszan             ! cosine of solar zenith angle at next cssp step
   REAL, INTENT(OUT) :: &
         eccfn              ! earth-sun distance factor <(1/r)^2> at next cssp step

   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme) :: &
#ifdef MHydro
         xfdshc,           &! decay constant of saturated hydraulic conductivity [1/m]
#endif
         PBLH,             &! planetary boundary layer height [m]
         xolrg,            &! outgoing long-wave radiation from ground+canopy [W/m2]
         xzwt,             &! water talble depth [m] !-choi2
         xwt,              &! total water storage (unsaturated soil water + groundwater) (mm)
         xtsw,             &! total soil water below the water table [mm]
         xrainc,           &! accumulated convective rainfall at previous CSSP step [mm]
         xrainnc,          &! accumulated resolved   rainfall at previous CSSP step [mm]
         xsnowc,           &! accumulated convective snowfall at previous CSSP step [mm]
         xsnownc,          &! accumulated resolved   snowfall at previous CSSP step [mm]
         xxerr,            &! accumulated water  balance error [mm]
         xzerr              ! accumulated energy balance error [W/m2]
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme):: xqflx_sub_snow  !sublimation rate from snow pack (mm H2O /s) [+]
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme):: xfrost_table    ! frost table depth (m)
   REAL, INTENT(OUT), DIMENSION(ims:ime,jms:jme):: xzwt_perched    ! perched water table depth (m)
!
#ifdef XOML
!^^^^^^^^^^^ FOR XOML ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   integer  ,intent(in):: couple_oml , oml_tskin  , oml_vgrdopt , oml_rlaxopt
   integer  ,intent(in):: oml_initopt, oml_immpene, oml_isrpene, oml_jrvtype
   real     ,intent(in):: oml_mdldepth,oml_stretch, oml_dt , oml_rlaxtau, oml_kbkspar
   integer  ,intent(in)::nl_oml
   real,dimension(ims:ime,       jms:jme),intent(in)::     &
         SSTB_1
   real,dimension(nl_oml),intent(inout):: omlzly ,omlzlv ,omldz
   integer,intent(inout)::sst_nstep
   real,dimension(ims:ime,nl_oml,jms:jme),intent(inout)::  &
         t_oml,s_oml,u_oml,v_oml,e_oml ,         &
         qcl_oml,buo_oml,gbo_oml,eps_oml,vtz_oml,ksz_oml,  &
         kqz_oml,eln_oml
   real,dimension(ims:ime,       jms:jme),intent(inout)::  &
         hobl   
   real,dimension(ims:ime,       jms:jme),intent(inout)::  &
         dtw1,  sstsk
   real,dimension(ims:ime,       jms:jme),intent(inout)::  &
         prevdsst, tempsst, obsprevdsst, obstempsst
 
!$$$$$$$$$$$ FOR XOML $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#endif

#ifdef NewLAKE
    real, intent(in ), dimension(ims:ime,jms:jme) :: &
        xfrlak,      &! lake area fraction [%]
        xdplak        ! depth of lake [m]
    real, intent(out), dimension(ims:ime,jms:jme) :: &
        xz0ma,       &!
        xstke1,      &!
        xustar        !
    !real, intent(out), dimension(ims:ime,1:nl_lake,jms:jme) :: &
    real, intent(out), dimension(ims:ime,1:nl_soil,jms:jme) :: &
        xzlake,      &!
        xdzlak,      &!
        xlktmp,      &!
        xlkice        !
#endif


   LOGICAL ::              &
         dosfclay,         &! true for sfclay calculation
         dosfcflx,         &! true for surface flux calculation
         doocndyn,         &! true for ocean dynamic calculation
         getveg             ! true for updating vegetation cover

#ifdef CalUCM
   LOGICAL :: dourban       ! true for calling urban canopy model
   INTEGER :: ucmdoit       ! 1 for actually calling UCM
#endif
!
   INTEGER ::              &
         schscw,           &! scheme id for surface water category
                            ! 1=SC_WATER
                            ! 2=IVGTYP
                            ! 3=XLAND
         on_lake,          &! 1=update TSK over in-land lake as well
         Ic_water,         &! in-land water category in IVGTYP (USGS=16)
         ic_ocean,         &! ocean category in IVGTYP (USGS=0)
         ic_ice             ! ice category in IVGTYP (USGS=24)

   REAL ::                 &
         fc_land,          &! threshold for XLAND to be a land grid (defult=1.5)
         fc_ice,           &! threshold for fractional cover to be an ice grid (defult=0.5)
         xice_threshold,   &! threshold for sea-ice grid
         dtimesf            ! time interval [seconds] calling this surface_driver

   REAL, DIMENSION(ims:ime,jms:jme) :: &
         sicthk             ! sea-ice thickness [m]
!
   integer i,j
   integer itf,jtf,ktf
!======================================================================
!
   itf = min(ite,ide-1)
   jtf = min(jte,jde-1)
   ktf = min(kte,kde-1)

!--Calculate SFC time step

   STEPSF = nint(SFDT*60./DT)
   STEPSF = max(STEPSF,1)
!
!----------------------------------------------------------------------
!  Turn on all switches
!----------------------------------------------------------------------
!
   IF (.NOT. (sf_land_physics/=0 .OR. sf_ocean_physics/=0) ) RETURN
!
   dtimesf = DT * max(STEPSF,1)
!
   dosfclay = .true.
   IF (sf_land_physics==   CSSPSCHEME)  dosfclay = .false.
   IF ((bl_pbl_physics== MYJPBLSCHEME) .OR. &
       (bl_pbl_physics== MYJURBSCHEME) .OR. &
       (bl_pbl_physics==QNSEPBLSCHEME)) dosfclay = .true.
       !-enable {MYJSFC,QNSESFC} to set (THZ0,QZ0,UZ0,VZ0)

   IF (dosfclay) THEN
#ifndef LaxSFC
      IF (bl_pbl_physics   == MYJPBLSCHEME .AND. &
          sf_sfclay_physics/= MYJSFCSCHEME) THEN  ! it is relaxed by dosfcvisc
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/pbl = ', &
                 sf_sfclay_physics,bl_pbl_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (bl_pbl_physics   == MYJURBSCHEME .AND. &
          sf_sfclay_physics/= MYJSFCSCHEME) THEN  ! it is relaxed by dosfcvisc
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/pbl = ', &
                 sf_sfclay_physics,bl_pbl_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (sf_land_physics  == RUCLSMSCHEME .AND. &
          sf_sfclay_physics/= MYJSFCSCHEME) THEN  ! it is relaxed
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/land = ', &
                 sf_sfclay_physics,sf_land_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (sf_land_physics  ==  PXLSMSCHEME .AND. &
          sf_sfclay_physics/=  PXSFCSCHEME) THEN  ! it is relaxed
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/land = ', &
                 sf_sfclay_physics,sf_land_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
#endif
      IF (sf_sfclay_physics==MYNNSFCSCHEME .AND. &
          bl_pbl_physics   /=MYNNPBLSCHEME) THEN  ! it is required to set (tsq,qsq,cov)
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/pbl = ', &
                 sf_sfclay_physics,bl_pbl_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
      IF (sf_sfclay_physics==0 .AND. &
         (sf_land_physics  ==   SLABSCHEME .OR. &
          sf_land_physics  ==     GFDLSLAB .OR. &
          sf_land_physics  == NMMLSMSCHEME .OR. &
          sf_land_physics  ==    LSMSCHEME)) THEN
          WRITE( wrf_err_message , * ) 'xxxx NOT match: sfclay/land = ', &
                 sf_sfclay_physics,sf_land_physics
          CALL wrf_error_fatal ( wrf_err_message )
      ENDIF
   ENDIF

   dosfcflx = ISFFLX/=0
   IF (sf_land_physics==CSSPSCHEME) dosfcflx = .true.

   doocndyn = .false.
   IF (sf_land_physics==CSSPSCHEME .AND. &
       sf_ocean_physics/=0) doocndyn = .true.

   getveg = setveg
   IF (sf_land_physics==SLABSCHEME .OR. &
       sf_land_physics==GFDLSLAB) getveg = .false.

   xice_threshold = 0.5
   IF (FRACTIONAL_SEAICE==1) xice_threshold = 0.02

   fc_ice   = xice_threshold
   fc_land  = 1.5
   ic_ice   = ISICE
   ic_water = ISWATER
   ic_ocean = ic_water   ! IVGTYP uses same ISWATER to define all water bodies
   on_lake  = 1  ! update TSK over in-land lake as well as ocean
   schscw   = 3  ! SLAB or PXLSM uses XLAND  to define surface water category
   IF (dosfcflx) THEN
      IF (sf_land_physics==CSSPSCHEME) THEN
#if (defined NewLAKE || defined OldLAKE)
         on_lake  = 0    ! CSSP has a built-in predictive lake model
#else
         on_lake  = 1    ! +xum
#endif
         schscw   = 1    ! CSSP uses SC_WATER to define surface water category
         ic_ocean = 0    ! SC_WATER uses 0   to define ocean
      ELSEIF (sf_land_physics==LSMSCHEME .OR. &
              sf_land_physics==NMMLSMSCHEME .OR. &
              sf_land_physics==RUCLSMSCHEME) THEN
         schscw   = 2    ! LSM uses IVGTYP   to define surface water category
      ENDIF
   ENDIF

#ifdef CalUCM
   ucmdoit = 0
   dourban = ( &
   !--optional io for urban canopy model
      (sf_urban_physics==UCMSCHEME     .AND. num_soil_layers==4               .AND. &
      PRESENT(num_roof_layers)         .AND. num_roof_layers==num_soil_layers .AND. &
      PRESENT(num_wall_layers)         .AND. num_wall_layers==num_soil_layers .AND. &
      PRESENT(num_road_layers)         .AND. num_road_layers==num_soil_layers .AND. &
      PRESENT(dzr)       .AND. PRESENT(dzb)       .AND. PRESENT(dzg)          .AND. &
      PRESENT(tr_urb2d)                .AND. PRESENT(tb_urb2d)                .AND. &
      PRESENT(tg_urb2d)                .AND. PRESENT(tc_urb2d)                .AND. &
      PRESENT(qc_urb2d)                .AND.                                        &
      PRESENT(xxxr_urb2d)              .AND. PRESENT(xxxb_urb2d)              .AND. &
      PRESENT(xxxg_urb2d)              .AND. PRESENT(xxxc_urb2d)              .AND. &
      PRESENT(trl_urb3d)               .AND. PRESENT(tbl_urb3d)               .AND. &
      PRESENT(tgl_urb3d)                                                   )  .OR.  &
   !--optional io for urban multi-layer model
      (sf_urban_physics==BEPSCHEME     .AND. num_urban_layers>=4              .AND. &
      PRESENT(trb_urb4d)               .AND. PRESENT(tw1_urb4d)               .AND. &
      PRESENT(tw2_urb4d)               .AND. PRESENT(tgb_urb4d)               .AND. &
      PRESENT(sfw1_urb3d)              .AND. PRESENT(sfw2_urb3d)              .AND. &
      PRESENT(sfr_urb3d)               .AND. PRESENT(sfg_urb3d)               .AND. &
      PRESENT(a_u_bep)                 .AND. PRESENT(b_u_bep)                 .AND. &
      PRESENT(a_v_bep)                 .AND. PRESENT(b_v_bep)                 .AND. &
      PRESENT(a_t_bep)                 .AND. PRESENT(b_t_bep)                 .AND. &
      PRESENT(a_q_bep)                 .AND. PRESENT(b_q_bep)                 .AND. &
      PRESENT(a_e_bep)                 .AND. PRESENT(b_e_bep)                 .AND. &
      PRESENT(dlg_bep)                 .AND. PRESENT(dl_u_bep)                .AND. &
      PRESENT(sf_bep)                  .AND. PRESENT(vl_bep)               ) ).AND. &
   !--optional io for all urban models
      PRESENT(utype_urb2d)             .AND. PRESENT(frc_urb2d)
   IF (sf_urban_physics==UCMSCHEME .OR. sf_urban_physics==BEPSCHEME) THEN
      IF (dourban) THEN
         ucmdoit = sf_urban_physics
      ELSE
         print*,'**** Warning: Inappropriate IO => Abort calling UCM'
      ENDIF
   ENDIF
#endif
!
!  ==================
   IF (dosfclay) THEN
!  ==================
!
!----------------------------------------------------------------------
!  Initialize surface layer schemes
!----------------------------------------------------------------------
!
   sfclay_select: SELECT CASE(sf_sfclay_physics)

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (SFCLAYSCHEME)
         CALL wrf_debug(300,'surface_init: sfclayinit')
         CALL sfclayinit(ZNT,UST,USTM,xtstar,QSFC,HFX,QFX,              &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (MYNNSFCSCHEME)
         CALL wrf_debug(300,'surface_init: mynnsfcinit')
         CALL mynnsfcinit(ZNT,UST,USTM,xtstar,HFX,QFX,                  &
                  restart, allowed_to_read ,mynn_closure_order,         &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (PXSFCSCHEME)
         CALL wrf_debug(300,'surface_init: pxsfclayinit')
         CALL pxsfclayinit(ZNT,UST,HFX,QFX,                             &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (MYJSFCSCHEME)
         CALL wrf_debug(300,'surface_init: myjsfcinit')
         CALL myjsfcinit(LOWLYR,UST,ZNT,XLAND,XICE,IVGTYP,              &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (QNSESFCSCHEME)
         CALL wrf_debug(300,'surface_init: qnsesfcinit')
         CALL qnsesfcinit(LOWLYR,UST,ZNT,XLAND,XICE,IVGTYP,             &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (GFSSFCSCHEME)
         CALL wrf_debug(300,'surface_init: sf_gfsinit')
         CALL sf_gfsinit(XLAND,ZNT,                                     &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte  )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE DEFAULT
         CALL wrf_debug(300,'xxxx The surface layer scheme does not exist')

   END SELECT sfclay_select

!  =====
   ENDIF
!  =====
!
!  ==================
   IF (dosfcflx) THEN
!  ==================
!
!----------------------------------------------------------------------
!  Initialize land surface models
!----------------------------------------------------------------------
!
   land_select: SELECT CASE(sf_land_physics)

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (CSSPSCHEME)
         CALL wrf_debug(300,'surface_init: csspinit')
         CALL csspinit (scheme_soilwater,                               &
                  cdate   ,dtimesf ,DX      ,DY      ,                  &
                  nl_soil ,nl_snow ,sw_band ,fdtcza  ,                  &
                  XLONG   ,XLAT    ,SC_WATER,SC_LANDU,                  &
                  BEDROCK ,FR_SAND ,FR_CLAY ,OC_3d   ,                  &
                  xfrost_table,xzwt_perched , xqflx_sub_snow,           &
#ifdef MHydro
                  xwtslope,xfdshc  ,                                    &
#endif
#ifdef GenSWC
                  xgswcA  ,xgswcB  ,xgswcC  ,xgswcD  ,xvwc10  ,xvwchs  ,&
#endif
                  xbsw    ,                                             &
                  xcsol   ,xporsl  ,xphi0   ,xdkmg   ,xdksatu ,         &
                  xdkdry  ,xhksati ,xhksm   ,xrootfr ,xfveg   ,xwdry   ,&
#ifdef DoVAST
                  xswdepth,xoutflow,                                    &
#endif
#ifdef LXZALB
                  xsalf   ,                                             &
#else
                  xalbsol ,                                             &
#endif

#ifdef NewLAKE
                  xfrlak, xdplak,                                       & ! +xum in
                  xz0ma,  xstke1, xustar,                               &
                  xzlake, xdzlak, xlktmp, xlkice,                       &
#endif

!CS
                  xtss    ,xwliq   ,xwice   ,xz      ,xdz     ,xzwt    ,&
!!+ygan
!                  xtss    ,xwliq   ,xwice   ,xz      ,xdz     ,         &
!!.ygan
!.CS
                  xtg     ,xtlsun  ,xtlsha  ,CANWAT  ,xxerr   ,xzerr   ,&
                  SNOWH   ,SNOW    ,xsag    ,SNOWF   ,xsigf   ,PBLH    ,&
                  xlai    ,xsai    ,xgreen  ,schveg  ,getveg  ,setalb  ,&
                  xalb    ,xssun   ,xssha   ,xtsw    ,xwt     ,         &
                  xfsun   ,xthermk ,xextkb  ,xextkd  ,                  &
                  EMISS   ,coszan  ,eccfn   ,xalbg   ,xalbv   ,xolrg   ,&
                  VEGFRA  ,                                             &
#ifdef CalUCM
                  ucmdoit ,                                             &
#endif
                  restart ,allowed_to_read  ,&
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,ite ,jts,jte ,kts,kte )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (NMMLSMSCHEME)
         CALL wrf_debug(300,'surface_init: nmmlsminit')
         CALL nmmlsminit(                                               &
                  VEGFRA,SNOW,CANWAT,SMSTAV,SMSTOT,                     &
                  SFRUNOFF,UDRUNOFF,ACSNOW,ACSNOM,                      &
                  IVGTYP,ISLTYP,TSLB,SMOIS,DZS,                         &
                  POTEVP,SNOPCX,SNOWH,TSK,                              &
                  TMN,XLAND,XICE,                                       &
                  num_soil_layers,           &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (LSMSCHEME)
         CALL wrf_debug(300,'surface_init: lsminit')
         CALL lsminit (                                                 &
                  VEGFRA,SNOW,SNOWF,SNOWH,CANWAT,SMSTAV,SMSTOT,         &
                  SFRUNOFF,UDRUNOFF,ACSNOW,ACSNOM,                      &
                  IVGTYP,ISLTYP,TSLB,SMOIS,SH2O,DZS,                    &
                  FNDSOILW,FNDSNOWH,                                    &
                  MMINLU,RDMAXALB,SNOALB,SNOTIME,xzerr,                 &
                  XLAND,XICE,SNOPCX,POTEVP,xice_threshold,              &
                  num_soil_layers,                                      &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (RUCLSMSCHEME)
         CALL wrf_debug(300,'surface_init: lsmrucinit')
         CALL lsmrucinit(                                               &
                  SMFR3D,TSLB,SMOIS,ISLTYP,MAVAIL,                      &
                  SFRUNOFF,UDRUNOFF,ACSNOW,ACSNOM,                      &
                  num_soil_layers,                                      &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (PXLSMSCHEME)

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE (SLABSCHEME)
         CALL wrf_debug(300,'surface_init: slabinit')
         CALL slabinit(                                                 &
                  TSK,TMN,TSLB,ZS,DZS,                                  &
                  num_soil_layers,                                      & 
                  omlcall, oml_hml0,tml,t0ml,hml,h0ml,huml,hvml,        &
                  restart, allowed_to_read,  &
                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte )

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CASE DEFAULT
         CALL wrf_debug(300,'xxxx The land surface model does not exist')

   END SELECT land_select

!  =====
   ENDIF
!  =====
!
#ifdef CalUCM
!  =================
   IF (dourban) THEN
!  =================
!
!----------------------------------------------------------------------
!  Initialize land surface models
!----------------------------------------------------------------------
!
         CALL wrf_debug(300,'surface_init: ucminit')
         CALL ucminit(                                                  &
                  IVGTYP,TSK,TSLB,                                      &
                  ucmdoit,                                              &
                  FRC_URB2D,UTYPE_URB2D,                                &
                  num_roof_layers,num_wall_layers,num_road_layers,      &
                  num_soil_layers,DZR,DZB,DZG,                          &
                  TR_URB2D,TB_URB2D,TG_URB2D,TC_URB2D,QC_URB2D,         &
                  XXXR_URB2D,XXXB_URB2D,XXXG_URB2D,XXXC_URB2D,          &
                  TRL_URB3D,TBL_URB3D,TGL_URB3D,                        &
                  num_urban_layers,                                     &
                  trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,              &
                  sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d,            &
                  a_u_bep,a_v_bep,a_t_bep,a_q_bep,a_e_bep,              &
                  b_u_bep,b_v_bep,b_t_bep,b_q_bep,b_e_bep,              &
                  dlg_bep,dl_u_bep,sf_bep,vl_bep,                       &
                  restart ,allowed_to_read  ,&
                  ids,ide ,jds,jde ,kds,kde ,&
                  ims,ime ,jms,jme ,kms,kme ,&
                  its,ite ,jts,jte ,kts,kte )

!  =====
   ENDIF
!  =====
#endif
!
!  ===============================================================
   IF (PRESENT(skin_sst) .AND. PRESENT(dtwsk) .AND. skin_sst) THEN
!  ===============================================================
      IF (.NOT.restart) THEN
         DO j = jts,jtf
         DO i = its,itf
            dtwsk(i,j) = 0.
         ENDDO
         ENDDO
      ENDIF
!  =====
   ENDIF
!  =====
!
#ifdef XOML
!^^^^^^^^^^^ FOR XOML ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   CALL XOML_init(                                                  &
        couple_oml,oml_tskin,oml_vgrdopt,oml_mdldepth,oml_stretch,  &
        oml_dt, oml_initopt,oml_rlaxopt,oml_rlaxtau,                &
        oml_immpene,oml_isrpene,oml_jrvtype,oml_kbkspar,            &
        restart, cdate,                                             &
        ids,ide, jds,jde, nl_oml,                                   &
        ims, ime, jms, jme,                                         &
        its, ite, jts, jte,                                         &
        sc_water,xlat ,                 &
        SSTB_1,                                                     &
  !grid_psfc,grid_swdns,grid_swups,grid_lwups,grid_lwdns,           &
  !grid_raincv,grid_rainncv, grid_qfx,  grid_hfx,grid_lh, grid_xtaux, grid_xtauy,     &
        omlzly ,omlzlv ,omldz  ,  sst_nstep ,                       &
        t_oml,s_oml,u_oml,v_oml,e_oml ,                   &
        qcl_oml,buo_oml,gbo_oml,eps_oml,vtz_oml,ksz_oml,            &
        kqz_oml,eln_oml,                    &
        hobl   ,                                   &
        prevdsst, tempsst, obsprevdsst, obstempsst,                 &
        nl_soil,nl_snow,  TSK      ,xtss   ,                &
        dtw1, sstsk )
!$$$$$$$$$$$ FOR XOML $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#endif

!  ==================
   IF (doocndyn) THEN
!  ==================
!
!----------------------------------------------------------------------
!  Initialize ocean surface model
!----------------------------------------------------------------------
!
   ocean_select: SELECT CASE(sf_ocean_physics)

      CASE (SOMSCHEME)
         CALL wrf_debug(300,'surface_init: sominit')
         CALL sominit(restart,allowed_to_read)

      CASE DEFAULT
         CALL wrf_debug(300,'xxxx The ocean surface model does not exist')

   END SELECT ocean_select

!  =====
   ENDIF
!  =====
!
   IF (restart) RETURN   ! exit for a restart (continuation) run
!
!----------------------------------------------------------------------
!  Initialize surface driver
!----------------------------------------------------------------------
!
   IF (do_sfpr_avg) dt_sfpr_avg = 0.

   DO j = jts,jtf
   DO i = its,itf
        xrainc(i,j) = 0.
       xrainnc(i,j) = 0.
        xsnowc(i,j) = 0.
       xsnownc(i,j) = 0.
        ACSNOW(i,j) = 0.
        ACSNOM(i,j) = 0.
        SNOPCX(i,j) = 0.
        SFCEVP(i,j) = 0.
        POTEVP(i,j) = 0.
      SFRUNOFF(i,j) = 0.
      UDRUNOFF(i,j) = 0.
           QCG(i,j) = 0.
   ENDDO
   ENDDO

   END SUBROUTINE surface_init


END MODULE module_surface_driver
