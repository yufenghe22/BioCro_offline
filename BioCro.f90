module BioCroX
implicit none

type :: biocro_outputs 
real,dimension(:),allocatable::cws 
real ::TTc,Grain,Leaf,Stem,Root,Rhizome,LeafN,Sp,&
       alpha,vmax,stomataWS,biolai,bioev,biotr,CA
end type biocro_outputs

type :: c4structure
real :: Assim,Gs,Ci,GrossAssim
end type c4structure

type :: CanACstructure
real :: Assim,Trans,GrossAssim,Idir,Idiff,Assim_Idir
end type CanACstructure

type :: lightMEstructure
real :: dir_frac,diff_frac
end type lightMEstructure

type :: sunMLstructure
real :: dir_irr,diff_irr,tot_irr,sunlit_frac,shade_frac,height
end type sunMLstructure

type :: nitroPstructure
real :: N=85.0, kLN=0.5, Vmaxb1=0.6938, Vmaxb0=-16.25, alphab1=0.000488,&
alphab0=0.02367,Rdb1=0.1247, Rdb0=-4.5917 
end type nitroPstructure

type :: ETstructure
real :: TransR,EPenman,EPriestly,Deltat,LayerCond
end type ETstructure

type :: soilMLstructure
real :: rcoefPhoto,drainage,Nleach,rcoefSpleaf,SoilEvapo
real,dimension(:),allocatable ::cws,hourlyWflux,rootdist
end type soilMLstructure

type :: soTexStructure
real::silt,clay,sand,air_entry,b,Ks,satur,fieldc,wiltp,bulkd
end type soTexStructure

type :: wsStructure
real::psim,rcoefPhoto,rcoefSpleaf,awc,runoff,Nleach
end type wsStructure

type :: biopartistructure
real ::kLeaf,kStem,kRoot,kRhiz,kGrain
end type biopartistructure

real,parameter :: pi=3.1415926

contains

real function AbiotEff(smoist,stemp)
    real :: smoist,stemp,TempEff,MoisEff
    if (stemp < 35.0) then
    TempEff = 1.0087 / (1.0 + (46.2 * exp(-0.1899 * stemp)))
    else
    TempEff = -0.0826 * stemp + 3.84
    end if
    MoisEff = 1.0267 / (1.0 + 14.7 * exp(-6.5 * smoist))
    AbiotEff=TempEff * MoisEff
end function AbiotEff

function ballBerry(Amu, Cappm, Temp, RelH, beta0, beta1) result(gsmol)
real,intent(in) :: Amu, Cappm, Temp, RelH, beta0, beta1
real :: aaa,bbb,ccc,ddd,hs,gswmol,acs,Cs,leafTk,pwi,pwaPa,Camf,assimn,wa,wi
real :: gsmol
	real,parameter:: gbw = 1.2 ! According to Collatz et al. (1992) pg. 526
	real,parameter :: ptotPa = 101325.0 ! Atmospheric pressure 
	leafTk = Temp + 273.15
	pwi = fnpsvp(leafTk)
	pwaPa = RelH * pwi
	Camf = Cappm * 1e-6
	assimn = Amu * 1e-6
  
 	wa  = pwaPa / ptotPa
	wi  = pwi / ptotPa

	if (assimn < 0.0) then
		gswmol = beta0 
    else
		Cs  = Camf - (1.4/gbw)*assimn
		if (Cs < 0.0) then
			print *,"Camf,assimn are",Camf,assimn
                        error stop "Cs is less than 0"
        end if
		acs = assimn/Cs

		if (acs < 1e-6) then
            acs = 1e-6		
        end if
		aaa = beta1 * acs
		bbb = beta0 + gbw - (beta1 * acs)
		ccc = -(wa / wi) * gbw - beta0

		ddd = bbb * bbb - 4.0*aaa*ccc

		hs  = (-bbb + sqrt(ddd)) / (2.0* aaa)

		gswmol = beta1 * hs * acs + beta0
    end if
	gsmol = gswmol * 1000.0 ! converting to mmol 

	if (gsmol <= 0.0) then
        gsmol = 1e-2
    end if
end function ballBerry

function fnpsvp(Tkelvin) result(esat)
	real, intent(in) :: Tkelvin
	real :: esat
	real :: tmp,u,v
! water boiling point = 373.16 oK
 ! This is the Arden Buck Equation 
 !http:!en.wikipedia.org/wiki/Arden_Buck_equation
	tmp = Tkelvin - 273.15
	u = (18.678 - tmp/234.5)*tmp
	v = 257.14 + tmp
	esat = (6.1121 * exp(u/v))/10.0
	!esat = esat/10.0
end function fnpsvp

function c4photoC(Qp, Tl, RH, vmax, alpha,& 
		       kparm, theta, beta,&
		       Rd, bb0, bb1, StomaWS, Ca, ws,upperT,lowerT) result(c4str)
			real, intent(in)::Qp, Tl, RH, vmax, alpha,kparm, theta, beta,&
		       Rd, bb0, bb1, StomaWS, Ca,upperT,lowerT
		    real:: Csurface,InterCellularCO2,KQ10,kT,Vtn,Vtd,VT,Rtn,Rtd,RT,&
		    b0,b1,b2,M1,M2,M,kT_IC_P,Quada,Quadb,Quadc,a21,a22,a2,Assim,csurfaceppm,&
		    Gs,diff,miC,GrossAssim,Assim0,Gs0,IntCO2,g,h,Assim1,Gs1,g1,dg,Ci2
			type(c4structure):: c4str
	real, parameter :: AP = 101325.0 !Atmospheric pressure According to wikipedia (Pa)*/
	real, parameter :: P = AP / 1e3 ! kPa */
	!const double PS = 38   Atmospheric pressure of CO2 */
	real,parameter :: Q10 = 2.0  ! Q10 increase in a reaction by 10 C temp */
	! Defining biochemical variables */
	real :: OldAssim, Tol
	integer :: iterCounter,ws
        OldAssim = 0.0
        Tol = 0.01

 ! When the method does not converge I do not use the iterative solution*/
 	! double Assim0 = 0.0 unused
 	! double Gs0 = 0.0 unused
 	! double IntCO2 = 0.0 unused
 	! partial pressure of CO2 at the leaf surface */
   
 	! if(StomaWS < 0.5) ws = 0 */

	Csurface = (Ca * 1e-6) * AP 
  
	InterCellularCO2 = Csurface * 0.4 ! Initial guestimate */

	KQ10 =  Q10**((Tl - 25.0) / 10.0)

	kT = kparm * KQ10

 	! First chunk of code see Collatz (1992) */
 !	Vtn = vmax * pow(2,((Tl-25.0)/10.0))
 !	Vtd = ( 1 + exp(0.3 * (3.0-Tl)) ) * (1 + exp( 0.3*(Tl-37.5) ))
 !	VT  = Vtn / Vtd
 
 !       This is the code implementing temperature limitations
      Vtn = vmax * 2.0**((Tl-25.0)/10.0)
      Vtd = ( 1.0 + exp(0.3 * (lowerT-Tl)) ) * (1.0 + exp( 0.3*(Tl-upperT) ))
      VT  = Vtn / Vtd


 	! Second chunk of code see Collatz (1992) */
	Rtn = Rd * 2.0**((Tl-25.0)/10.0)
	Rtd =  1.0 + exp( 1.3 * (Tl-55.0) ) 
	RT = Rtn / Rtd  

 	! Third chunk of code again see Collatz (1992) */
	b0 = VT * alpha  * Qp 
	b1 = VT + alpha  * Qp 
	b2 = theta 

 	! Calculate the 2 roots */
	M1 = (b1 + sqrt(b1*b1 - (4.0 * b0 * b2)))/(2.0*b2) 
	M2 = (b1 - sqrt(b1*b1 - (4.0 * b0 * b2)))/(2.0*b2) 

 	! This piece of code selects the smaller root */
    M=min(M1,M2)
 	! Here the iterations will start */
	iterCounter = 0
if(isnan(M)) print *, "M too large and T and Q are",M,Tl,Qp
	do while (iterCounter < 50)
		kT_IC_P = kT * (InterCellularCO2 / P*1000.0) !Yufeng
		Quada = M * kT_IC_P
		Quadb = M + kT_IC_P
		Quadc = beta 

		a21 = (Quadb - sqrt(Quadb*Quadb - (4.0 * Quada * Quadc))) / (2.0 * Quadc)
        a22 = (Quadb + sqrt(Quadb*Quadb - (4.0 * Quada * Quadc))) / (2.0 * Quadc)
        a2  =  min(a21,a22)
		Assim = a2 - RT !Yufeng

		if (ws == 0) then
            Assim = Assim * StomaWS 
        end if
!print *,"a2 M InterCellularCO2 Tl RH Qp StomaWS are",a2,M,InterCellularCO2,Tl,RH,Qp,StomaWS
 		! milimole per meter square per second*/
		csurfaceppm = Csurface * 10.0 
!        print *,"Assim,InterCellularCO2,a21,a22",Assim,InterCellularCO2,a21,a22,iterCounter
 		! Need to create the Ball-Berry function */
		Gs =  ballBerry(Assim,csurfaceppm, Tl, RH, bb0, bb1) 
                !Gs is [mmol m-2 s-1]
        if(isnan(Gs)) then
       print *,"Assim is too large",Assim,iterCounter
        print *,"csurfaceppm is",csurfaceppm
        print *,"Tl is",Tl
        print *, "RH is",RH
        print *, "bb0, bb1",bb0, bb1
       print *,"Tl,Qp,M are",Tl,Qp,M
        print *,"InterCellularCO2,a21,a22",InterCellularCO2,a21,a22 
        endif 
        !Stomata conductance.
	if (ws == 1) then 
            Gs = Gs * StomaWS 
        end if
        g = (Csurface - (Assim * 1.6 * AP) / (Gs * 1e3))-InterCellularCO2 !pa
        if (abs(g) < 1e-6) then
            diff = 0.
            exit
        else
        h = g
        kT_IC_P = kT * ((InterCellularCO2 + h) / P * 1000.)!kT: mol m-2 s-1.1e6 -> umol
		Quada = M * kT_IC_P !M:umol m-2 s-1
		Quadb = M + kT_IC_P
		Quadc = beta 
		a21 = (Quadb - sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a22 = (Quadb + sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a2  =  min(a21,a22)
	Assim1 = a2 - RT
        Gs1 =  ballBerry(Assim1,csurfaceppm, Tl, RH, bb0, bb1) 
        if (ws == 1) then
            Gs1 = Gs1 * StomaWS 
        endif
        g1=(Csurface - (Assim1 * 1.6 * AP) / (Gs1 * 1e3))-InterCellularCO2-h
        dg=(g1-g)/h ! dg: derivative of g(x) at Ci
        if (abs(dg) < 1e-6) then
            diff = 0.
           exit 
        endif
        Ci2 = InterCellularCO2 - g / dg
!		InterCellularCO2 = Csurface - (Assim * 1e-6 * 1.6 * AP) / (Gs * 0.001)

!		if (InterCellularCO2 < 0.0) then
!			InterCellularCO2 = 1e-5
!        end if

!         if(iterCounter == 0) then
!                        Assim0 = Assim
!                        Gs0 = Gs
!                        IntCO2 = InterCellularCO2
!         endif
!		diff = OldAssim - Assim
!		if (diff < 0.0) then
!            diff = -diff
!        end if
!		if (diff < Tol) then
!			exit
!        else
!            OldAssim = Assim
!        end if
         diff = Ci2 - InterCellularCO2
	if (abs(diff) < Tol) then
		exit
        else
            InterCellularCO2 = Ci2
        endif
		iterCounter=iterCounter+1
        endif !if else g==0
    end do !while

       
!Yufeng: This may not converge if stomaWS is near zero,which causes InterCO2 too large!
!For this, I added a condition for StomaWS (no less than 1e-2) and also when
!it's not converged, assign Assim to -RT only. Feel like the stomata is closed
!and only respiration exists.

!if (diff >Tol) then
!Assim=Assim0
!Gs = Gs0
!InterCellularCO2 = IntCO2 ! These assignments are problematic
!print *,"c4photoC did not converge,Assim0,Gs0,IntCO2 are",Assim0,Gs0,IntCO2
!endif
        kT_IC_P = kT * (InterCellularCO2 / P*1000.) !kT: mol m-2 s-1.1e6 -> umol
	Quada = M * kT_IC_P!M:umol m-2 s-1
	Quadb = M + kT_IC_P
	Quadc = beta 
	a21 = (Quadb - sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a22 = (Quadb + sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a2  =  min(a21,a22)
	Assim = a2 - RT
        if (ws == 0) then
            Assim = Assim * StomaWS 
        endif
        Gs =  ballBerry(Assim,csurfaceppm, Tl, RH, bb0, bb1)
	if (ws == 1) then
            Gs = Gs * StomaWS 
        endif
        InterCellularCO2 = Csurface - (Assim * 1.6 * AP) / (Gs * 1e3)

	miC = (InterCellularCO2 / AP) * 1e6
 
    if (Gs > 600.0) then
	    Gs = 600.0
    end if
    GrossAssim=Assim+RT
	c4str%Assim = Assim
	c4str%Gs = Gs
	c4str%Ci = miC
    c4str%GrossAssim=GrossAssim
end function c4photoC

function CanAC(LAI,DOY,cosz,solar_dir,solar_diff,Temp,RH,WindSpeed,lat,nlayers,&
    Vmax,Alpha,Kparm,beta,Rd,Catm,b0,b1,theta,kd,chil,heightf,leafN,kpLN,&
    lnb0,lnb1,lnfun,upperT,lowerT,nitroP,leafwidth,eteq,StomataWS,ws) result(CanACstr)
	real, parameter :: cf = 3600 * 1e-6 * 44 * 1e-6 * 10000 
!      For Assimilation 
!      3600 - seconds per hour 
!      1e-6 - micromoles per mole 
!      30 - grams per mole for CO2 ???? YH: 44???
!      1e-6 - megagrams per gram 
!      10000 - meters squared per hectare
!      Need a different conversion factor for transpiration 
	real, parameter :: cf2 = 3600 * 1e-3 * 18 * 1e-6 * 10000 
!      For Transpiration 
!      3600 - seconds per hour 
!      1e-3 - millimoles per mole 
!      18 - grams per mole for H2O 
!      1e-6 - megagrams per  gram 
!      10000 - meters squared per hectare
real, intent(in) :: LAI,solar_dir,solar_diff,Temp,RH,WindSpeed,lat,&
    Vmax,Kparm,beta,Catm,b0,b1,theta,kd,chil,heightf,leafN,kpLN,&
    lnb0,lnb1,upperT,lowerT,leafwidth,StomataWS,cosz
    real, intent(inout) :: Alpha,Rd
    integer, intent(in) :: DOY,nlayers,lnfun,eteq,ws
    type(CanACstructure) :: CanACstr
    type(lightMEstructure) :: lightME_tmp
    type(sunMLstructure),dimension(nlayers) :: light_profile 
! nitroP is a structure, which contains Vmaxb0,Vmaxb1,alphab0,alphab1,Rdb0,Rdb1
	type(nitroPstructure) :: nitroP
	type(c4structure) :: temp_photo_results
	type(ETstructure) :: tmp5_ET,tmp6_ET
    real, dimension(1:nlayers):: relative_humidity_profile, wind_speed_profile,leafN_profile
    real :: Idir1,Idiff,cosTh,LAIc,leafN_lay,vmax1,rhi,layerWindSpeed,Idir2,Itot,pLeafsun,&
            CanHeight,Leafsun,TempIdir,AssIdir,GAssIdir,Idiff2,pLeafshade,Leafshade,TempIdiff,&
            AssIdiff,GAssIdiff,solarR
    integer :: current_layer,i
    real :: CanopyA,CanopyT, GCanopyA,Assim_Idir !initilization for the loop interation
    CanopyA=0.0
    CanopyT=0.0
    GCanopyA=0.0
    Assim_Idir=0.0 !initilization for the loop interation
!   if (.not.allocated(light_profile)) then 
!    allocate(light_profile(nlayers))
!   end if
    !solarR=solar_dir+solar_diff
    !lightME_tmp = lightME(lat, DOY,cosz)
    Idir1 = solar_dir!lightME_tmp%dir_frac * solarR
    Idiff = solar_diff!lightME_tmp%diff_frac * solarR
    cosTh = cosz!cosine_zenith_angle(lat,DOY,hr)

    light_profile = sunML(Idir1, Idiff, LAI, nlayers, cosTh, kd, chil, heightf)
        !print *, "Idir1, Idiff, LAI, nlayers, cosTh, kd, chil, heightf",Idir1,&
        !                        Idiff, LAI, nlayers, cosTh, kd, chil, heightf
       ! print *,"light_profile%dir_irr is",light_profile%dir_irr 
!      results from multilayer model 
    LAIc = LAI / nlayers
    
!      /* Next I need the RH and wind profile */
    relative_humidity_profile=RHprof(RH, nlayers)

    wind_speed_profile=WINDprof(WindSpeed, LAI, nlayers)

    leafN_profile=LNprof(leafN, LAI, nlayers, kpLN)
    
    do i=1,nlayers
        current_layer = nlayers - i+1
        leafN_lay = leafN_profile(current_layer)
        if (lnfun == 0) then
            vmax1 = Vmax
        else
            vmax1=nitroP%Vmaxb1*leafN_lay+nitroP%Vmaxb0
            if (vmax1<0.0) then
                vmax1=0.0
            end if
		    if (vmax1>Vmax) then 
                vmax1=Vmax
            end if
            Alpha=nitroP%alphab1*leafN_lay+nitroP%alphab0
            Rd=nitroP%Rdb1*leafN_lay+nitroP%Rdb0
        end if

        rhi = relative_humidity_profile(current_layer)
        layerWindSpeed = wind_speed_profile(current_layer)

        Idir2 = light_profile(current_layer)%dir_irr
        Itot = light_profile(current_layer)%tot_irr
        pLeafsun = light_profile(current_layer)%sunlit_frac
        CanHeight = light_profile(current_layer)%height
        Leafsun = LAIc * pLeafsun
        temp_photo_results = c4photoC(Idir2, Temp, rhi, vmax1, Alpha, Kparm, theta, beta, &
        	Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
        tmp5_ET = EvapoTrans2(Idir2, Itot, Temp, rhi, layerWindSpeed, LAIc, CanHeight, temp_photo_results%Gs, leafwidth, eteq)

        TempIdir = Temp + tmp5_ET%Deltat
        !if(isnan(TempIdir)) then
        !print *, "Here is CanAC,Temp too large, Deltat is",tmp5_ET%Deltat 
        !print *, "Idir1 is",Idir1
        !print *, "Idiff is",Idiff
        !print *, "solarR is",solarR
        !print *, "Gs is",temp_photo_results%Gs
        !print *, "Idir2 is",Idir2
        !print *, "StomataWS is",StomataWS
        !print *, "Itot is",Itot
        !print *, "rh is ",rhi,i
        !print *, "layerWindSpeed is",layerWindSpeed
        !print *, "LAI is ",LAI
        !print *, "Assim is ",temp_photo_results%Assim
        !print *, "Idir1, Idiff, LAI, nlayers, cosTh, kd, chil, heightf",Idir1,&
!Idiff, LAI, nlayers, cosTh, kd, chil, heightf
 !       endif
        temp_photo_results = c4photoC(Idir2, TempIdir, rhi, vmax1, Alpha, Kparm, &
        	theta, beta, Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
        AssIdir  = temp_photo_results%Assim
        GAssIdir = temp_photo_results%GrossAssim

        Idiff2 = light_profile(current_layer)%diff_irr
        pLeafshade = light_profile(current_layer)%shade_frac
        Leafshade = LAIc * pLeafshade

        temp_photo_results = c4photoC(Idiff2, Temp, rhi, vmax1, Alpha, Kparm, theta, &
        	beta, Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
        tmp6_ET = EvapoTrans2(Idiff2, Itot, Temp, rhi, layerWindSpeed, LAIc, CanHeight, temp_photo_results%Gs, leafwidth, eteq)
        TempIdiff = Temp + tmp6_ET%Deltat
        temp_photo_results = c4photoC(Idiff2, TempIdiff, rhi, vmax1, Alpha, Kparm, theta, &
        	beta, Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
        AssIdiff = temp_photo_results%Assim
        GAssIdiff = temp_photo_results%GrossAssim

        CanopyA = CanopyA + Leafsun * AssIdir + Leafshade * AssIdiff
        CanopyT = CanopyT + Leafsun * tmp5_ET%TransR + Leafshade * tmp6_ET%TransR
        GCanopyA = GCanopyA + Leafsun * GAssIdir + Leafshade * GAssIdiff
        Assim_Idir=Assim_Idir+AssIdir
    end do
!print *,"CanopyA,CanopyT is ",CanopyA,CanopyT 
! ## These are micromoles of CO2 per m2 per sec for Assimilation
! ## and mili mols of H2O per m2 per sec for Transpiration
! ## Need to convert to 
! ## 3600 converts seconds to hours
! ## 10^-6 converts micro mols to mols
! ## 30 converts mols of CO2 to grams
! ## (1/10^6) converts grams to Mg
! ## 10000 scales up to ha 
!      A similar conversion is made for water but
!        replacing 30 by 18 and mili mols are converted to
!        mols (instead of micro) 
!   deallocate(light_profile)
    CanACstr%Assim = cf * CanopyA
    CanACstr%Trans = cf2 * CanopyT 
    CanACstr%GrossAssim = cf * GCanopyA
    CanACstr%Idir=Idir1
    CanACstr%Idiff=Idiff
    CanACstr%Assim_Idir=Assim_Idir
end function CanAC

! /* New EvapoTrans function */
function EvapoTrans2(Rad,Iave,Airtemperature,RH,WindSpeed,LeafAreaIndex,&
        CanopyHeight,stomatacond,leafw,eteq) result(ETstr)
    real :: Rad,Iave,Airtemperature,RH,WindSpeed,LeafAreaIndex,&
        CanopyHeight,stomatacond,leafw
    real :: Tair,DdryA,LHV,SlopeFS,SWVP,SWVC,PsycParam,DeltaPVa,ActualVaporPressure,&
    totalradiation,Ja,Ja2,LayerWindSpeed,gvs,Deltat,ChangeInLeafTemp,OldDeltaT,rlc,ga,&
    PhiN2,TopValue,BottomValue,PhiN,TransR,EPen,EPries
    integer :: Counter,eteq
!     /* creating the structure to return */
    type(ETstructure)::ETstr
real::WindSpeedHeight,tau,LeafReflectance,SpecificHeat,StefanBoltzmann
!     //  kappa = 0.41 /* von Karmans ant */ unused
     WindSpeedHeight = 2.0 !/* This is the height at which the wind speed was measured */
!     //  dCoef = 0.77  unused
    tau = 0.2 !/* Leaf transmission coefficient */
    !//  ZetaCoef = 0.026 unused
    !//  ZetaMCoef = 0.13 unused
     LeafReflectance = 0.2 
     SpecificHeat = 1010.0 !/* J kg-1 K-1 */
     StefanBoltzmann = 5.67037e-8 !/* J m^-2 s^-1 K^-4 */

    !gvs /* Conductance to vapor from stomata same as stomatacond (input variable) */ 

    !rlc /* Long wave radiation for iterative calculation */
    !// WindSpeedTopCanopy = WindSpeed // unused.
    Tair = Airtemperature

    if (CanopyHeight < 0.1) then
        CanopyHeight = 0.1
    end if
! 
!     /* When the height at which wind was measured is lower than the canopy height */
!     /* There can be problems with the calculations */
!     /* This is a very crude way of solving this problem */
    if (CanopyHeight + 1 > WindSpeedHeight) then
        WindSpeedHeight = CanopyHeight + WindSpeedHeight
    end if

    DdryA = 1.295163636 + (-0.004258182) * Tair !/* Density of dry air, kg / m^3 */

!     /* In the original code in WIMOVAC this is used in J kg-1
!        but Thornley and Johnson use it as MJ kg-1  */
    LHV = 2.501 + (-0.002372727) *Tair !/* This should be MJ kg^-1*/
    LHV = LHV * 1e6 !/* Now it is converted to Joules kg^-1*/
    SlopeFS = TempToSFS(Tair) * 1e-3 !/* kg m^-3 K^-1 */
    SWVP = TempToSWVC(Tair) !/* this is hecto Pascals */
    !/* Convert to kg/m3 */
    SWVC = (DdryA * 0.622 * SWVP)/1013.25 !/* This last number is atmospheric pressure in hecto pascals */
   ! /* SWVC is saturated water vapor concentration (or density) in kg/m3 */

    PsycParam =(DdryA * SpecificHeat) / LHV !/* This is in kg m-3 K-1 */

    DeltaPVa = SWVC * (1.0 - RH) !/* kg/m3 */

    ActualVaporPressure = RH * SWVP !/* hecto Pascals */

!     /* SOLAR RADIATION COMPONENT*/
! 
!     /* First calculate the Radiation term */
!     /*' Convert light assuming 1 micromole PAR photons = 0.235 J */
!     /* The next step converts from PAR photons to Joules */
!     /* There are 2.35 x 10^5 Joules in a mol */
!     /* or 0.235 Joules in a micro mol */
!     /* A Watt is equivalent to J/s */
    totalradiation = Rad * 0.235 !/* This is essentially Watts m^-2 */
!     /* On a clear sky it may exceed 1000 in some parts of the world 
!        Thornley and Johnson pg 400 */
!     /* This values can not possibly be higher than 650 */
    if (totalradiation > 650.0) then
        error stop "total radiation too high"
    end if

!     /* Ja = (2 * totalradiation * ((1 - LeafReflectance - tau) / (1 - tau))) * LeafAreaIndex */
!     /* It seems that it is not correct to multiply by the leaf area index. The previous
!        version was used in WIMOVAC (check) */
    Ja = (2.0 * totalradiation * ((1.0 - LeafReflectance - tau) / (1.0 - tau)))

!     /* The value below is only for leaf temperature */
    Ja2 = (2.0 * Iave * 0.235 * ((1.0 - LeafReflectance - tau) / (1.0 - tau)))

!     /* AERODYNAMIC COMPONENT */
    if (WindSpeed < 0.5) then
        WindSpeed = 0.5
    end if

    LayerWindSpeed = WindSpeed

!     /* Rprintf("Gs !.3f \n", stomatacond) */
!     /* Leaf Conductance */
    gvs = stomatacond 
!     /* Convert from mmol H20/m2/s to m/s */
    gvs = gvs * (1.0/41000.0)
!     /* 1/41000 is the same as 24.39 * 1e-6 */
!     /* Thornley and Johnson use m s^-1 on page 418 */

!     /* prevent errors due to extremely low Layer conductance */
    if (gvs <= 0.001) then
        gvs = 0.001
    end if

!     /* This is the original from WIMOVAC*/
    Deltat = 0.01
    ChangeInLeafTemp = 10.0

    Counter = 0
    do while (ChangeInLeafTemp > 0.5 .and. Counter <= 10)
    
        OldDeltaT = Deltat

        rlc = 4.0 * StefanBoltzmann * (273.0 + Tair)**3.0 * Deltat  

!         /* rlc=net long wave radiation emittted per second =radiation emitted per second - radiation absorbed per second=sigma*(Tair+deltaT)^4-sigma*Tair^4 */
! 
!         /* Then you do a Taylor series about deltaT = 0 and keep only the zero and first order terms. */
! 
!         /* or rlc=sigma*Tair^4+deltaT*(4*sigma*Tair^3)-sigma*Tair^4=4*sigma*Tair^3*deltaT */
! 
!         /* where 4*sigma*Tair^3 is the derivative of sigma*(Tair+deltaT)^4 evaluated at deltaT=0, */

        ga = leafboundarylayer(LayerWindSpeed, leafw,& 
                Airtemperature, Deltat,&
                gvs, ActualVaporPressure)
!         /* This returns leaf-level boundary layer conductance */ 
!         /* In WIMOVAC this was added to the canopy conductance */
!         /* ga = (ga * gbcW)/(ga + gbcW)  */

        PhiN2 = (Ja2 - rlc)  !/* * LeafAreaIndex  */

        TopValue = PhiN2 * (1.0 / ga + 1.0 / gvs) - LHV * DeltaPVa
        BottomValue = LHV * (SlopeFS + PsycParam * (1.0 + ga / gvs))
        Deltat = TopValue / BottomValue !/* This equation is from Thornley and Johnson pg. 418 */
        if (Deltat > 10.0) then
            Deltat = 10.0
        end if
        if (Deltat < -10.0) then
            Deltat = -10.0
        end if
        ChangeInLeafTemp = abs(OldDeltaT - Deltat)

        Counter=Counter+1
    end do
     

!     /* Net radiation */
    PhiN = Ja - rlc

    if (PhiN < 0.0) then
        PhiN = 0.0
    end if

    TransR = (SlopeFS * PhiN + (LHV * PsycParam * ga * DeltaPVa)) / (LHV * (SlopeFS + PsycParam * (1.0 + ga / gvs)))

!     /* Penman will use the WIMOVAC conductance */
    EPen = (((SlopeFS * PhiN) + LHV * PsycParam * ga * DeltaPVa)) / (LHV * (SlopeFS + PsycParam))

    EPries = 1.26 * ((SlopeFS * PhiN) / (LHV * (SlopeFS + PsycParam)))

!     /* Choose equation to report */
    if (eteq == 1 ) then
        TransR = EPen
    end if
    if (eteq == 2) then 
        TransR = EPries
    end if

!     /* This values need to be converted from Kg/m2/s to
!        mmol H20 /m2/s according to S Humphries */
!     /* 1e3 - kgrams to grams  */
!     /* 1e3 - mols to mmols */
!     /* grams to mols - 18g in a mol */
!     /* Let us return the structure now */

    ETstr%TransR = TransR * 1e6 / 18.0 
    ETstr%EPenman = EPen * 1e6 / 18.0 
    ETstr%EPriestly = EPries * 1e6 / 18.0 
    ETstr%Deltat = Deltat
    ETstr%LayerCond = gvs * 41000.0   
end function EvapoTrans2

function leafboundarylayer( windspeed,  leafwidth,  AirTemp,&
         deltat,  stomcond,  vappress) result(gbv)
real :: windspeed,  leafwidth,  AirTemp,&
         deltat,  stomcond,  vappress
!     /* This is the leaf boundary layer computed using the approach in MLcan
!        which is based on (Nikolov, Massman, Schoettle),         !
!        Ecological Modelling, 80 (1995), 205-235 */
 real,parameter ::     Pa = 101325.0,  cf = 1.6361e-3
 real :: gbv,leaftemp,gsv,Tak,Tlk,ea,ws,lw,esTl,gbv_forced,gbv_free,eb,Tvdiff

     leaftemp = AirTemp + deltat
     gsv = stomcond !/* input is in m/s */
     Tak = AirTemp + 273.15 !/* Converts from C to K */
     Tlk = leaftemp + 273.15  !/* Converts from C to K */
     ea = vappress * 1e2 !/* From hPa to Pa */
     ws = windspeed !/* m s^-1 */
     lw = leafwidth !/* meters */

    esTl = TempToSWVC(leaftemp) * 100.0 !/* The function returns hPa, but need Pa */

!     /* Forced convection */ 
    gbv_forced = cf *  Tak**0.56 * ((Tak+120.0)*((ws/lw)/Pa))**0.5
    gbv_free = gbv_forced
    eb = (gsv * esTl + gbv_free * ea)/(gsv + gbv_free) !/* Eq 35 */
    Tvdiff = (Tlk / (1.0 - 0.378 * eb/Pa)) - (Tak / (1.0-0.378*ea/Pa)) !/* Eq 34*/

    if (Tvdiff < 0.0) then
        Tvdiff = -Tvdiff
    end if

    gbv_free = cf * Tlk**0.56 * ((Tlk+120)/Pa)**0.5 * (Tvdiff/lw)**0.25

    if (gbv_forced > gbv_free) then
        gbv = gbv_forced
    else
        gbv = gbv_free
    end if
end function leafboundarylayer
!     // gbh = 0.924 * gbv // set but not used

function TempToSFS(Temp) result(SlopeFS)
real::Temp, SlopeFS

    SlopeFS = 0.338376068 +  0.011435897 * Temp +  0.001111111 * Temp**2

end function TempToSFS

function TempToSWVC(Temp) result(SWVC)
real :: Temp,SWVC,a,b
    a = (18.678 - Temp/234.5) * Temp
    b = 257.14 + Temp
!     /* SWVC = (6.1121 * exp(a/b))/10; */
    SWVC = (6.1121 * exp(a/b))
end function TempToSWVC

function RHprof(RH, nlayers) result(relative_humidity_profile)
	real :: RH
	integer :: nlayers,i
	real,dimension(1:nlayers) ::relative_humidity_profile
	real :: kh,temp_rh
kh = 1.0 - RH
!     /* kh = 0.2; */
!     /*kh = log(1/RH);*/
    do i = 1, nlayers
        temp_rh = RH * exp(kh * (real(i)/nlayers))
!         // temp_rh = RH * exp(-kh * (j/nlayers));  // new simpler version from Joe Iverson*
        if (temp_rh > 1 ) then
            temp_rh = 0.99
        end if
        relative_humidity_profile(i) = temp_rh
    end do
!     /* It should return values in the 0-1 range */
end function RHprof

function WINDprof(WindSpeed, LAI, nlayers) result(wind_speed_profile)
	integer::nlayers,i
	real::LAI,WindSpeed
	real,dimension(1:nlayers) :: wind_speed_profile
	real::k,LI,CumLAI
 k = 0.7
 LI = LAI / nlayers
    do i=1,nlayers
        CumLAI = LI * i
        wind_speed_profile(i) = WindSpeed * exp(-k * (CumLAI-LI))
    end do
end function WINDprof

function LNprof(LeafN, LAI, nlayers, kpLN) result(leafN_profile)
    integer::nlayers,i
	real::LeafN, LAI, kpLN
	real,dimension(1:nlayers) :: leafN_profile
	real::LI,CumLAI
    LI = LAI / nlayers
    do i = 1, nlayers
        CumLAI = LI * i
        leafN_profile(i) = LeafN * exp(-kpLN * (CumLAI-LI))
    end do
end function LNprof

function cosine_zenith_angle(lat,doy,hod) result(cza)
	real::lat,cza,omega,delta,phi,tau,hod
	integer::doy
omega=360.0*(doy+10.0)/365.0*pi/180.0
delta=-23.5*pi/180.0*cos(omega)
phi=lat*pi/180.0
tau=(hod-12.0)*15.0*pi/180.0
cza=sin(delta)*sin(phi)+cos(delta)*cos(phi)*cos(tau)
end function cosine_zenith_angle

function lightME(lat,doy,cza) result(lightMEstr)
type(lightMEstructure)::lightMEstr
real::lat
integer::doy
real::atmo_trans,atmo_press,local_atmo,press_r,prop_scat,dir_frac,diff_frac,cza
atmo_trans=0.85
atmo_press=1e5
local_atmo=1e5
press_r=local_atmo/atmo_press
prop_scat=0.3
dir_frac=atmo_trans**(press_r/cza)
diff_frac=prop_scat*(1.0-dir_frac)*cza
if(cza<=0)then
dir_frac=0.0
diff_frac=1.0
end if 
lightMEstr%dir_frac=dir_frac/(dir_frac+diff_frac)
lightMEstr%diff_frac=diff_frac/(dir_frac+diff_frac)
end function lightME

function sunML(soldir,soldiff,&
    lai,nlayers,cosTheta,kd,chil,heightf) result(sunMLstr)
integer::nlayers,i
type(sunMLstructure),dimension(nlayers) :: sunMLstr
real::soldir,soldiff,lai,cosTheta,kd,chil,heightf
real::theta,k0,k1,k,laii,Ibeam,Isolar,cumlai,Iscat,Idiff,&
Ls,Ld,Fsun,Fshade,Iave
real,parameter::alphascatter = 0.8
theta=acos(min(max(cosTheta,-1.0),1.0))!Yufeng: make sure acos takes a value between -1 and 1
k0=sqrt(chil**2.0+(tan(theta))**2.0)
k1=chil+1.744*(chil+1.183)**(-0.733)
k=k0/k1!k Foliar absorption coefficent
if (k < 0.0)then
  k = -k
end if
!print *, "k0,k1,cosTheta,chil,theta are",k0,k1,cosTheta,chil,theta
laii=lai/nlayers
Ibeam=soldir*cosTheta
Isolar=Ibeam*k
do i=1,nlayers
    if (cosTheta<=1e-10) then
        Isolar=soldir/k1
        Idiff=soldiff*exp(-kd*cumlai)
        Fsun=0.0
        Fshade=0.0
        Iave=0.0
    else
        cumlai=laii*(i-1.0+0.5)
        Iscat=Ibeam*(exp(-k*sqrt(alphascatter)*cumlai)-exp(-k*cumlai))
        Idiff=soldiff*exp(-kd*cumlai)+Iscat !kd: Light extinction coefficient for diffuse light
        Ls=(1.0-exp(-k*laii))*exp(-k*cumlai)/k
        Ld=laii-Ls
        Fsun=Ls/(Ls+Ld)
        Fshade=Ld/(Ls+Ld)
        Iave=(Fsun*(Isolar+Idiff)+Fshade*Idiff)*(1.0-exp(-k*laii))/k    
    end if
    sunMLstr(i)%dir_irr=Isolar+Idiff
    sunMLstr(i)%diff_irr=Idiff
    sunMLstr(i)%tot_irr=Iave
    sunMLstr(i)%sunlit_frac=Fsun
    sunMLstr(i)%shade_frac=Fshade
    sunMLstr(i)%height=(lai-cumlai)/heightf
end do
end function sunML
! /* Function to simulate the multilayer behavior of soil water. In the
!    future this could be coupled with Campbell (BASIC) ideas to
!    esitmate water potential. */
function soilML( precipit,  transp,  cws,  soildepth,  depths,&
         fieldc,  wiltp,  phi1,  phi2, soTexS,  wsFun,&
         layers,  rootDB,  LAI,  k,  AirTemp,  IRad,  winds,&
         RelH,  hydrDist,  rfl,  rsec,  rsdf) result(soilMLstr)
real::precipit,  transp,  soildepth,&
         fieldc,  wiltp,  phi1,  phi2,&
         rootDB,  LAI,  k,  AirTemp,  IRad,  winds,&
         RelH,  rfl,  rsec,  rsdf
integer::layers, wsFun,hydrDist,i,j
real, parameter :: g = 9.8 !/* m / s-2  ##  http://en.wikipedia.org/wiki/Standard_gravity */
type(soilMLstructure)::soilMLstr
type(soTexStructure) :: soTexS
real::theta_s,rootDepth,waterIn,psim1,psim2,dPsim,K_psim,J_w,aw,diffw,rootATdepth,paw,awc,awc2,&
Sevap,Ctransp,EvapoTra,Newpawha,slp,intcpt,wsSpleaf,layerDepth,Nleach,phi10,theta,wsPhoto,wsPhotoCol,&
cw,runoff
real,dimension(layers)::root_distribution,  cws
real,dimension(layers+1)::depths
real:: drainage,oldWaterIn,oldEvapoTra,wsSpleafCol
!     /* Variables */
drainage = 0.0
oldWaterIn=0.0
oldEvapoTra = 0.0
wsSpleafCol=0.0
J_w=0.0
Nleach=0.0
wsPhotoCol=0.0
!     /* Here is a convention aw is available water in volume and awc
!        is available water content as a fraction of the soil section being investigated.
!        paw is plant available water aw - wiltp */
if (.not.allocated(soilMLstr%cws)) then
allocate(soilMLstr%cws(layers))
allocate(soilMLstr%hourlyWflux(layers))
allocate(soilMLstr%rootDist(layers))
end if
!     /* Specify the soil type */

    if (fieldc < 0.0) then
        fieldc = soTexS%fieldc
    end if
    if (wiltp < 0) then
        wiltp = soTexS%wiltp
    end if

    theta_s = soTexS%satur
!     /* rooting depth */
!     /* Crude empirical relationship between root biomass and rooting depth*/
    rootDepth = rootDB * rsdf
    if (rootDepth > soildepth) then
        rootDepth = soildepth
    end if

    root_distribution = rootDist(layers,rootDepth,depths,rfl)

     ! unit conversion for precip 
    waterIn = precipit * 1e-3 !/* convert precip in mm to m*/
    if(waterIn>0) then
    do i=1,layers
            layerDepth = depths(i+1) - depths(i)
           if(i == 1) then
 !             !/* I only add the water to the first layer */
              !/* This model does not really consider the infiltration rate and
!therefore runoff */
            cw = (cws(i) * layerDepth) + waterIn
           else
             cw = (cws(i) * layerDepth) + oldWaterIn
           endif
            cws(i) = cw / layerDepth 
     diffw = theta_s * layerDepth - cw
    if (diffw<0) then
        if (i==1) then
                runoff=-diffw  !Yufeng: surface runoff!
                cws(i)=theta_s
                oldWaterIn=0.0
        else
                cws(i)=theta_s ! Yufeng
                oldWaterIn=-diffw
        endif
    else
        oldWaterIn=0.0
    endif
    enddo
    endif
        print *,"before HydrDist,cws is ",cws
        do i=1,layers
        cws(i)=min(cws(i),theta_s)
        cws(i)=max(cws(i),1e-2)
        enddo
    !i=layers
    do i=1,layers
!         /* This supports unequal depths. */
            layerDepth = depths(i+1) - depths(i)

      !  print *,"cws(i+1),cws(i),layerDepth are ",cws(i+1),cws(i),layerDepth
        if (hydrDist > 0) then
!             /* For this section see Campbell and Norman "Environmental BioPhysics" Chapter 9*/
!             /* First compute the matric potential */
            psim1 = soTexS%air_entry * ((cws(i)/theta_s))**(-soTexS%b) !/* This is matric potential of current layer */
            if (i == layers) then  !last layer
               ! psim2=soTexS%air_entry !Assumes below the last layer there is saturation
                dPsim=0.0
            else
                psim2 = soTexS%air_entry * ((cws(i+1)/theta_s))**(-soTexS%b)!/* This is matric potential of next layer */
                dPsim = psim2 - psim1 !J/kg 
            end if
            K_psim = soTexS%Ks * (soTexS%air_entry/psim1)**(2.0+3.0/soTexS%b) !kg s /m3, This is hydraulic conductivity 
            J_w = K_psim * (dPsim/layerDepth) + g * K_psim !/*  Campbell, pg 129 do not ignore the graviational effect*/
!             /* This last result should be in kg/(m2 * s)*/
            J_w = J_w * 3600* 0.9882 * 1e-3 !  This is flow in m3 / (m^2 * hr). 
            !Yufeng : cws is in 1e-3 kg /m3 or  m3/m3
             if (i == layers) then
                cws(i) = cws(i) - J_w / layerDepth
                drainage = J_w
            else
                cws(i)=cws(i) - J_w/layerDepth
                cws(i+1)=cws(i+1) + J_w/layerDepth
            end if
        end if
	!if (i==1) &
      print *,"cws(i+1),cws(i),K_psim,dPsim,J_w,psim1,psim2,i are ",cws(i+1),cws(i),K_psim,dPsim,J_w,psim1,psim2,i

        if (cws(i) > theta_s) then
            cws(i) = theta_s 
        end if
!         /* if(cws[i+1] > fieldc) cws[i+1] = fieldc */
        if (cws(i) < wiltp) then
            cws(i) = wiltp 
        end if
!         /* if(cws[i+1] < wiltp) cws[i+1] = wiltp  */

        aw = cws(i) * layerDepth
!         /* Available water (for this layer) is the current water status times the layer depth */

!         /* Root Biomass */
        rootATdepth = rootDB * root_distribution(i)
        soilMLstr%rootDist(i) = rootATdepth
!         /* Plant available water is only between current water status and permanent wilting po */
!         /* Plant available water */
        paw = aw - wiltp * layerDepth
        if (paw < 0) then
            paw = 0.0 
        end if

        if (i == 1) then
!             /* Only the first layer is affected by soil evaporation */
            awc2 = aw / layerDepth
!             /* SoilEvapo function needs soil water content  */
            Sevap = SoilEvapo(LAI,k,AirTemp,IRad,awc2,fieldc,wiltp,winds,RelH,rsec)
!             /* I assume that crop transpiration is distributed simlarly to
!                root density.  In other words the crop takes up water proportionally
!                to the amount of root in each respective layer.*/
            Ctransp = transp*root_distribution(1)
            EvapoTra = Ctransp + Sevap
            Newpawha = (paw * 1e4) - EvapoTra / 0.9982 !/* See the watstr function for this last number 0.9882 */
!             /* The first term in the rhs (paw * 1e4) is the m3 of water available in this layer.
!                EvapoTra is the Mg H2O ha-1 of transpired and evaporated water. 1/0.9882 converts from Mg to m3 */
        else
            Ctransp = transp*root_distribution(i)
            EvapoTra = Ctransp
            Newpawha = (paw * 1e4) - (EvapoTra + oldEvapoTra)
        end if

        if (Newpawha < 0) then
!             /* If the Demand is not satisfied by this layer. This will be stored and added to subsequent layers*/
            oldEvapoTra = -Newpawha
            aw = wiltp * layerDepth 
        awc=wiltp
        !----Yufeng
        !Newpawha=0.0
        !------
        else
        oldEvapoTra = 0.0
        paw = Newpawha / 1e4
        awc = paw / layerDepth + wiltp  
        end if
    !            print *, "after ET calc, awc is ",awc
!         /* This might look like a weird place to populate the structure, but is more convenient*/
        soilMLstr%cws(i) = awc
        soilMLstr%hourlyWflux(i) =J_w
        if (wsFun == 0) then
            slp = 1.0/(fieldc - wiltp)
            intcpt = 1.0 - fieldc * slp
            wsPhoto = slp * awc + intcpt
        else if (wsFun == 1) then
            phi10 = (fieldc + wiltp)/2.0
            wsPhoto = 1.0/(1.0 + exp((phi10 - awc)/ phi1))
        else if (wsFun == 2) then
            slp = (1.0 - wiltp)/(fieldc - wiltp)
            intcpt = 1.0 - fieldc * slp
            theta = slp * awc + intcpt
            wsPhoto = (1.0 - exp(-2.5 * (theta - wiltp)/(1.0 - wiltp))) / (1.0 - exp(-2.5))
        else if (wsFun == 3) then
            wsPhoto = 1.0
        end if

        if (wsPhoto <= 0.0) then 
            wsPhoto = 1e-2 !/* This can be mathematically lower than zero in some cases but I should prevent that. */
            !Yufeng: I think wsPhoto is better not too close to zero. Affecting c4photoC !
        end if
       if (wsPhoto > 1.0) then  !Yufeng
           wsPhoto = 1.0
       end if
           wsPhotoCol = wsPhotoCol + &
              wsPhoto!*root_distribution(i)/sum(root_distribution)
  !print *,"awc,wsPhoto,root_distribution(i),i are ",awc,wsPhoto,root_distribution(i),i
       
       wsSpleaf = awc**phi2 * 1.0/(fieldc**phi2) 
      !print *,"i,awc,fieldc and phi2 are",i,awc,fieldc,phi2  
      if (wsFun == 3) then
            wsSpleaf = 1.0
        end if
       if (wsSpleaf > 1.0) then  !Yufeng
          wsSpleaf = 1.0
       end if
          wsSpleafCol = wsSpleafCol + &
                wsSpleaf!*root_distribution(i)/sum(root_distribution)
       !Yufeng weighted sum based on root_dist 
    end do

!         /* Need to convert to units used in the Parton et al 1988 paper. */
!         /* The data comes in mm/hr and it needs to be in cm/month */
        Nleach = drainage * 0.1 * (1.0/24.0*30.0) / (18.0 * (0.2 + 0.7 * soTexS%sand))
!     /* Apparently wsPhoto and wsSpleaf can be greater than 1 */
!    if (wsPhoto > 1) then
!        wsPhoto = 1.0
!    end if
!    if (wsSpleaf > 1) then
!        wsSpleaf = 1.0
!    end if

!     /* returning the structure */
    soilMLstr%rcoefPhoto = (wsPhotoCol/layers)
    soilMLstr%drainage = drainage
    soilMLstr%Nleach = Nleach
    soilMLstr%rcoefSpleaf = (wsSpleafCol/layers)
    soilMLstr%SoilEvapo = Sevap
end function soilML

function soilML_yh(cws, transp, soildepth, depths,fieldc,  wiltp,  phi1,  phi2,  wsFun,&
         layers,  rootDB,  LAI,  k,  AirTemp,  IRad,  winds,&
         RelH,  rfl,  rsec,  rsdf) result(soilMLstr)
real::   transp, soildepth, phi1,  phi2,&
         rootDB,  LAI,  k,  AirTemp,  IRad,  winds,&
         RelH,  rfl,  rsec,  rsdf
integer::layers, wsFun, i, nlayergtbedrock
type(soilMLstructure)::soilMLstr
real::theta_s,rootDepth,Sevap,Ctransp,EvapoTra,Newpawha,&
slp,intcpt,layerDepth,phi10,theta,wsPhoto,wsPhotoCol,&
aw,awc,paw,runoff
real,dimension(layers)::root_distribution,  cws,fieldc,wiltp
real,dimension(layers+1)::depths
real:: oldEvapoTra,wsSpleaf,wsSpleafCol
if (.not.allocated(soilMLstr%cws)) then
allocate(soilMLstr%cws(layers))
end if
wsPhotoCol=0.0
wsSpleafCol=0.0
root_distribution = 0.0
    rootDepth = rootDB * rsdf
    if (rootDepth > soildepth) then
        rootDepth = soildepth
    end if
nlayergtbedrock=0
do i=1,layers
	if (depths(i+1)<=soilDepth) then
	nlayergtbedrock=nlayergtbedrock+1
	endif
enddo
    root_distribution(1:nlayergtbedrock) = rootDist(nlayergtbedrock,rootDepth,&
                                            depths(1:nlayergtbedrock+1),rfl)
!print *, "root dist is ",root_distribution
	do i=1,nlayergtbedrock
	layerDepth=depths(i+1)-depths(i)
	aw = cws(i) * layerDepth
	paw = aw - wiltp(i) * layerDepth
	if (paw < 0) then
            paw = 0.0 
        end if
	if (i == 1) then
!             /* SoilEvapo function needs soil water content  */
            Sevap = SoilEvapo(LAI,k,AirTemp,IRad,cws(i),fieldc(i),wiltp(i),winds,RelH,rsec)
!             /* I assume that crop transpiration is distributed simlarly to
!                root density.  In other words the crop takes up water proportionally
!                to the amount of root in each respective layer.*/
	!	print *, "Sevap is",Sevap
            Ctransp = transp*root_distribution(i)
            EvapoTra = Ctransp + Sevap
            Newpawha = (paw * 1.e4) - EvapoTra / 0.9982 !/* See the watstr function for this last number 0.9882 */
!             /* The first term in the rhs (paw * 1e4) is the m3 of water available in this layer.
!                EvapoTra is the Mg H2O ha-1 of transpired and evaporated water. 1/0.9882 converts from Mg to m3 */
        else
            Ctransp = transp*root_distribution(i)
            EvapoTra = Ctransp
            Newpawha = (paw * 1.e4) - (EvapoTra / 0.9982 + oldEvapoTra)
        end if
        if (Newpawha < 0) then
!       /* If the Demand is not satisfied by this layer. This will be stored and added to subsequent layers*/
        oldEvapoTra = -Newpawha
        aw  = wiltp(i) * layerDepth 
        awc = wiltp(i)
        !----Yufeng
        !Newpawha=0.0
        !------
        else
        oldEvapoTra = 0.0
        paw = Newpawha / 1e4
        awc = paw / layerDepth + wiltp(i)  
        end if
	soilMLstr%cws(i) = awc
!	print *, "awc and i are",awc,i
	if (wsFun == 0) then
            slp = 1.0/(fieldc(i) - wiltp(i))
            intcpt = 1.0 - fieldc(i) * slp
            wsPhoto = slp * awc + intcpt
        else if (wsFun == 1) then
            phi10 = (fieldc(i) + wiltp(i))/2.0
            wsPhoto = 1.0/(1.0 + exp((phi10 - awc)/ phi1))
        else if (wsFun == 2) then
            slp = (1.0 - wiltp(i))/(fieldc(i) - wiltp(i))
            intcpt = 1.0 - fieldc(i) * slp
            theta = slp * awc + intcpt
            wsPhoto = (1.0 - exp(-2.5 * (theta - wiltp(i))/(1.0 - wiltp(i)))) / (1.0 - exp(-2.5))
        else if (wsFun == 3) then
            wsPhoto = 1.0
        end if

        if (wsPhoto <= 0.0) then 
            wsPhoto = 1e-2 !/* This can be mathematically lower than zero in some cases but I should prevent that. */
            !Yufeng: I think wsPhoto is better not too close to zero. Affecting c4photoC !
        end if
       if (wsPhoto > 1.0) then  !Yufeng
           wsPhoto = 1.0
       end if
           wsPhotoCol = wsPhotoCol + &
              wsPhoto!*root_distribution(i)/sum(root_distribution)
  !print *,"awc,wsPhoto,root_distribution(i),i are ",awc,wsPhoto,root_distribution(i),i
       
       wsSpleaf = awc**phi2 * 1.0/(fieldc(i)**phi2) 
      !print *,"i,awc,fieldc and phi2 are",i,awc,fieldc,phi2  
      if (wsFun == 3) then
            wsSpleaf = 1.0
        end if
       if (wsSpleaf > 1.0) then  !Yufeng
          wsSpleaf = 1.0
       end if
          wsSpleafCol = wsSpleafCol + &
                wsSpleaf!*root_distribution(i)/sum(root_distribution)
       !Yufeng weighted sum based on root_dist 
    enddo
	!     /* returning the structure */
    soilMLstr%rcoefPhoto = (wsPhotoCol)/nlayergtbedrock
    soilMLstr%rcoefSpleaf = (wsSpleafCol)/nlayergtbedrock
    soilMLstr%SoilEvapo = Sevap
end function soilML_yh

function rootDist(layer, rootDepth, depthsp, rfl) result(rootDistri)
    integer :: layer,i,j,k
    real::rootDepth, rfl,layerDepth,CumLayerDepth,a
    real,dimension(layer)::rootDistri, depthsp
integer :: CumRootDist
real::ca
CumRootDist=1
ca=0.0
CumLayerDepth=0.0
    do i=1,layer
            layerDepth = depthsp(i+1) - depthsp(i)

        CumLayerDepth = CumLayerDepth + layerDepth

        if (rootDepth > CumLayerDepth) then
            CumRootDist=CumRootDist+1
        end if
    end do

    do j=1,layer
        if (j <= CumRootDist) then
            a = poisson_density(j,CumRootDist*rfl)
            rootDistri(j) = a
            ca = ca + a
        else
            rootDistri(j) = 0.0
        end if
    end do

    do k=1,layer
        rootDistri(k) = rootDistri(k)/ ca
    end do
end function rootDist

function poisson_density(x, lambda) result(log_result)
real::log_result,lambda,factorial_x
integer::x
real,parameter::e=2.71828
factorial_x = sqrt(2.0 * pi * x) * (x / e)**x ! Stirling's approximation for n!.
log_result = -lambda + x * log(lambda) - log(factorial_x)
end function poisson_density

function SoilEvapo( LAI,  k,  AirTemp,  IRad,&
         awc,  fieldc,  wiltp,  winds,  RelH,  rsec ) result(Evaporation)
      !k = extinction coefficient 
integer,parameter::method = 1
real::LAI,  k,  AirTemp,  IRad,&
     awc,  fieldc,  wiltp,  winds,  RelH,  rsec
real::SoilArea,SoilTemp,rawc,Up,TotalRadiation,DdryA,LHV,SlopeFS,SWVC,PsycParam,&
DeltaPVa,BoundaryLayerThickness,DiffCoef,SoilBoundaryLayer,Ja,rlc,PhiN,Evaporation
real,parameter::cf2 = 3600 * 1e-3 * 18 * 1e-6 * 10000,SoilClodSize = 0.04,&
SoilReflectance = 0.2,SoilTransmission = 0.01,SpecificHeat = 1010.0,StefanBoltzman = 5.67e-8 

    ! For Transpiration 
    ! 3600 converts seconds to hours 
    ! 1e-3 converts mili mols to mols 
    ! 18 is the grams in one mol of H20 
    ! 1e-6 converts g to Mg 
    ! 10000 scales from meter squared to hectare 

    ! Let us assume a simple way of calculating the proportion of the
    !   soil with direct radiation 
    SoilArea = exp(-k * LAI)

    ! For now the temperature of the soil will be the same as the air.
   !    At a later time this can be made more accurate. I looked at the
    !   equations for this and the issue is that it is strongly dependent on
    !   depth. Since the soil model now has a single layer, this cannot be
    !   implemented correctly at the moment.  

    SoilTemp = AirTemp

    ! Let us use an idea of Campbell and Norman. Environmental
     !  Biophysics. 
    ! If relative available water content is 
    rawc = (awc - wiltp)/(fieldc - wiltp)
        if(rawc < 1e-7) rawc = 1e-7
        if(rawc > 1.0) rawc = 1.0

    ! Page 142 
    ! Maximum Dimensionless Uptake Rate 
    !Up = 1.0 - ((1.0 + 1.3 * rawc))**(-5.0)  
    ! This is a useful idea because dry soils evaporate little water when dry
     if(awc > fieldc) then
                Up =  1.0
     else
               ! Up = exp(rawc * 5.0) / exp(5.0)
                Up = 1.0 - (1.0 + 1.3 * rawc)**(-5.0)
     endif
       ! /* This value will be close to 1 with saturated soil and close to zero
       !         as the soil dries*/
       ! /* This is an empirical relationship that I made up */

    ! Total Radiation 
    !' Convert light assuming 1 micromole PAR photons = 0.235 J/s Watts
    ! At the moment soil evaporation is grossly overestimated. In WIMOVAC
!        the light reaching the last layer of leaves is used. Here instead
!        of calculating this again, I will for now assume a 10! as a rough
!        estimate. Note that I could maybe get this since layIdir and
!        layIDiff in sunML are external variables.  Rprintf("IRad
!        !.5f",layIdir[0],"\n") Update: 03-13-2009. I tried printing this
!        value but it is still too high and will likely overestimate soil
!        evaporation. However, this is still a work in progress.
       
    IRad = IRad * rsec ! Radiation soil evaporation coefficient  

    TotalRadiation = IRad * 0.235

    DdryA = 1.295163636 + (-0.004258182) *(AirTemp)
    LHV = (2.501 + (-0.002372727) * AirTemp )* 1e6 
    ! Here LHV is given in MJ kg-1 and this needs to be converted
!        to Joules kg-1  
    SlopeFS = TempToSFS(AirTemp) * 1e-3
    SWVC = TempToSWVC(AirTemp) * 1e-3
    
    PsycParam = (DdryA * SpecificHeat) / LHV
    DeltaPVa = SWVC * (1.0 - RelH) !Yufeng

    BoundaryLayerThickness = 4e-3 * sqrt(SoilClodSize / winds) 
    DiffCoef = 2.126e-5 * 1.48e-7 * SoilTemp
    SoilBoundaryLayer = DiffCoef / BoundaryLayerThickness

    Ja = 2.0 * TotalRadiation * ((1.0 - SoilReflectance - SoilTransmission) / (1.0 - SoilTransmission))

    rlc = 4.0 * StefanBoltzman * ((273.0 + SoilTemp))**3.0 * 0.005
    ! the last term should be the difference between air temperature and
!        soil. This is not actually calculated at the moment. Since this is
!        mostly relevant to the first soil layer where the temperatures are
!        similar. I will leave it like this for now. 

    PhiN = Ja - rlc ! Calculate the net radiation balance
    if (PhiN < 0) then
        PhiN = 1e-7
    end if

    ! Priestly-Taylor 
    if (method == 0) then
        Evaporation = 1.26 * (SlopeFS * PhiN) / (LHV * (SlopeFS + PsycParam))
     else 
        ! Penman-Monteith 
        Evaporation = (SlopeFS * PhiN + LHV * PsycParam * SoilBoundaryLayer * DeltaPVa) / (LHV * (SlopeFS + PsycParam))
    end if
    !  Report back the soil evaporation rate in Units mmoles/m2/s 
    !     Evaporation = Evaporation * 1000:   ' Convert Kg H20/m2/s to g H20/m2/s 
    !     Evaporation = Evaporation / 18:     ' Convert g H20/m2/s to moles H20/m2/s 
    !     Evaporation = Evaporation * 1000:   ' Convert moles H20/m2/s to mmoles H20/m2/s 

    !     If Evaporation <= 0 Then Evaporation = 0.00001: 
!            ' Prevent any odd looking values which might get through at very low light levels 

    Evaporation = Evaporation * 1.e6/18.
    ! Adding the area dependence and the effect of drying 
    ! Converting from m2 to ha (times 1e4) 
    ! Converting to hour 
    Evaporation = Evaporation * SoilArea * Up * cf2 
    if (Evaporation < 0 ) then
        Evaporation = 1e-6
        !else if (Evaporation > 100 .or. isnan(Evaporation)==.true.) then
        !print *, "Evaporation too large",Evaporation
    end if
end function SoilEvapo

! Function to select the correct dry biomass partitioning coefficients */
! It should take a vector of size 24 as an argument and return a structure with four numbers */
function sel_dbp_coef(coefs,  TherPrds,  TherTime) result(biomassPartistr)
!coefs: length 25
!TherPrds: length 6
real,dimension(25)::coefs
real,dimension(5):: TherPrds
real :: TherTime
type(biopartistructure)::biomassPartistr
    biomassPartistr%kLeaf = 0.0
    biomassPartistr%kStem = 0.0
    biomassPartistr%kRoot = 0.0
    biomassPartistr%kRhiz = 0.0
    biomassPartistr%kGrain = 0.0 ! kGrain is always zero except for the last thermal period */

    if (TherTime <= TherPrds(1)) then
    
        biomassPartistr%kStem = coefs(1)
        biomassPartistr%kLeaf = coefs(2)
        biomassPartistr%kRoot = coefs(3)
        biomassPartistr%kRhiz = coefs(4)

    else if (TherTime <= TherPrds(2)) then
    
        biomassPartistr%kStem = coefs(5)
        biomassPartistr%kLeaf = coefs(6)
        biomassPartistr%kRoot = coefs(7)
        biomassPartistr%kRhiz = coefs(8)

    else if (TherTime <= TherPrds(3)) then
    
        biomassPartistr%kStem = coefs(9)
        biomassPartistr%kLeaf = coefs(10)
        biomassPartistr%kRoot = coefs(11)
        biomassPartistr%kRhiz = coefs(12)

    else if (TherTime <= TherPrds(4)) then
    
        biomassPartistr%kStem = coefs(13)
        biomassPartistr%kLeaf = coefs(14)
        biomassPartistr%kRoot = coefs(15)
        biomassPartistr%kRhiz = coefs(16)

    else if (TherTime <= TherPrds(5)) then
        
        biomassPartistr%kStem = coefs(17)
        biomassPartistr%kLeaf = coefs(18)
        biomassPartistr%kRoot = coefs(19)
        biomassPartistr%kRhiz = coefs(20)

    else 

        biomassPartistr%kStem = coefs(21)
        biomassPartistr%kLeaf = coefs(22)
        biomassPartistr%kRoot = coefs(23)
        biomassPartistr%kRhiz = coefs(24)
        biomassPartistr%kGrain = coefs(25)

    end if
end function sel_dbp_coef

function leaf_n_limitation(kLn, leaf_n_0, current_leaf,current_stem) result(leafN)
real::kLn, leaf_n_0, current_leaf,current_stem,leaf_n
real::leafN
    leaf_n = leaf_n_0 * (current_leaf + current_stem)**(-kLn)
    if (leaf_n>leaf_n_0) then
        leafN=leaf_n_0
    else
        leafN=leaf_n
    end if
end function leaf_n_limitation

! /* Respiration. It is assumed that some of the energy produced by the
!    plant has to be used in basic tissue maintenance. This is handled
!    quite empirically by some relationships developed by McCree (1970)
!    and Penning de Vries (1972) */
real function resp( comp,  mrc,  temp)
real::comp,  mrc,  temp
    resp = comp *  (1.0 - (mrc * 2.0**(temp/10.0)))
    if (resp<0) then 
        resp = 0.
    end if
end function resp

function GroX(lat,doy,cosz,solar_dir,solar_diff,temp,rh,windspeed,precip,&
    kd,chil,leafwidth,et_equation,heightf,nlayers,initial_biomass,sencoefs,&
    timestep,iSp,SpD,dbpcoefs,thermalp,tbase,vmax1,alpha1,kparm,theta,&
    beta,Rd,Catm,b0,b1,soilcoefs,ileafn,kLN,vmaxb1,alphab1,mresp,wsFun,&
    ws,soilLayers,soilDepths,cws,&
    secs,kpLN,lnb0,lnb1,lnfun,upperT,lowerT,nitroP,StomataWS,&
    TTc,Grain,Leaf,Stem,Root,Rhizome,LeafN,Sp,vmax,alpha,&
    doy0,doy1,iWatCont,hour,fieldc,wiltp,bedrock) result(bio_out)
    real::lat,solar_dir,solar_diff,temp,rh,windspeed,precip,kd,chil,leafwidth,heightf,&
    iSp,SpD,tbase,vmax1,alpha1,kparm,theta,&
    beta,Rd,Catm,b0,b1,ileafn,kLN,vmaxb1,alphab1,&
    kpLN,lnb0,lnb1,upperT,lowerT,StomataWS,hour,cosz,bedrock
    real,dimension(4):: initial_biomass,sencoefs
    real,dimension(3):: secs
    real,dimension(25)::dbpcoefs
    real,dimension(5)::thermalp
    real,dimension(9)::soilcoefs
    real,dimension(2)::mresp
    integer::doy,et_equation,nlayers,wsFun,ws,lnfun,timestep,soilLayers,doy0,doy1
    integer::hydrDist,i,k,q,m,n,ri
    real,dimension(soilLayers)::water_status,root_distribution,psi,cws,iWatCont,fieldc,wiltp
    real,dimension(soilLayers+1)::soilDepths
type(nitroPstructure) :: nitroP
type(CanACstructure) :: Canopy
!type(soTexStructure) :: soTexS
type(biopartistructure)::dbpS
type(soilMLstructure)::soilMLS
type(wsStructure)::WaterS
type(biocro_outputs)::bio_out
real::Rhizome,Stem,Leaf,Root,Grain,LeafLitter,RootLitter,Remob,LeafN_0,LeafN,vmax,alpha,CanopyA,CanopyT,&
Sp,cwsVecSum,LeafPsim,LAI,rfl,rsec,rsdf,scsf,transpRes,leafPotTh,kLeaf,kStem,kRoot,&
kRhizome,kGrain,TTc,mrc1,mrc2,phi1,phi2,soilDepth,waterCont,seneLeaf,seneStem,&
seneRoot,seneRhizome,newStem,newLeaf,newRoot,newRhizome,newLeafcol,newRootcol,newRhizomecol,MinNitro,&
current_leaf,current_stem,LeafWS,newGrain,newStemcol,soilEvap,TotEvap,seneRate
cwsVecSum=0.0
if (.not.allocated(soilMLS%cws)) then
allocate(soilMLS%cws(soilLayers))
allocate(soilMLS%hourlyWflux(soilLayers))
allocate(soilMLS%rootDist(soilLayers))
end if

if (.not.allocated(bio_out%cws)) then
allocate(bio_out%cws(soilLayers))
end if
!if (.not.allocated(water_status)) then
!allocate(water_status(soilLayers))
!allocate(root_distribution(soilLayers))
!allocate(psi(soilLayers))
!allocate(cws(soilLayers))
!end if

!vectors: centcoefs
!    Variables needed for collecting litter, index +1 in MATLAB to Match C
    ! LeafLitter = centcoefs(21)
    ! StemLitter = centcoefs(22)
    ! RootLitter = centcoefs(23)
    ! RhizomeLitter = centcoefs(24)
    LeafN_0 = ileafn
    !Yufeng: 0.0005 is an estimate against Miguez et al., 2012, hourly sene rate 
    seneRate = 0.0005 
!     Maintenance respiration
    mrc1 = mresp(1)
    mrc2 = mresp(2)

    !FieldC = soilcoefs(1)
    !WiltP = soilcoefs(2)
    phi1 = soilcoefs(3)
    phi2 = soilcoefs(4)
    soilDepth = soilcoefs(5)
    soilDepth = min(soilDepth,bedrock)
    !waterCont = soilcoefs(6)
    
    seneLeaf = sencoefs(1)
    seneStem = sencoefs(2)
    seneRoot = sencoefs(3)
    seneRhizome = sencoefs(4)
    
    !soTexS = soilTchoose(soilType)
    
! Some soil related empirical coefficients
    rfl = secs(1)  ! root factor lambda 
    rsec = secs(2) ! radiation soil evaporation coefficient 
    rsdf = secs(3) ! root soil depth factor 
    scsf = soilcoefs(7) ! stomatal conductance sensitivity factor */ /* Rprintf("scsf !.2f",scsf) */
    transpRes = soilcoefs(8) ! Resistance to transpiration from soil to leaf 
    leafPotTh = soilcoefs(9) ! Leaf water potential threshold   
    
 !    water_status=zeros(soilLayers*vecsize)
 !    root_distribution=zeros(soilLayers*vecsize)
 ! psi=zeros(soilLayers*vecsize)
    
!Yufeng : make sure the condition for initialization is correct!
    if(doy==doy0 .and. (abs(hour-NINT(hour)) < 0.05)) then
    TTc=0.0
    Grain=0.0
    Rhizome = initial_biomass(1)
    Stem = initial_biomass(2)
    Leaf = initial_biomass(3)
    Root = initial_biomass(4)
    LeafN = ileafn ! Need to set it because it is used by CanA before it is computed
    Sp=iSp
    vmax = vmax1
    alpha = alpha1
    StomataWS=1.0
    !cws=iWatCont
    ! k=1
    ! q=1
    ! m=1
    ! n=1
    ! ri=1
    !else
    ! The ones to be iterated
    !TTc=TTc_ij
    !Grain=Grain_ij
    !Leaf=Leaf_ij
    !Stem=Stem_ij
    !Root=Root_ij
    !Rhizome=Rhizome_ij
    !LeafN=LeafN_ij
    !Sp=Sp_ij
    !vmax=vmax_ij
    !alpha=alpha_ij
    !elseif (doy>doy0) then
    endif
     LAI = min(10.0,Leaf * Sp)    

 !   for i = 1:vecsize
       !First calculate the elapsed Thermal Time
        if (temp > tbase) then
            TTc = TTc + (temp-tbase) / (24./timestep) 
        end if
        ! Do the magic! Calculate growth

        Canopy = CanAC(LAI, doy, cosz,&
                solar_dir,solar_diff, temp, rh, windspeed,&
                lat, nlayers, vmax, alpha, kparm, beta,&
                Rd, Catm, b0, b1, theta, kd, chil,&
                heightf, LeafN, kpLN, lnb0, lnb1, lnfun, upperT, lowerT,&
                nitroP, leafwidth, et_equation, StomataWS, ws)

        CanopyA = Canopy%Assim * timestep
        CanopyT = Canopy%Trans * timestep
      !print *,"CanopyA is",CanopyA
	!if(isnan(CanopyA)) then
        ! print *,"CanopyA too large and the CanopyT is",CanopyA,CanopyT
        ! print *,"LAI is",LAI
        ! print *,"rh is",rh
        ! print *,"temp is",temp
        ! print *,"date is",doy,hour
        ! print *,"StomataWS is",StomataWS
        ! print *,"solar is",solar_dir,solar_diff
        ! print *,"CanACstr%Idir,CanACstr%Idiff,CanACstr%Assim_Idir",Canopy%Idir,Canopy%Idiff,Canopy%Assim_Idir
	! print *,"cosz is",cosz
        !endif
        ! Inserting the multilayer model
!print *,"Yufeng debug 1,CanopyA and T are",CanopyA,CanopyT 
             soilMLS = soilML_yh(cws, CanopyT,soilDepth,soilDepths, fieldc, wiltp, phi1, phi2, wsFun,&
                              soilLayers, Root, LAI, 0.68, temp, solar_dir+solar_diff, windspeed,&
                              rh, rfl, rsec, rsdf)
             cws = soilMLS%cws
             StomataWS = soilMLS%rcoefPhoto
             LeafWS = soilMLS%rcoefSpleaf
             soilEvap =soilMLS%SoilEvapo
!print *,"StomataWS,LeafWS,solar_dir,solar_diff are ",StomataWS,LeafWS,solar_dir,solar_diff
!         An alternative way of computing water stress is by doing the leaf
!          water potential. This is done if the wsFun is equal to 4
!print *,"Yufeng debug 2"
        if (wsFun == 4) then
!             /* Calculating the leaf water potential */
!             /* From Campbell E = (Psim_s - Psim_l)/R or
!              * evaporation is equal to the soil water potential
!              * minus the leaf water potential divided by the resistance.
!              * This can be rearranged to Psim_l = Psim_s - E x R   */
!             /* It is assumed that total resistance is 5e6 m^4 s^-1
!              * kg^-1 
!              * Transpiration is in Mg ha-2 hour-1
!              * Multiply by 1e3 to go from Mg to kg
!              * Multiply by 1e-4 to go from ha to m^2 
!              * This needs to go from hours to seconds that's
!              * why the conversion factor is (1/3600).*/
            LeafPsim = WaterS%psim - (CanopyT * 1e3 * 1e-4 * 1.0/3600.0) * transpRes
! 
!             /* From WIMOVAVC the proposed equation to simulate the effect of water
!              * stress on stomatal conductance */
            if (LeafPsim < leafPotTh) then
!                 /* StomataWS = 1 - ((LeafPsim - leafPotTh)/1000 *
!                  * scsf) In WIMOVAC this equation is used but
!                  * the absolute values are taken from the
!                  * potentials. Since they both should be
!                  * negative and leafPotTh is greater than
!                  * LeafPsim this can be rearranged to*/ 
                StomataWS = 1 - ((leafPotTh - LeafPsim)/1000.0 * scsf)
!                 /* StomataWS = 1 */
                if (StomataWS < 0.1) then
                    StomataWS = 0.1
                end if
            else 
                StomataWS = 1.0
            end if
        else 
            LeafPsim = 0.0
        end if

!         /* Picking the dry biomass partitioning coefficients */
        dbpS = sel_dbp_coef(dbpcoefs, thermalp, TTc)

        kLeaf = dbpS%kLeaf
        kStem = dbpS%kStem
        kRoot = dbpS%kRoot
        kGrain = dbpS%kGrain
        kRhizome = dbpS%kRhiz

!         Here I can insert the code for Nitrogen limitations on photosynthesis
!            parameters. This is taken From Harley et al. (1992) Modelling cotton under
!            elevated CO2. PCE. This is modeled as a simple linear relationship between
!            leaf nitrogen and vmax and alpha. Leaf Nitrogen should be modulated by N
!            availability and possibly by the Thermal time accumulated.

        current_leaf=Leaf
        current_stem=Stem
        LeafN = leaf_n_limitation(kLN, LeafN_0, current_leaf,current_stem)

        vmax = (LeafN_0 - LeafN) * vmaxb1 + vmax1
        alpha = (LeafN_0 - LeafN) * alphab1 + alpha1

!         The crop demand for nitrogen is the leaf concentration times the amount of biomass.
!            This modifies the amount of N available in the soil. 
!            MinNitro is the available amount of N (kg/m2). 
!            The demand is in Mg/ha. I need a conversion factor of 
!            multiply by 1000, divide by 10000.

        MinNitro = MinNitro - LeafN * (Stem + Leaf) * 1e-1
        if (MinNitro < 0) then 
            MinNitro = 1e-3
        end if

        if (kLeaf > 0) then
            newLeaf = CanopyA * kLeaf * LeafWS
            newLeaf = resp(newLeaf, mrc1, temp)
            newLeafcol = newLeaf
        else
            newLeaf = Leaf * kLeaf
            Rhizome = Rhizome + kRhizome * (-newLeaf) * 0.9 !0.9 is the efficiency of retranslocation
            Stem = Stem + kStem * (-newLeaf) * 0.9
            Root = Root + kRoot * (-newLeaf) * 0.9
            Grain = Grain + kGrain * (-newLeaf) * 0.9
        end if

        if (TTc < seneLeaf) then!(doy <= doy0+0.6*(doy1-doy0)) then 
            Leaf = Leaf + newLeaf
        else 
	!	if ((nint(hour)==24) .and. (hour-nint(hour) < 0.05)) then
		    Leaf = (Leaf + newLeaf)- seneRate * Leaf
                    Remob = Leaf * seneRate*0.6
                    Rhizome = Rhizome + kRhizome * Remob
                    Stem = Stem + kStem * Remob
                    Root = Root + kRoot * Remob
                    Grain = Grain + kGrain * Remob
	!	endif
        !     Leaf = Leaf + newLeaf - newLeafcol(k) ! This means that the new value of leaf is
        !                                 !   the previous value plus the newLeaf
        !                                 !   (Senescence might start when there is
        !                                 !   still leaf being produced) minus the leaf
        !                                 !  produced at the corresponding k.
        !     Remob = newLeafcol(k) * 0.6
        !     LeafLitter = LeafLitter+ newLeafcol(k) * 0.4 !Collecting the leaf litter  
        !     Rhizome = Rhizome + kRhizome * Remob
        !     Stem = Stem + kStem * Remob
        !     Root = Root + kRoot * Remob
        !     Grain = Grain + kGrain * Remob
        !     k=k+1
        end if

!         The specific leaf area declines with the growing season at least in
!            Miscanthus.  See Danalatos, Nalianis and Kyritsis "Growth and Biomass
!            Productivity of Miscanthus sinensis "Giganteus" under optimum cultural
!            management in north-eastern greece

        if (mod(NINT(hour),24) == 0) then
            Sp = iSp - (doy - doy0) * SpD
        end if

        LAI = Leaf * Sp

      !New Stem
        if (kStem >= 0) then
            newStem = CanopyA * kStem
            newStem = resp(newStem, mrc1, temp)
            newStemcol = newStem
        else 
            error stop "kStem should be positive"
        end if

        if (TTc < seneStem) then !(doy <= doy0+0.6*(doy1-doy0)) then 
            Stem = Stem + newStem
	else
	!	if ((nint(hour)==24) .and. (hour-nint(hour) < 0.05)) then
		    Stem = (Stem + newStem) - seneRate * Stem
	!	endif
        end if
        
        if (kRoot > 0) then 
            newRoot = CanopyA * kRoot
            newRoot = resp(newRoot, mrc2, temp)
            newRootcol = newRoot
        else
            newRoot = Root * kRoot
            Rhizome = Rhizome + kRhizome * (-newRoot) * 0.9
            Stem = Stem + kStem * (-newRoot) * 0.9
            Leaf = Leaf + kLeaf * (-newRoot) * 0.9
            Grain = Grain + kGrain * (-newRoot) * 0.9
        end if

        if  (TTc < seneRoot) then 
            Root = Root + newRoot
        else
        !    if ((nint(hour)==24) .and. (hour-nint(hour) < 0.05)) then
               Root=(Root+newRoot)- seneRate * Root
        !    endif
        end if

        if (kRhizome > 0) then
            newRhizome = CanopyA * kRhizome
            newRhizome = resp(newRhizome, mrc2, temp)
            newRhizomecol = newRhizome
        else
            if (Rhizome < 0) then
                Rhizome = 1e-4
                print *, "warning:Rhizome became negative"
            end if

            newRhizome = Rhizome * kRhizome
            Root = Root + kRoot * (-newRhizome)
            Stem = Stem + kStem * (-newRhizome)
            Leaf = Leaf + kLeaf * (-newRhizome)
            Grain = Grain + kGrain * (-newRhizome)
        end if

        if  (TTc < seneRhizome) then 
            Rhizome = Rhizome + newRhizome
        else
          !  if ((nint(hour)==24) .and. (hour-nint(hour) < 0.05)) then
               Rhizome=(Rhizome+newRhizome) - seneRate * Rhizome
         !   endif
        end if
        
        if (kGrain < 1e-10 .or. TTc < thermalp(5)) then
            newGrain = 0.0
            Grain = Grain + newGrain
        else
            newGrain = CanopyA * kGrain
            Grain = Grain + newGrain
        end if
    
    bio_out%TTc=TTc
    bio_out%Grain=Grain
    bio_out%Leaf=Leaf
    bio_out%Stem=Stem
    bio_out%Root=Root
    bio_out%Rhizome=Rhizome
    bio_out%LeafN=LeafN
    bio_out%Sp=Sp
    bio_out%vmax=vmax
    bio_out%alpha=alpha
    bio_out%stomataWS=StomataWS
    bio_out%cws=cws         !m3/m3
    bio_out%biolai=LAI
    bio_out%bioev=soilEvap  !Mg ha-1 hr-1
    bio_out%biotr=CanopyT   !Mg ha-1 hr-1
    bio_out%CA = CanopyA
end function GroX

function stp_c2b(cwrf_soiltype) result(biocro_soiltype)
integer :: cwrf_soiltype,biocro_soiltype
if(cwrf_soiltype>12) then
biocro_soiltype=6
else
select case(cwrf_soiltype)
case (1)
biocro_soiltype=0
case (2)
biocro_soiltype=1
case (3)
biocro_soiltype=2
case (4)
biocro_soiltype=4
case (6)
biocro_soiltype=3
case (7)
biocro_soiltype=5
case (9)
biocro_soiltype=6
case (8)
biocro_soiltype=7
case (10)
biocro_soiltype=8
case (11)
biocro_soiltype=9
case (12)
biocro_soiltype=10
end select
end if
end function stp_c2b

end module BiocroX
