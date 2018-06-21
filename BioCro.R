# Please keep BioCro.R updated referring to its Fortran version
#
# log: May 24 2018 (EST), construct all functions that
# are currently used in CWRF-BioCro
# 

#pi=3.1415926

AbiotEff <- function(smoist,stemp){
    if (stemp < 35.0){ 
    TempEff = 1.0087 / (1.0 + (46.2 * exp(-0.1899 * stemp)))
    } else {TempEff = -0.0826 * stemp + 3.84
    }
    MoisEff = 1.0267 / (1.0 + 14.7 * exp(-6.5 * smoist))
    Abiot=TempEff * MoisEff
    return(Abiot)
} 

ballBerry <- function(Amu, Cappm, Temp, RelH, beta0, beta1){
gbw = 1.2 # According to Collatz et al. (1992) pg. 526
ptotPa = 101325.0 # Atmospheric pressure 
leafTk = Temp + 273.15
pwi = fnpsvp(leafTk)
pwaPa = RelH * pwi
Camf = Cappm * 1e-6
assimn = Amu * 1e-6
wa  = pwaPa / ptotPa
wi  = pwi / ptotPa

if (assimn < 0.0){ 
	gswmol = beta0 
}else{
		Cs  = Camf - (1.4/gbw)*assimn
	if (Cs < 0.0){
		print(c("Camf,assimn are",Camf,assimn))
                stop("Cs is less than 0")
        }
		acs = assimn/Cs

	if (acs < 1e-6){
            acs = 1e-6		
        }
	aaa = beta1 * acs
	bbb = beta0 + gbw - (beta1 * acs)
	ccc = -(wa / wi) * gbw - beta0

	ddd = bbb * bbb - 4.0*aaa*ccc

	hs  = (-bbb + sqrt(ddd)) / (2.0* aaa)

	gswmol = beta1 * hs * acs + beta0
}
	gsmol = gswmol * 1000.0 # converting to mmol 

	if (gsmol <= 0.0){
        gsmol = 1e-2
        }
return(gsmol)
}

fnpsvp <- function(Tkelvin){
	tmp = Tkelvin - 273.15
	u = (18.678 - tmp/234.5)*tmp
	v = 257.14 + tmp
	esat = (6.1121 * exp(u/v))/10.0
return(esat)
}

c4photoC <- function(Qp, Tl, RH, vmax, alpha,kparm, theta, beta,Rd, bb0, bb1, StomaWS, Ca, ws,upperT,lowerT){
AP = 101325.0 #Atmospheric pressure According to wikipedia (Pa)*/
P = AP / 1e3 # kPa */
Q10 = 2.0  # Q10 increase in a reaction by 10 C temp */
	# Defining biochemical variables */
OldAssim = 0.0
Tol = 0.01

Csurface = (Ca * 1e-6) * AP 
  
InterCellularCO2 = Csurface * 0.4 # Initial guestimate */

KQ10 =  Q10**((Tl - 25.0) / 10.0)

KT = kparm * KQ10
 
# This is the code implementing temperature limitations
Vtn = vmax * 2.0**((Tl-25.0)/10.0)
Vtd = ( 1.0 + exp(0.3 * (lowerT-Tl)) ) * (1.0 + exp( 0.3*(Tl-upperT) ))
VT  = Vtn / Vtd

# Second chunk of code see Collatz (1992) */
Rtn = Rd * 2.0**((Tl-25.0)/10.0)
Rtd =  1.0 + exp( 1.3 * (Tl-55.0) ) 
RT = Rtn / Rtd  

# Third chunk of code again see Collatz (1992) */
b0 = VT * alpha  * Qp 
b1 = VT + alpha  * Qp 
b2 = theta 

# Calculate the 2 roots */
M1 = (b1 + sqrt(b1*b1 - (4.0 * b0 * b2)))/(2.0*b2) 
M2 = (b1 - sqrt(b1*b1 - (4.0 * b0 * b2)))/(2.0*b2) 

# This piece of code selects the smaller root */
M=min(M1,M2)
# Here the iterations will start */
iterCounter = 0
#if(isnan(M)) print *, "M too large and T and Q are",M,Tl,Qp
while (iterCounter < 50){
kT_IC_P = kT * (InterCellularCO2 / P*1000.0)
Quada = M * kT_IC_P
Quadb = M + kT_IC_P
Quadc = beta 

a21 = (Quadb - sqrt(Quadb*Quadb - (4.0 * Quada * Quadc))) / (2.0 * Quadc)
a22 = (Quadb + sqrt(Quadb*Quadb - (4.0 * Quada * Quadc))) / (2.0 * Quadc)
a2  =  min(a21,a22)
Assim = a2 - RT #Yufeng

	if (ws == 0) { 
            Assim = Assim * StomaWS 
        }
 		# milimole per meter square per second*/
	csurfaceppm = Csurface * 10.0 
 		# Need to create the Ball-Berry function */
	Gs =  ballBerry(Assim,csurfaceppm, Tl, RH, bb0, bb1) 
                #Gs is [mmol m-2 s-1]
        #Stomata conductance.
	if (ws == 1) { 
            Gs = Gs * StomaWS 
        }
        g = (Csurface - (Assim * 1.6 * AP) / (Gs * 1e3))-InterCellularCO2 #pa
        if (abs(g) < 1e-6) { 
            diff = 0.
            break 
        }else{
        h = g
        kT_IC_P = kT * ((InterCellularCO2 + h) / P * 1000.)
		Quada = M * kT_IC_P !M:umol m-2 s-1
		Quadb = M + kT_IC_P
		Quadc = beta 
		a21 = (Quadb - sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a22 = (Quadb + sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a2  =  min(a21,a22)
	Assim1 = a2 - RT
        Gs1 =  ballBerry(Assim1,csurfaceppm, Tl, RH, bb0, bb1) 
        if (ws == 1) {
            Gs1 = Gs1 * StomaWS 
        }
        g1=(Csurface - (Assim1 * 1.6 * AP) / (Gs1 * 1e3))-InterCellularCO2-h
        dg=(g1-g)/h # dg: derivative of g(x) at Ci
        if (abs(dg) < 1e-6) {
           diff = 0.
           break 
        }
        Ci2 = InterCellularCO2 - g / dg
        diff = Ci2 - InterCellularCO2
	if (abs(diff) < Tol) {
		break
        }else{
            InterCellularCO2 = Ci2
        }
	iterCounter=iterCounter+1
        } #if else g==0
    } #while

        kT_IC_P = kT * (InterCellularCO2 / P*1000.) #kT: mol m-2 s-1.1e6 -> umol
	Quada = M * kT_IC_P #M:umol m-2 s-1
	Quadb = M + kT_IC_P
	Quadc = beta 
	a21 = (Quadb - sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a22 = (Quadb + sqrt(Quadb*Quadb - (4. * Quada * Quadc))) / (2. * Quadc)
        a2  =  min(a21,a22)
	Assim = a2 - RT
        if (ws == 0){ 
            Assim = Assim * StomaWS 
        }
        Gs =  ballBerry(Assim,csurfaceppm, Tl, RH, bb0, bb1)
	if (ws == 1) {
            Gs = Gs * StomaWS 
        }
        InterCellularCO2 = Csurface - (Assim * 1.6 * AP) / (Gs * 1e3)

	miC = (InterCellularCO2 / AP) * 1e6
 
    if (Gs > 600.0){ 
	    Gs = 600.0
    }
    GrossAssim=Assim+RT
    c4str=list(Assim=Assim,Gs=Gs,Ci=miC,GrossAssim=GrossAssim)
return(c4str)
}

CanAC <- function(LAI,DOY,cosz,solar_dir,solar_diff,Temp,RH,WindSpeed,lat,nlayers,Vmax,Alpha,Kparm,beta,Rd,Catm,b0,b1,theta,kd,chil,heightf,leafN,kpLN,lnb0,lnb1,lnfun,upperT,lowerT,nitroP,leafwidth,eteq,StomataWS,ws){
cf = 3600 * 1e-6 * 44 * 1e-6 * 10000 
cf2 = 3600 * 1e-3 * 18 * 1e-6 * 10000 
CanopyA=0.0
CanopyT=0.0
GCanopyA=0.0
Assim_Idir=0.0 #initilization for the loop interation
    #solarR=solar_dir+solar_diff
    #lightME_tmp = lightME(lat, DOY,cosz)
Idir1 = solar_dir#lightME_tmp%dir_frac * solarR
Idiff = solar_diff#lightME_tmp%diff_frac * solarR
cosTh = cosz#cosine_zenith_angle(lat,DOY,hr)

light_profile = sunML(Idir1, Idiff, LAI, nlayers, cosTh, kd, chil, heightf)
LAIc = LAI / nlayers
    
relative_humidity_profile=RHprof(RH, nlayers)

wind_speed_profile=WINDprof(WindSpeed, LAI, nlayers)

leafN_profile=LNprof(leafN, LAI, nlayers, kpLN)
    
for (i in 1:nlayers){
    current_layer = nlayers - i+1
    leafN_lay = leafN_profile(current_layer)
    if (lnfun == 0){ 
            vmax1 = Vmax
    }else{
            vmax1=nitroP$Vmaxb1*leafN_lay+nitroP$Vmaxb0
         if (vmax1<0.0){
             vmax1=0.0
         }
	 if (vmax1>Vmax){ 
             vmax1=Vmax
         }
         Alpha=nitroP$alphab1*leafN_lay+nitroP$alphab0
         Rd=nitroP$Rdb1*leafN_lay+nitroP$Rdb0
    }

    rhi = relative_humidity_profile(current_layer)
    layerWindSpeed = wind_speed_profile(current_layer)

    Idir2 = light_profile(current_layer)%dir_irr
    Itot = light_profile(current_layer)%tot_irr
    pLeafsun = light_profile(current_layer)%sunlit_frac
    CanHeight = light_profile(current_layer)%height
    Leafsun = LAIc * pLeafsun
    temp_photo_results = c4photoC(Idir2, Temp, rhi, vmax1, Alpha, Kparm, theta, beta,Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
    tmp5_ET = EvapoTrans2(Idir2, Itot, Temp, rhi, layerWindSpeed, LAIc, CanHeight, temp_photo_results$Gs, leafwidth, eteq)

    TempIdir = Temp + tmp5_ET$Deltat
    temp_photo_results = c4photoC(Idir2, TempIdir, rhi, vmax1, Alpha, Kparm, theta, beta, Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
    AssIdir  = temp_photo_results$Assim
    GAssIdir = temp_photo_results$GrossAssim

    Idiff2 = light_profile(current_layer)%diff_irr
    pLeafshade = light_profile(current_layer)%shade_frac
    Leafshade = LAIc * pLeafshade

    temp_photo_results = c4photoC(Idiff2, Temp, rhi, vmax1, Alpha, Kparm, theta, beta, Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
    tmp6_ET = EvapoTrans2(Idiff2, Itot, Temp, rhi, layerWindSpeed, LAIc, CanHeight, temp_photo_results$Gs, leafwidth, eteq)
    TempIdiff = Temp + tmp6_ET$Deltat
    temp_photo_results = c4photoC(Idiff2, TempIdiff, rhi, vmax1, Alpha, Kparm, theta, beta, Rd, b0, b1, StomataWS, Catm, ws, upperT, lowerT)
    AssIdiff = temp_photo_results$Assim
    GAssIdiff = temp_photo_results$GrossAssim

    CanopyA = CanopyA + Leafsun * AssIdir + Leafshade * AssIdiff
    CanopyT = CanopyT + Leafsun * tmp5_ET$TransR + Leafshade * tmp6_ET$TransR
    GCanopyA = GCanopyA + Leafsun * GAssIdir + Leafshade * GAssIdiff
    Assim_Idir=Assim_Idir+AssIdir
} 
    CanACstr=list(Assim=cf*CanopyA,Trans=cf2*CanopyT,GrossAssim=cf*GCanopyA,Idir=Idir1,Idiff=Idiff,Assim_Idir=Assim_Idir)
return(CanACstr)
}

EvapoTrans2 <- function(Rad,Iave,Airtemperature,RH,WindSpeed,LeafAreaIndex,CanopyHeight,stomatacond,leafw,eteq){
    WindSpeedHeight = 2.0 #/* This is the height at which the wind speed was measured */
    tau = 0.2 #/* Leaf transmission coefficient */
    LeafReflectance = 0.2 
    SpecificHeat = 1010.0 #/* J kg-1 K-1 */
    StefanBoltzmann = 5.67037e-8 #/* J m^-2 s^-1 K^-4 */
    Tair = Airtemperature

    if (CanopyHeight < 0.1){
        CanopyHeight = 0.1
    }
    if (CanopyHeight + 1 > WindSpeedHeight){
        WindSpeedHeight = CanopyHeight + WindSpeedHeight
    }

    DdryA = 1.295163636 + (-0.004258182) * Tair #/* Density of dry air, kg / m^3 */
    LHV = 2.501 + (-0.002372727) *Tair #/* This should be MJ kg^-1*/
    LHV = LHV * 1e6 #/* Now it is converted to Joules kg^-1*/
    SlopeFS = TempToSFS(Tair) * 1e-3 #/* kg m^-3 K^-1 */
    SWVP = TempToSWVC(Tair) #/* this is hecto Pascals */
    SWVC = (DdryA * 0.622 * SWVP)/1013.25 #/* This last number is atmospheric pressure in hecto pascals */

    PsycParam =(DdryA * SpecificHeat) / LHV #/* This is in kg m-3 K-1 */

    DeltaPVa = SWVC * (1.0 - RH) #/* kg/m3 */

    ActualVaporPressure = RH * SWVP #/* hecto Pascals */

    totalradiation = Rad * 0.235 #/* This is essentially Watts m^-2 */
    if (totalradiation > 650.0) { 
    #    error stop "total radiation too high"
    }

    Ja = (2.0 * totalradiation * ((1.0 - LeafReflectance - tau) / (1.0 - tau)))

    Ja2 = (2.0 * Iave * 0.235 * ((1.0 - LeafReflectance - tau) / (1.0 - tau)))

    if (WindSpeed < 0.5){ 
        WindSpeed = 0.5
    }

    LayerWindSpeed = WindSpeed
    gvs = stomatacond 
    gvs = gvs * (1.0/41000.0)
    if (gvs <= 0.001){ 
        gvs = 0.001
    }

    Deltat = 0.01
    ChangeInLeafTemp = 10.0

    Counter = 0
    while (ChangeInLeafTemp > 0.5 & Counter <= 10){
    
        OldDeltaT = Deltat

        rlc = 4.0 * StefanBoltzmann * (273.0 + Tair)**3.0 * Deltat  
        ga = leafboundarylayer(LayerWindSpeed, leafw, Airtemperature, Deltat, gvs, ActualVaporPressure)

        PhiN2 = (Ja2 - rlc)  #/* * LeafAreaIndex  */

        TopValue = PhiN2 * (1.0 / ga + 1.0 / gvs) - LHV * DeltaPVa
        BottomValue = LHV * (SlopeFS + PsycParam * (1.0 + ga / gvs))
        Deltat = TopValue / BottomValue #/* This equation is from Thornley and Johnson pg. 418 */
        if (Deltat > 10.0) {
            Deltat = 10.0
        }
        if (Deltat < -10.0){
            Deltat = -10.0
        }
        ChangeInLeafTemp = abs(OldDeltaT - Deltat)

        Counter=Counter+1
    }
     
    PhiN = Ja - rlc

    if (PhiN < 0.0) { 
        PhiN = 0.0
    }

    TransR = (SlopeFS * PhiN + (LHV * PsycParam * ga * DeltaPVa)) / (LHV * (SlopeFS + PsycParam * (1.0 + ga / gvs)))

    EPen = (((SlopeFS * PhiN) + LHV * PsycParam * ga * DeltaPVa)) / (LHV * (SlopeFS + PsycParam))

    EPries = 1.26 * ((SlopeFS * PhiN) / (LHV * (SlopeFS + PsycParam)))

    if (eteq == 1 ) { 
        TransR = EPen
    }
    if (eteq == 2) { 
        TransR = EPries
    }

    ETstr=list(TransR = TransR * 1e6 / 18.0,EPenman = EPen * 1e6 / 18.0,EPriestly = EPries * 1e6 / 18.0,Deltat = Deltat,LayerCond = gvs * 41000.0)
return(ETstr)
}

leafboundarylayer <- function(windspeed,  leafwidth,  AirTemp, deltat,  stomcond,  vappress) {
     Pa = 101325.0,  cf = 1.6361e-3
     leaftemp = AirTemp + deltat
     gsv = stomcond #/* input is in m/s */
     Tak = AirTemp + 273.15 #/* Converts from C to K */
     Tlk = leaftemp + 273.15  #/* Converts from C to K */
     ea = vappress * 1e2 #/* From hPa to Pa */
     ws = windspeed #/* m s^-1 */
     lw = leafwidth #/* meters */

    esTl = TempToSWVC(leaftemp) * 100.0 #/* The function returns hPa, but need Pa */

    gbv_forced = cf *  Tak**0.56 * ((Tak+120.0)*((ws/lw)/Pa))**0.5
    gbv_free = gbv_forced
    eb = (gsv * esTl + gbv_free * ea)/(gsv + gbv_free) #/* Eq 35 */
    Tvdiff = (Tlk / (1.0 - 0.378 * eb/Pa)) - (Tak / (1.0-0.378*ea/Pa)) #/* Eq 34*/

    if (Tvdiff < 0.0){ 
        Tvdiff = -Tvdiff
    }

    gbv_free = cf * Tlk**0.56 * ((Tlk+120)/Pa)**0.5 * (Tvdiff/lw)**0.25

    if (gbv_forced > gbv_free){
        gbv = gbv_forced
    }else{
        gbv = gbv_free
    }
return(gbv)
}

TempToSFS <- function(Temp){
    SlopeFS = 0.338376068 +  0.011435897 * Temp +  0.001111111 * Temp**2
return(SlopeFS)
}

TempToSWVC <- function(Temp){
    a = (18.678 - Temp/234.5) * Temp
    b = 257.14 + Temp
    SWVC = (6.1121 * exp(a/b))
return(SWVC)
}

RHprof <- function(RH, nlayers){
kh = 1.0 - RH
    for (i in 1:nlayers){
        temp_rh = RH * exp(kh * (i/nlayers))
        if (temp_rh > 1){
            temp_rh = 0.99
        }
        relative_humidity_profile[i] = temp_rh
    }
return(relative_humidity_profile)
}

WINDprof <- function(WindSpeed, LAI, nlayers){
 k = 0.7
 LI = LAI / nlayers
    for (i in 1:nlayers){
        CumLAI = LI * i
        wind_speed_profile[i] = WindSpeed * exp(-k * (CumLAI-LI))
    }
return(wind_speed_profile)
}

LNprof <- function(LeafN, LAI, nlayers, kpLN){
    LI = LAI / nlayers
    for (i in 1:nlayers){
        CumLAI = LI * i
        leafN_profile[i]=LeafN * exp(-kpLN * (CumLAI-LI))
    }
return(leafN_profile)
}

lightME <- function(lat,doy,cza){
atmo_trans=0.85
atmo_press=1e5
local_atmo=1e5
press_r=local_atmo/atmo_press
prop_scat=0.3
dir_frac=atmo_trans**(press_r/cza)
diff_frac=prop_scat*(1.0-dir_frac)*cza
if(cza<=0){
dir_frac=0.0
diff_frac=1.0
}
lightMEstr=list(dir_frac=dir_frac/(dir_frac+diff_frac),diff_frac=diff_frac/(dir_frac+diff_frac))
return(lightMEstr)
}

sunML <- function(soldir,soldiff,lai,nlayers,cosTheta,kd,chil,heightf){ 
alphascatter = 0.8
theta=acos(min(max(cosTheta,-1.0),1.0))#Yufeng: make sure acos takes a value between -1 and 1
k0=sqrt(chil**2.0+(tan(theta))**2.0)
k1=chil+1.744*(chil+1.183)**(-0.733)
k=k0/k1#k Foliar absorption coefficent
if (k < 0.0){
  k = -k
}
laii=lai/nlayers
Ibeam=soldir*cosTheta
Isolar=Ibeam*k
for (i in 1:nlayers){
    if (cosTheta<=1e-10){
        Isolar=soldir/k1
        Idiff=soldiff*exp(-kd*cumlai)
        Fsun=0.0
        Fshade=0.0
        Iave=0.0
    }else{
        cumlai=laii*(i-1.0+0.5)
        Iscat=Ibeam*(exp(-k*sqrt(alphascatter)*cumlai)-exp(-k*cumlai))
        Idiff=soldiff*exp(-kd*cumlai)+Iscat #kd: Light extinction coefficient for diffuse light
        Ls=(1.0-exp(-k*laii))*exp(-k*cumlai)/k
        Ld=laii-Ls
        Fsun=Ls/(Ls+Ld)
        Fshade=Ld/(Ls+Ld)
        Iave=(Fsun*(Isolar+Idiff)+Fshade*Idiff)*(1.0-exp(-k*laii))/k    
    }
    dir_irr[i]=Isolar+Idiff
    diff_irr[i]=Idiff
    tot_irr[i]=Iave
    sunlit_frac[i]=Fsun
    shade_frac[i]=Fshade
    height[i]=(lai-cumlai)/heightf
}
sunMLstr=list(dir_irr=dir_irr,diff_irr=diff_irr,tot_irr=tot_irr,sunlit_frac=sunlit_frac,shade_frac=shade_frac,height=height)
return(sunMLstr)
}

 soilML_yh <- function(cws, transp, soildepth, depths,fieldc,  wiltp,  phi1,  phi2,  wsFun,
         layers,  rootDB,  LAI,  k,  AirTemp,  IRad,  winds,
         RelH,  rfl,  rsec,  rsdf){ 
wsPhotoCol=0.0
wsSpleafCol=0.0
root_distribution = 0.0
rootDepth = rootDB * rsdf
if (rootDepth > soildepth){
   rootDepth = soildepth
}
nlayergtbedrock=0
for (i in 1:layers){
    if (depths[i+1]<=soilDepth) {
	nlayergtbedrock=nlayergtbedrock+1
     }
}
    root_distribution[1:nlayergtbedrock] = rootDist(nlayergtbedrock,rootDepth,depths[1:nlayergtbedrock+1],rfl)
    for (i in 1:nlayergtbedrock){
	layerDepth=depths[i+1]-depths[i]
	aw = cws[i] * layerDepth
	paw = aw - wiltp[i] * layerDepth
      if (paw < 0){
            paw = 0.0 
      }
	if (i == 1){
          Sevap = SoilEvapo(LAI,k,AirTemp,IRad,cws[i],fieldc[i],wiltp[i],winds,RelH,rsec)
          Ctransp = transp*root_distribution[i]
          EvapoTra = Ctransp + Sevap
          Newpawha = (paw * 1.e4) - EvapoTra / 0.9982 #/* See the watstr function for this last number 0.9882 */
        }else{
            Ctransp = transp*root_distribution[i]
            EvapoTra = Ctransp
            Newpawha = (paw * 1.e4) - (EvapoTra / 0.9982 + oldEvapoTra)
        }
        if (Newpawha < 0) {
        oldEvapoTra = -Newpawha
        aw  = wiltp[i] * layerDepth 
        awc = wiltp[i]
        }else{
        oldEvapoTra = 0.0
        paw = Newpawha / 1e4
        awc = paw / layerDepth + wiltp[i]  
        }
	soilMLstr$cws[i] = awc
	if (wsFun == 0){
            slp = 1.0/(fieldc[i] - wiltp[i])
            intcpt = 1.0 - fieldc[i] * slp
            wsPhoto = slp * awc + intcpt
        }else if (wsFun == 1){
            phi10 = (fieldc[i] + wiltp[i])/2.0
            wsPhoto = 1.0/(1.0 + exp((phi10 - awc)/ phi1))
        }else if (wsFun == 2){
            slp = (1.0 - wiltp[i])/(fieldc[i] - wiltp[i])
            intcpt = 1.0 - fieldc[i] * slp
            theta = slp * awc + intcpt
            wsPhoto = (1.0 - exp(-2.5 * (theta - wiltp[i])/(1.0 - wiltp[i]))) / (1.0 - exp(-2.5))
        }else if (wsFun == 3){
            wsPhoto = 1.0
        }

        if (wsPhoto <= 0.0){ 
            wsPhoto = 1e-2 #/* This can be mathematically lower than zero in some cases but I should prevent that. */
        }
       if (wsPhoto > 1.0) {  #Yufeng
           wsPhoto = 1.0
       }
           wsPhotoCol = wsPhotoCol +  wsPhoto
       
       wsSpleaf = awc**phi2 * 1.0/(fieldc[i]**phi2) 
      if (wsFun == 3){
            wsSpleaf = 1.0
       }
       if (wsSpleaf > 1.0) {  
          wsSpleaf = 1.0
       }
          wsSpleafCol = wsSpleafCol + wsSpleaf#*root_distribution(i)/sum(root_distribution)
  }
    soilMLstr$rcoefPhoto = (wsPhotoCol)/nlayergtbedrock
    soilMLstr$rcoefSpleaf = (wsSpleafCol)/nlayergtbedrock
    soilMLstr$SoilEvapo = Sevap
}

SoilEvapo <- function(LAI,  k,  AirTemp,  IRad,
         awc,  fieldc,  wiltp,  winds,  RelH,  rsec ) {
method = 1
cf2 = 3600 * 1e-3 * 18 * 1e-6 * 10000
SoilClodSize = 0.04
SoilReflectance = 0.2
SoilTransmission = 0.01
SpecificHeat = 1010.0
StefanBoltzman = 5.67e-8 
    SoilArea = exp(-k * LAI)
    SoilTemp = AirTemp
    rawc = (awc - wiltp)/(fieldc - wiltp)
        if(rawc < 1e-7) {rawc = 1e-7}
        if(rawc > 1.0) {rawc = 1.0}

     if(awc > fieldc) {
                Up =  1.0
     }else{
                Up = 1.0 - (1.0 + 1.3 * rawc)**(-5.0)
     }
       
    IRad = IRad * rsec # Radiation soil evaporation coefficient  

    TotalRadiation = IRad * 0.235

    DdryA = 1.295163636 + (-0.004258182) *(AirTemp)
    LHV = (2.501 + (-0.002372727) * AirTemp )* 1e6 
    SlopeFS = TempToSFS(AirTemp) * 1e-3
    SWVC = TempToSWVC(AirTemp) * 1e-3
    
    PsycParam = (DdryA * SpecificHeat) / LHV
    DeltaPVa = SWVC * (1.0 - RelH) 

    BoundaryLayerThickness = 4e-3 * sqrt(SoilClodSize / winds) 
    DiffCoef = 2.126e-5 * 1.48e-7 * SoilTemp
    SoilBoundaryLayer = DiffCoef / BoundaryLayerThickness

    Ja = 2.0 * TotalRadiation * ((1.0 - SoilReflectance - SoilTransmission) / (1.0 - SoilTransmission))

    rlc = 4.0 * StefanBoltzman * ((273.0 + SoilTemp))**3.0 * 0.005

    PhiN = Ja - rlc # Calculate the net radiation balance
    if (PhiN < 0){ 
        PhiN = 1e-7
    }

    if (method == 0){ 
        Evaporation = 1.26 * (SlopeFS * PhiN) / (LHV * (SlopeFS + PsycParam))
     }else{
        # Penman-Monteith 
        Evaporation = (SlopeFS * PhiN + LHV * PsycParam * SoilBoundaryLayer * DeltaPVa) / (LHV * (SlopeFS + PsycParam))
    }
    Evaporation = Evaporation * 1.e6/18.
    Evaporation = Evaporation * SoilArea * Up * cf2 
    if (Evaporation < 0 ){
        Evaporation = 1e-6
    }
}

rootDist <- function(layer, rootDepth, depthsp, rfl){
CumRootDist=1
ca=0.0
CumLayerDepth=0.0
   for(i in 1:layer){
            layerDepth = depthsp[i+1] - depthsp[i]

        CumLayerDepth = CumLayerDepth + layerDepth

        if (rootDepth > CumLayerDepth){
            CumRootDist=CumRootDist+1
        }
    }

    for (j in 1:layer){
        if (j <= CumRootDist) {
            a = poisson_density(j,CumRootDist*rfl)
            rootDistri[j] = a
            ca = ca + a
        }else{
            rootDistri[j] = 0.0
        }
    }

    for (k in 1:layer){
        rootDistri[k] = rootDistri[k]/ ca
    }
}

poisson_density <- function(x, lambda){
e=2.71828
factorial_x = sqrt(2.0 * pi * x) * (x / e)**x
log_result = -lambda + x * log(lambda) - log(factorial_x)
}

sel_dbp_coef <- function(coefs,  TherPrds,  TherTime){
#coefs: length 25
#TherPrds: length 6
    biomassPartistr$kLeaf = 0.0
    biomassPartistr$kStem = 0.0
    biomassPartistr$kRoot = 0.0
    biomassPartistr$kRhiz = 0.0
    biomassPartistr$kGrain = 0.0 

    if (TherTime <= TherPrds[1]){
    
        biomassPartistr$kStem = coefs[1]
        biomassPartistr$kLeaf = coefs[2]
        biomassPartistr$kRoot = coefs[3]
        biomassPartistr$kRhiz = coefs[4]

    }else if (TherTime <= TherPrds[2]){
    
        biomassPartistr$kStem = coefs[5]
        biomassPartistr$kLeaf = coefs[6]
        biomassPartistr$kRoot = coefs[7]
        biomassPartistr$kRhiz = coefs[8]

    }else if (TherTime <= TherPrds[3]) {
    
        biomassPartistr$kStem = coefs[9]
        biomassPartistr$kLeaf = coefs[10]
        biomassPartistr$kRoot = coefs[11]
        biomassPartistr$kRhiz = coefs[12]

    }else if (TherTime <= TherPrds[4]) {
    
        biomassPartistr$kStem = coefs[13]
        biomassPartistr$kLeaf = coefs[14]
        biomassPartistr$kRoot = coefs[15]
        biomassPartistr$kRhiz = coefs[16]

    }else if (TherTime <= TherPrds[5]) {
        
        biomassPartistr$kStem = coefs[17]
        biomassPartistr$kLeaf = coefs[18]
        biomassPartistr$kRoot = coefs[19]
        biomassPartistr$kRhiz = coefs[20]

    }else{ 

        biomassPartistr$kStem = coefs[21]
        biomassPartistr$kLeaf = coefs[22]
        biomassPartistr$kRoot = coefs[23]
        biomassPartistr$kRhiz = coefs[24]
        biomassPartistr$kGrain = coefs[25]
    }
}

leaf_n_limitation <- function(kLn, leaf_n_0, current_leaf,current_stem){ 
    leaf_n = leaf_n_0 * (current_leaf + current_stem)**(-kLn)
    if (leaf_n>leaf_n_0) {
        leafN=leaf_n_0
    }else{
        leafN=leaf_n
    }
}

resp <- function( comp,  mrc,  temp){
    resp = comp *  (1.0 - (mrc * 2.0**(temp/10.0)))
    if (resp<0){ 
        resp = 0.
    }
}

GroX <- function(lat,doy,cosz,solar_dir,solar_diff,temp,rh,windspeed,precip,kd,chil,leafwidth,et_equation,heightf,nlayers,initial_biomass,sencoefs,timestep,iSp,SpD,dbpcoefs,thermalp,tbase,vmax1,alpha1,kparm,theta,beta,Rd,Catm,b0,b1,soilcoefs,ileafn,kLN,vmaxb1,alphab1,mresp,wsFun,ws,soilLayers,soilDepths,cws,secs,kpLN,lnb0,lnb1,lnfun,upperT,lowerT,nitroP,StomataWS,TTc,Grain,Leaf,Stem,Root,Rhizome,LeafN,Sp,vmax,alpha,doy0,doy1,iWatCont,hour,fieldc,wiltp,bedrock){
    LeafN_0 = ileafn
    seneRate = 0.0005 
    mrc1 = mrespr[1]
    mrc2 = mresp[2]

    phi1 = soilcoefs[3]
    phi2 = soilcoefs[4]
    soilDepth = soilcoefs[5]
    soilDepth = min(soilDepth,bedrock)
    
    seneLeaf = sencoefs[1]
    seneStem = sencoefs[2]
    seneRoot = sencoefs[3]
    seneRhizome = sencoefs[4]
    
    
    rfl = secs[1]  # root factor lambda 
    rsec = secs[2] # radiation soil evaporation coefficient 
    rsdf = secs[3] # root soil depth factor 
    scsf = soilcoefs[7] # stomatal conductance sensitivity factor */ /* Rprintf("scsf !.2f",scsf) */
    transpRes = soilcoefs[8] # Resistance to transpiration from soil to leaf 
    leafPotTh = soilcoefs[9] # Leaf water potential threshold   
    
    if(doy==doy0 .and. (abs(hour-NINT(hour)) < 0.05)) {
    TTc=0.0
    Grain=0.0
    Rhizome = initial_biomass[1]
    Stem = initial_biomass[2]
    Leaf = initial_biomass[3]
    Root = initial_biomass[4]
    LeafN = ileafn # Need to set it because it is used by CanA before it is computed
    Sp=iSp
    vmax = vmax1
    alpha = alpha1
    StomataWS=1.0
    }     
   LAI = min(10.0,Leaf * Sp)    

        if (temp > tbase){ 
            TTc = TTc + (temp-tbase) / (24./timestep) 
        }

        Canopy = CanAC(LAI, doy, cosz,
                solar_dir,solar_diff, temp, rh, windspeed,
                lat, nlayers, vmax, alpha, kparm, beta,
                Rd, Catm, b0, b1, theta, kd, chil,
                heightf, LeafN, kpLN, lnb0, lnb1, lnfun, upperT, lowerT,
                nitroP, leafwidth, et_equation, StomataWS, ws)

        CanopyA = Canopy$Assim * timestep
        CanopyT = Canopy$Trans * timestep
        soilMLS = soilML_yh(cws, CanopyT,soilDepth,soilDepths, fieldc, wiltp, phi1, phi2, wsFun,
                              soilLayers, Root, LAI, 0.68, temp, solar_dir+solar_diff, windspeed,
                              rh, rfl, rsec, rsdf)
             cws = soilMLS$cws
             StomataWS = soilMLS$rcoefPhoto
             LeafWS = soilMLS$rcoefSpleaf
             soilEvap =soilMLS$SoilEvapo
        
        dbpS = sel_dbp_coef(dbpcoefs, thermalp, TTc)

        kLeaf = dbpS$kLeaf
        kStem = dbpS$kStem
        kRoot = dbpS$kRoot
        kGrain = dbpS$kGrain
        kRhizome = dbpS$kRhiz


        current_leaf=Leaf
        current_stem=Stem
        LeafN = leaf_n_limitation(kLN, LeafN_0, current_leaf,current_stem)

        vmax = (LeafN_0 - LeafN) * vmaxb1 + vmax1
        alpha = (LeafN_0 - LeafN) * alphab1 + alpha1


        MinNitro = MinNitro - LeafN * (Stem + Leaf) * 1e-1
        if (MinNitro < 0){ 
            MinNitro = 1e-3
         }

        if (kLeaf > 0) {
            newLeaf = CanopyA * kLeaf * LeafWS
            newLeaf = resp(newLeaf, mrc1, temp)
            newLeafcol = newLeaf
        }else{
            newLeaf = Leaf * kLeaf
            Rhizome = Rhizome + kRhizome * (-newLeaf) * 0.9 #0.9 is the efficiency of retranslocation
            Stem = Stem + kStem * (-newLeaf) * 0.9
            Root = Root + kRoot * (-newLeaf) * 0.9
            Grain = Grain + kGrain * (-newLeaf) * 0.9
        }

        if (TTc < seneLeaf) { 
            Leaf = Leaf + newLeaf
        }else{ 
		    Leaf = (Leaf + newLeaf)- seneRate * Leaf
                    Remob = Leaf * seneRate*0.6
                    Rhizome = Rhizome + kRhizome * Remob
                    Stem = Stem + kStem * Remob
                    Root = Root + kRoot * Remob
                    Grain = Grain + kGrain * Remob
         }

        if (mod(NINT(hour),24) == 0) { 
            Sp = iSp - (doy - doy0) * SpD
        }

        LAI = Leaf * Sp

        if (kStem >= 0){ 
            newStem = CanopyA * kStem
            newStem = resp(newStem, mrc1, temp)
            newStemcol = newStem
        }else{
          break("kStem should be positive")
        }
        if (TTc < seneStem){ 
            Stem = Stem + newStem
	}else{
		    Stem = (Stem + newStem) - seneRate * Stem
        } 
        
        if (kRoot > 0){
            newRoot = CanopyA * kRoot
            newRoot = resp(newRoot, mrc2, temp)
            newRootcol = newRoot
        }else{
            newRoot = Root * kRoot
            Rhizome = Rhizome + kRhizome * (-newRoot) * 0.9
            Stem = Stem + kStem * (-newRoot) * 0.9
            Leaf = Leaf + kLeaf * (-newRoot) * 0.9
            Grain = Grain + kGrain * (-newRoot) * 0.9
        }

        if  (TTc < seneRoot){ 
            Root = Root + newRoot
        }else{
               Root=(Root+newRoot)- seneRate * Root
        }

        if (kRhizome > 0) {
            newRhizome = CanopyA * kRhizome
            newRhizome = resp(newRhizome, mrc2, temp)
            newRhizomecol = newRhizome
        }else{
            if (Rhizome < 0){
                Rhizome = 1e-4
                print("warning:Rhizome became negative")
            }

            newRhizome = Rhizome * kRhizome
            Root = Root + kRoot * (-newRhizome)
            Stem = Stem + kStem * (-newRhizome)
            Leaf = Leaf + kLeaf * (-newRhizome)
            Grain = Grain + kGrain * (-newRhizome)
        }

        if  (TTc < seneRhizome){ 
            Rhizome = Rhizome + newRhizome
        }else{
               Rhizome=(Rhizome+newRhizome) - seneRate * Rhizome
        }       
 
        if (kGrain < 1e-10 .or. TTc < thermalp[5]){
            newGrain = 0.0
            Grain = Grain + newGrain
        }else{
            newGrain = CanopyA * kGrain
            Grain = Grain + newGrain
        }
    
    bio_out$TTc=TTc
    bio_out$Grain=Grain
    bio_out$Leaf=Leaf
    bio_out$Stem=Stem
    bio_out$Root=Root
    bio_out$Rhizome=Rhizome
    bio_out$LeafN=LeafN
    bio_out$Sp=Sp
    bio_out$vmax=vmax
    bio_out$alpha=alpha
    bio_out$stomataWS=StomataWS
    bio_out$cws=cws         #!m3/m3
    bio_out$biolai=LAI
    bio_out$bioev=soilEvap  #Mg ha-1 hr-1
    bio_out$biotr=CanopyT   #Mg ha-1 hr-1
    bio_out$CA = CanopyA

} 
