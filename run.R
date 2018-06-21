library(ncdf4)

et_equation = 0
nlayers = 10
timestep = 1
soilLayers=11
wsFun=0
ws=1
lnfun=0
kd=0.1
chil=1.0
leafwidth=0.04
iSp=1.7
SpD=0.0
heightf=3.0
tbase=0.0
vmax1=39.0
alpha1=0.04
kparm=0.7
theta=0.83
beta=0.93
Rd=0.8
Catm=380.0
b0=0.08
b1=3.0
ileafn=2.0
kLN=0.5
vmaxb1=0.0
alphab1=0.0
kpLN=0.2
lnb0=-5.0
lnb1=18.0
upperT=37.5
lowerT=3.0                      
initial_biomass = c(7.0,7.0*1e-3,7.0*1e-4,7.0*1e-3)      
sencoefs = c(3000.0,3500.0,4000.0,4000.0) 
dbpcoefs=c(0.37,0.33,0.3,-8e-4,0.85,0.14,0.01,-5e-4,0.63,0.01,0.01,0.35,0.63,0.01,0.01,0.35,0.63,0.01,0.01,0.35,0.63,0.01,0.01,0.35,0.0)
thermalp=c(562.0,1312.0,2063.0,3211.0,7000.0)
soilDepths=c(0.0,1.75e-2,4.51e-2,9.05e-2,0.1655,0.2891,0.4929,0.8288,1.3828,2.296,3.8018,5.676)        
soilcoefs=c(9999.0,9999.0,0.01,10.0,soilDepths[soilLayers+1],9999.0,1.0,5.e6,-800.0)
mresp=c(0.02,0.03)
secs=c(0.2,0.2,0.44)
doy0_coefs=c(-0.8567264318,-1.79315234,0.0490230963, -0.006801856,0.012252882)
doy1_coefs=c(12.874761643, -2.734511704, -0.138748429, -0.005066692,0.052052422)
source("BioCro.R")

years=2000:2007
for (i in 1:length(years)){
year=years[i]
sday=paste(year,"/1/1",sep="")
eday=paste(year,"/12/31",sep="")
dates=seq(as.Date(sday),as.Date(eday),by="days")
t2m=array(0,dim=c(195,138,length(dates)*8))
pr =t2m

for (j in 1:length(dates)){
date=dates[j]
fname=paste("wrfout_d01_",date,"_00:00:00",sep="")
path=paste("~/scratch/cwrfgo_ctl/ctl/",fname,sep="")
f=nc_open(path)
nc_close(f)

if (doy >= doy0 & doy <= doy1){ 
bio_out=GroX(XLAT[i,j],doy,xcoszn[i,j],solar_dir,solar_diff,temp,rh,windspeed,precip,
    kd,chil,leafwidth,et_equation,heightf,nlayers,initial_biomass,sencoefs,
    timestep,iSp,SpD,dbpcoefs,thermalp,tbase,vmax1,alpha1,kparm,theta,
    beta,Rd,Catm,b0,b1,soilcoefs,ileafn,kLN,vmaxb1,alphab1,mresp,wsFun,
    ws,soilLayers,soilDepths,cws[i,,j],
    secs,kpLN,lnb0,lnb1,lnfun,upperT,lowerT,nitroP,stomataWS[i,j],
    TTc[i,j],Grain[i,j],Leaf[i,j],Stem[i,j],Root[i,j],Rhizome[i,j],
    LeafN[i,j],Sp[i,j],vmax[i,j],alpha[i,j],
    doy0,doy1,iWatCont,hour,fieldc,wiltp,BEDROCK[i,j])

    }
}
}
