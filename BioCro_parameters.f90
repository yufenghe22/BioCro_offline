module biocro_parameters

integer, parameter :: et_equation = 0, & !0='Penman-Monteith' 1='Penman' 2='Priestly'
                      nlayers = 10, &
                      timestep = 1, &
                      soilLayers=11, &
                      wsFun=0, &
                      ws=1, &
                      lnfun=0
real, parameter :: kd=0.1, &
                   chil=1.0, & 
                   leafwidth=0.04, &
                   iSp=1.7, &
                   SpD=0.0, &
                   heightf=3.0, &
                   tbase=0.0, &
                   vmax1=39.0, &
                   alpha1=0.04, &
                   kparm=0.7, &
                   theta=0.83, &
                   beta=0.93, &
                   Rd=0.8, &
                   Catm=380.0, &
                   b0=0.08, &
                   b1=3.0, &
                   ileafn=2.0, &
                   kLN=0.5, &
                   vmaxb1=0.0, &
                   alphab1=0.0, &
                   kpLN=0.2, &
                   lnb0=-5.0, &
                   lnb1=18.0, &
                   upperT=37.5, &
                   lowerT=3.0                      
real, dimension(4),parameter  :: initial_biomass = (/7.0,7.0*1e-3,7.0*1e-4,7.0*1e-3/), &      
                                sencoefs = (/3000.0,3500.0,4000.0,4000.0/)                   
real, dimension(25),parameter ::  dbpcoefs=(/0.37,0.33,0.3,-8e-4,0.85,0.14,0.01,-5e-4,0.63,0.01,0.01,0.35,0.63,&
 0.01,0.01,0.35,0.63,0.01,0.01,0.35,0.63,0.01,0.01,0.35,0.0/)         
real, dimension(5),parameter  :: thermalp=(/562.0,1312.0,2063.0,3211.0,7000.0/)

real, dimension(soilLayers+1),parameter  :: soilDepths=(/0.0,1.75e-2,4.51e-2,9.05e-2,0.1655,0.2891,0.4929,0.8288,& 
 1.3828,2.296,3.8018,5.676/) ! used for multiple-layer soil ML        
real, dimension(9),parameter  :: soilcoefs=(/9999.0,9999.0,0.01,10.0,soilDepths(soilLayers+1),9999.0,1.0,5.e6,-800.0/)!filedC,WiltP                                          

real, dimension(2),parameter  :: mresp=(/0.02,0.03/)

real, dimension(3),parameter  :: secs=(/0.2,0.2,0.44/)

real, dimension(5),parameter  :: doy0_coefs=(/-0.8567264318,-1.79315234,0.0490230963, -0.006801856,0.012252882/),&
                                 doy1_coefs=(/12.874761643, -2.734511704, -0.138748429, -0.005066692,0.052052422/)
end module biocro_parameters
