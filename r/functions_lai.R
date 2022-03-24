lai.drop.opt.func <- function (VPD,#kPa
                               E,#starting soil availibale water; mm
                               PAR, #Annual PAR; mj par m-2 ground d-1
                               TMAX, #mean Tmax of each month; Celsius
                               lai.initial, #starting lai; m2 m-2
                               g1 = 4, #slope of A/E
                               r.turnover = 0.1,# leaf turnover rate; fraction
                               beta.lai = 1,#power of beta function
                               beta.g1 = 2,#power of beta function
                               beta.jmax = 2,#power of beta function
                               sm.cap = 250, #soil water avalibility at capacity; mm
                               # sm.wilt =C 50, #wilting pointmm
                               theta = 0.7, #medlyn 2007
                               budyco = F,#when true, use budyco curve and E is a fraction of MAP; Otherwise E=MAP
                               gammaS25 = 42.75, #ppm; 45 as Mcmurtire or 42.75 as Bernacchi 2001
                               k = 0.5, #light extinction coeffcient; sands 1995,zhang 2014
                               Ca = 375, #atmospheric co2 concentration;ppm
                               # Jmax25  = 1.8 * 70.0, # mumol s-1 based on leaf N content = 1.8 g N m-2 leaf
                               Jmax25 = 160, #synthetic data set
                               m = 0, #leaf transmittion coefficient
                               P = 100, #atmospheric pressure to convert VPD unit, Kpa
                               Rlight = 1.0, #mumol m-2 s-1
                               Rdark = 2.0, #mumol m-2 s-1
                               SLA = 52.6,#specific leaf are; cm2 leaf g-1 dry mass
                               R.test=F,#when true Rdark is a function of MAP otherwise uses ginven value
                               lma.test=F,#when true LMA is a function of PAR otherwise based on SLA given
                               r.min = 0.17 #minimum leaf respiration at abient temp; mumol C m-2 s-1; Atkin 2015
){
  mole.weight.c <- 12.0
  g2ummole.c <- 1.0/12.0*10^6
  mole.weight.ch <- 1.6
  mumole2mole <- 10^-6
  # days <- 365.25
  h  <- 12.0 * 60.0 * 60.0 # second of day length in a year
  mole2mm.water <- 1.8 * 10^-2
  construction.r.constant <- 1.3
  convertor <- 2.0 * 10^6 # 2 mumol PAR MJ-1 (Sands,1996);
  MOL2MJ <- 4.0 #converts mol of PAR to MJ; 4 MJ mol-1
  Q <- 2.0 * PAR  #Daily irradiance; Mj total irradiance m-2 ground d-1
  
  # fraction based on Budyco curve (Zhang et al, 2001)
  # fraction of MAP available to plant
  ifelse(budyco == TRUE,
         fraction <- (E + 2.0*1410.0) / (1.0 + 2.0*1410.0 / E + E/1410.0) / E,
         fraction <- 1)
  
  # rdard as measured as a function of log(map) also from globresp dataset
  ifelse(R.test==TRUE,Rdark <- -0.84 * log(E) + 7.25,Rdark <- Rdark) #rd based on wright 2006
  # set minimum r dark
  rmin <- r.min #minimum leaf respiration at abient temp; mol C m-2 s-1; Atkin 2015
  
  if(is.na(Rdark) == FALSE){
    if(Rdark < rmin) Rdark <- rmin
  }
  
  # gammaS value as func of TMAX (Bernacchi,2001)
  gammaS <- exp(19.02-37.83/(8.314*10^-3*(TMAX + 273.15)))
  
  # functions to get LUE ####
  CiFunc <- function (g1) {
    Ca * g1 / (g1 + VPD^0.5)
  }
  
  
  #photo rate at saturate light on the top of canopy
  AxFunc <- function(g1){
    Jmax / 4.0 * (CiFunc(g1) - gammaS) / (CiFunc(g1) + 2.0 * gammaS)
  }
  
  # quntum effect (Ci function cite Belinda's 2000 CJFR paper;medlyn 2007)
  alphaFunc <- function (x){
    0.26/4.0 * (Ca - gammaS) / (Ca + 2.0 * gammaS)
  }
  
  #q=(π∙k∙α∙Q*gamma)/(2h(1-m)A_x ) Sands 1995 Eq 13
  qFunc <- function(x){
    (pi * k * alphaFunc(x) * Q * convertor) / (2.0 * h * (1.0-m) * AxFunc(x))
  }
  
  # mast mode lue
  Epsilon <- function(g1){
    qq <- qFunc(g1)
    
    sin1 <- sin(pi / 24)
    sin2 <- sin(pi * 3 / 24)
    sin3 <- sin(pi * 5 / 24)
    sin4 <- sin(pi * 7 / 24)
    sin5 <- sin(pi * 9 / 24)
    sin6 <- sin(pi * 11 / 24)
    g1 <- sin1 / (1 + qq * sin1 + sqrt((1 + qq * sin1) ^ 2 - 4 * theta * qq * sin1))
    g2 <- sin2 / (1 + qq * sin2 + sqrt((1 + qq * sin2) ^ 2 - 4 * theta * qq * sin2))
    g3 <- sin3 / (1 + qq * sin3 + sqrt((1 + qq * sin3) ^ 2 - 4 * theta * qq * sin3))
    g4 <- sin4 / (1 + qq * sin4 + sqrt((1 + qq * sin4) ^ 2 - 4 * theta * qq * sin4))
    g5 <- sin5 / (1 + qq * sin5 + sqrt((1 + qq * sin5) ^ 2 - 4 * theta * qq * sin5))
    g6 <- sin6 / (1 + qq * sin6 + sqrt((1 + qq * sin6) ^ 2 - 4 * theta * qq * sin6))
    
    #Trapezoidal rule - seems more accurate
    gg <- 1.0/6.0 * (g1 + g2 + g3 + g4 + g5 + g6)
    eps <- alphaFunc(g1) * gg * pi
    return(eps)
  }
  
  LUE <- function(g1) {
    Epsilon(g1) * MOL2MJ * mole.weight.c
  }
  
  # instrinsic tranpiration effciency as in Medlyn 2011
  ITE <- function(g1){
    1 / mole2mm.water * mole.weight.c * mumole2mole / mole.weight.ch * Ca *P / (VPD + g1 * VPD^0.5)
  }
  
  # lai.vec <- c()
  apar.func <- function(PAR,LAI,k=0.5){
    PAR * (1 - exp(-k*LAI))
  }
  # water use
  tran.func <- function(LAI,g1,k=0.5){
    apar <-  apar.func(PAR,LAI)
    return(LUE(g1) * apar  / ITE(g1))
  }
  

  # basiclly mol of water times ite
  GPPfunc <- function(g1,LAI){
    LUE(g1) * apar.func(PAR,LAI)
  }
  # npp
  netfunc <- function(g1,LAI){
    GPPfunc(g1,LAI) - mole.weight.c * mumole2mole  * 2 * h * LAI* Rdark
  }
  # 
  i = 2
  lai.tmp <- lai.initial
  swc.norm <- E  / sm.cap 
  water.remain = E
  npp.vec <- c()
  g1.tmp <- g1*swc.norm^beta.g1
  Jmax <- Jmax25 * swc.norm^ beta.jmax
  lai.tmp <- lai.initial
  npp.vec[1] <- netfunc(LAI = lai.initial,g1 = g1*swc.norm^beta.g1)
  # tracking 
  gpp.vec <- c()
  gpp.vec[1] <- GPPfunc(g1 = g1,LAI = lai.initial)
  g1.vec <- c()
  g1.vec[1] <- g1
  lai.vec <- c()
  lai.vec[1] <- lai.initial
  soil.water <- c()
  soil.water[1] <- E
  # take a g1 threshold of 88%
  while(swc.norm^beta.g1>0.12&lai.tmp>0.1){
    
    tran.tmp <- tran.func(LAI = lai.tmp,g1 = g1.tmp)
    # normlised swc
    water.remain <- water.remain - tran.tmp
    swc.norm <- water.remain  / sm.cap 
    # 
    g1.tmp <- g1*swc.norm^beta.g1
    Jmax <- Jmax25 * swc.norm^ beta.jmax
    lai.tmp <- lai.tmp - lai.tmp*r.turnover*(1-swc.norm)^beta.lai
    
    npp.vec[i] <- netfunc(g1 = g1.tmp,LAI = lai.tmp)
    gpp.vec[i] <- GPPfunc(g1 = g1.tmp,LAI = lai.tmp)
    g1.vec[i] <- g1.tmp
    lai.vec[i] <- lai.tmp
    soil.water[i] <- water.remain
    # day move on
    i=i+1
  }

 return(data.frame(npp = npp.vec,
                   gpp = gpp.vec,
                   g1 = g1.vec,
                   lai = lai.vec,
                   soil.water=soil.water
                   ))
  
  
  
  # # get LAI based on E*ITE and LUE####
  # LAI <- function(g1){
  #   -1.0/k * log(1.0 - fraction * E * ITE(g1) / LUE(g1) / PAR)
  # }
  # 
  # # basiclly mol of water times ite
  # GPPfunc <- function(g1){
  #   LUE(g1) * PAR * (1-exp(-k*LAI(g1)))
  # }
  # 
  # # gs  mol H2O m-2 leaf s-1
  # gs <- function(g1){
  #   mole.weight.ch * (1 + g1/VPD^0.5) * GPPfunc(g1) * g2ummole.c  / Ca / LAI(g1) / h 
  # }
  # 
  # # construction based on SLA following Manzonni et al. 2015;
  # # SLA based on Duursma et al. 2015:
  # LMA <- 1.0/SLA*10^4 #g Dry mass m-2 leaf
  # # LMA LL based on Wright 2005
  # # RAD is W m-2 (daily average) and par is MJ m-2 yr-1 (annual total)
  # # to convert PAR to RAD (1w = 1*10^-6 MJ/s):
  # mj2w <- 10^6
  # RAD <- PAR / days / h * mj2w
  # ifelse(lma.test==TRUE, LMA <- 10^(0.0039 * RAD + 1.43), LMA <- LMA) #lma based on wright 2005
  # 
  # # leaf lifespan in year
  # LL <- 10^(1.14*log10(LMA)-2.45) #yr
  # LCA <- 0.5*LMA  #g C m-2 leaf
  # #optimasation target
  # netfunc <- function(g1){
  #   GPPfunc(g1) -  #construction.r.constant * LCA * LAI(g1) / LL -
  #     mole.weight.c * mumole2mole  * 2 * h * LAI(g1)* Rdark
  # }
  # 
  # # define outputs
  # if (is.na(E) == TRUE | is.na(VPD) == TRUE | is.na(PAR) == TRUE | is.na(TMAX) == TRUE)
  #   optimal.g1 <- NA
  # else
  #   optimal.g1 <- optimise(netfunc,interval = c(0.5,25),maximum=TRUE)$maximum
  # 
  # if (is.na(optimal.g1)== TRUE){
  #   optimal.LAI <- NA
  #   optimal.ITE <- NA
  #   optimal.gs <- NA
  #   optimal.NCG <- NA
  #   optimal.LUE <- NA
  #   optimal.GPP <- NA
  #   optimal.WUE <- NA
  #   
  # }else{
  #   optimal.LAI <- LAI(optimal.g1)
  #   optimal.ITE <- ITE(optimal.g1)
  #   optimal.gs <- gs(optimal.g1)
  #   optimal.NCG <- netfunc(optimal.g1)
  #   optimal.LUE <- LUE(optimal.g1)
  #   optimal.GPP <- GPPfunc(optimal.g1)
  #   optimal.WUE <- ITE(optimal.g1)
  # }
  # 
  # # output vector
  # v <- c(g1=optimal.g1,
  #        LAI=optimal.LAI,
  #        gs=optimal.gs,
  #        NCE=optimal.NCG,
  #        LUE = optimal.LUE,
  #        GPP = optimal.GPP,
  #        WUE = optimal.WUE)
  # return(v)
}


