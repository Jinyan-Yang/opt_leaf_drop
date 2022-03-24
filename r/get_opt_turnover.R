tmp.df <- lai.drop.opt.func(VPD = 2,E = 200,PAR = 8,TMAX = 30,
                            lai.initial = 1.5,g1 = 4,r.turnover = 0.03,beta.lai = 0.2 
                  ,beta.g1 = 0.38,beta.jmax = 0.38)
plot(tmp.df$lai)
plot(tmp.df$gpp)
library(DEoptim)
model.de.func <- function(pars,
                          VPD = 1.5,E = 200,PAR = 8,TMAX = 30,
                          lai.initial = 2,g1 = 4,beta.g1 = 1,beta.jmax =1){
  
 npp.df <- lai.drop.opt.func(VPD = VPD,
                              E = E,
                              PAR = PAR,
                              TMAX = TMAX,
                              lai.initial = lai.initial,
                              g1 = g1,
                              beta.g1 = beta.g1,
                              beta.jmax =beta.g1,
                              r.turnover = pars[1],
                              beta.lai = pars[2])
  
  return(-sum(npp.df$npp))
}



# func to use DEoptim to calaulate initial values
get.turnover.func <- function(...,
                              NPmax=100,
                              maxiter=200){
  # on.exit(stopCluster(cl = NULL))
  # setting control parameters and limits to values
  lower <- c(0.01,0.1)
  upper <- c(0.3,10)

  # 
  set.seed(1935)
  OptBB.de.fit <- DEoptim(fn=model.de.func,lower=lower,upper=upper,
                          DEoptim.control(VTR = -2000,
                                          NP = NPmax,itermax=maxiter,trace=1,parallelType = 1,
                                          parVar = list('lai.drop.opt.func')),
                          ...)
  # Sys.sleep(10)
  initial.vec <- unname(OptBB.de.fit$optim$bestmem)
  return(data.frame(r.turnover= initial.vec[1],
                    beta.lai = initial.vec[2],
                    npp = -OptBB.de.fit$optim$bestval))
}
out.df <- get.turnover.func(VPD = 1.5,E = 200,PAR = 8,TMAX = 30,
                  lai.initial = 2,g1 = 4,beta.g1 = 0.38,beta.jmax = 0.38)




