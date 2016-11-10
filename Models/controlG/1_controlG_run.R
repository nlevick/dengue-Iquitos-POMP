## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
    keep.source=TRUE,
    stringsAsFactors=FALSE,
    encoding="UTF-8",
    digits = 16
)

## ----prelims,echo=F,cache=F----------------------------------------------
set.seed(594709947L)
library(plyr)
library(reshape2)
library(magrittr)
library(readr)
library(pomp)
library(ggplot2)
theme_set(theme_bw())
stopifnot(packageVersion("pomp")>="1.4.5")

## ----load-external-data-----------------------------------------------------

dengue = read.csv(file = "Data/dengue.csv", header = TRUE, sep = ",")

control = read.csv(file = "Data/control_covar.csv", header = TRUE, sep = ",")

read_file("Models/controlG/controlG_rprocess.txt") %>%
    Csnippet() -> seir_rprocess

# change this line below to change data import
"Models/controlG/controlG" %>%
    paste(".rds", sep = "") %>%
    readRDS() -> bakedRDS
## ----dengue-data-----------------------------------------------------------

# import preprocessed dengue from 3_pomp_model
dengue = dengue[-c(1)] # remove first column

## ----covar---------------------------------------------------------
# Import 
control = control[-c(1)] # remove first column

## ----seir_names---------------------------------------------------------

seir_statenames <- c("C","S_H","E_H","I_H","R_H","N_H",
                     "J_G","S_G","E_G","I_G","N_G","noise",
                     "R_0","R_eff"
)

## ----start_params--------------------------------------------------------

seir_paramnames <- c("S_0","rho","psi","b0","b1","b2","eta","muG","epsilon")

# remove failed particles
bakedRDS %<>%
    subset(nfail.max==0) %>%
#         mutate(eta=exp(signif(log(eta),5))) %>%
#         ddply(~eta,subset,rank(-loglik)<=20)
    mutate(S_0=exp(signif(log(S_0),5))) %>%
    ddply(~S_0,subset,rank(-loglik)<=20)

# order from least negative loglik to most and take top line/best particle
bakedRDS = bakedRDS[order(bakedRDS$loglik, decreasing = TRUE) ,]
ordered = bakedRDS[1,]
# remove rownames for input into pomp
rownames(ordered) = NULL

# assign particle from rds (with integer pop values)
theta <- c(
    S_0=ceiling(ordered$S_0),
    # E_0=ceiling(ordered$E_0),
    # I_0=ceiling(ordered$I_0),
    rho=ordered$rho,
    psi=ordered$psi,
    b0=ordered$b0,
    b1=ordered$b1,
    b2=0,
    eta=ordered$eta,
    muG=ordered$muG,
    epsilon=ordered$epsilon
)
## ----seir-construct-------------------------------------------------------

seir_init <- Csnippet("
    // HUMAN
    S_H = nearbyint(S_0);
    E_H = 6;
    I_H = 2;
    R_H = 457864-nearbyint(S_0);
    C = 0;
    N_H = nearbyint(S_H + E_H + I_H + R_H);
    
    // VECTOR
    J_G = 3226269; 
    //J_G= 2599533; // reulermultinom model
    S_G = 800000;
    E_G = 0;
    I_G = 0;
    N_G = nearbyint(S_G + E_G + I_G);
    // phase = 0;
    noise = 0;
")


seir_dmeasure <- Csnippet("
    // Overdispersed binomial reporting
    double m = rho*C;
    double v = m*(1.0-rho+psi*psi*m);
    double tol = 1.0e-18;
    if (cases > 0.0) {
        lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
    } else {
        lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
    }
    
")

seir_rmeasure <- Csnippet("
    // Poisson
    // total_cases = rpois(rho*R+1e-6);
    
    // Overdispersed binomial reporting
    double m = rho*C;
    double v = m*(1.0-rho+psi*psi*m);
    double tol = 1.0e-18;
    cases = rnorm(m,sqrt(v)+tol);
    if (cases > 0.0) {
        cases = nearbyint(cases);
    } else {
        cases = 0.0;
    }
")

seir_toEstimationScale <- Csnippet("
    TS_0 = log(S_0);
    Trho = logit(rho);
    Tepsilon = log(epsilon);
    Tb0 = log(b0);
    Tb1 = logit(b1);
    Tb2 = logit(b2);
    Tpsi = log(psi);
    Teta = log(eta);
    TmuG = log(muG);
")

seir_fromEstimationScale <- Csnippet("
    TS_0 = exp(S_0);
    Trho = expit(rho);
    Tepsilon = exp(epsilon);
    Tb0 = exp(b0);
    Tb1 = expit(b1);
    Tb2 = expit(b2);
    Tpsi = exp(psi);
    Teta = exp(eta);
    TmuG = exp(muG);
")


# save.image(file = "pomp_model.RData")

## ----prof-design-------------------------------------------------

# set lower and upper bounds 
# theta is ( 1-S_0, 2-rho, 3-psi, 4-b0, 5-b1, 6-b2, 7-eta, 8-muG, 9-epsilon )

theta.lo = replace(theta, c(2,3,4,5,7,8,9), c(0,0,0, 0,0, 0,0))
theta.hi = replace(theta, c(2,3,4,5,7,8,9), c(1,1,1,1,10,10,1))

profileDesign(
    S_0=seq(from=0,to=350000,length=15),
    lower=theta.lo,upper=theta.hi,nprof=15
) -> pd

# pairs(~S_0+rho+psi+b0+b1+eta+muG+epsilon,data=pd)

## ----prof-round1,eval=TRUE,cache=FALSE---------------------------
library(foreach)
library(doParallel)

registerDoParallel(cores=8)

# change name of save file
bake("S0-profile13.rds",{
    foreach (p=iter(pd,"row"),
             .combine=rbind,
             .errorhandling="remove",
             .packages=c("pomp","magrittr","reshape2","plyr"),
             .inorder=FALSE,
             .options.mpi=list(chunkSize=1,seed=1598260027L,info=TRUE)
    ) %dopar% {
        
        tic <- Sys.time()
        
        options(stringsAsFactors=FALSE)
        
        dengue %>% 
            pomp(
                times="week",
                t0= with(dengue,2*week[1]-week[2]),
                rprocess=euler.sim(
                    step.fun=seir_rprocess,
                    delta.t=1/(365)
                ),
                rmeasure=seir_rmeasure,
                dmeasure=seir_dmeasure,
                fromEstimationScale=seir_fromEstimationScale,
                toEstimationScale=seir_toEstimationScale,
                covar=control,
                tcovar="time",
                zeronames=c("C","noise"),
                statenames=seir_statenames,
                paramnames=seir_paramnames,
                initializer=seir_init
            ) %>% 
            # Comment out profiled variables (will not be added to random walk)
            mif2(start = unlist(p),
                 Nmif = 50, 
                 rw.sd = rw.sd(
                    # S_0=ivp(0.02),
                    # E_0=ivp(0.02),
                    # I_0=ivp(0.02),
                    rho=0.02,psi=0.02,
                    b0=0.02,b1=0.02,
                    muG=0.02,
                    eta=0.02,
                    epsilon=0.02
                 ),
                 Np = 1000,
                 cooling.type = "geometric",
                 cooling.fraction.50 = 0.1,
                 transform = TRUE
            ) %>%
            mif2() -> mf
        
        ## Runs 10 particle filters to assess Monte Carlo error in likelihood
        pf <- replicate(10, pfilter(mf, Np = 2000))
        ll <- sapply(pf,logLik)
        ll <- logmeanexp(ll, se = TRUE)
        nfail <- sapply(pf,getElement,"nfail")
        
        toc <- Sys.time()
        etime <- toc-tic
        units(etime) <- "hours"
        
        data.frame(as.list(coef(mf)),
                   loglik = ll[1],
                   loglik.se = ll[2],
                   nfail.min = min(nfail),
                   nfail.max = max(nfail),
                   etime = as.numeric(etime))
    }
}) -> profile
