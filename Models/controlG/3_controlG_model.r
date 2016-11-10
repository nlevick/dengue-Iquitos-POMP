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
library(latex2exp)
theme_set(theme_bw())
stopifnot(packageVersion("pomp")>="1.4.5")
source("multiplot.R")
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

## ----names---------------------------------------------------------

seir_statenames <- c("C","S_H","E_H","I_H","R_H","N_H",
                     "J_G","S_G","E_G","I_G","N_G",
                     "noise","R_0","R_eff"
)

## ----start_params--------------------------------------------------------

seir_paramnames <- c("S_0","rho","psi","b0","b1","b2","eta","muG","epsilon")

# remove failed particles
bakedRDS %<>%
    subset(nfail.max==0) %>%
    mutate(S_0=exp(signif(log(S_0),5))) %>%
    # change number of best particles for each profile value
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

## ----construct-------------------------------------------------------


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
    //phase = 0;
    noise = 0;
")


seir_dmeasure <- Csnippet("
    // Poisson
    // lik = dpois(total_cases,rho*R+1e-6,give_log);
    
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


## ----pomp----------------------------------------------------------

m1 <- pomp(
    data=dengue,
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
)

# # ----init_sim------------------------------------------------------------

plot = TRUE

if (plot){

    #--------------------------------------------------------------------------
    # Figure control_sims
    #--------------------------------------------------------------------------
    ## PLOT SIMULATION WITH DATA
    sims <- simulate(m1,params=theta,nsim=3,as=TRUE,include=TRUE)
    p1 = ggplot(sims, mapping = aes(x=time, y=cases, group=sim,
            # color=sim=="data" # runs vs data color scheme
            color=factor(sim)   # different color each run
            )
        ) + 
        labs(x = "", y = "Reported cases (c)") +         
        theme(
            axis.title.y = element_text(size=34, vjust=0.35),
            axis.text.x = element_text(size=24, vjust=-0.85),
            axis.text.y = element_text(size=24, vjust=0.35)
        ) +
        geom_line() + guides(color=FALSE)
    # print(p1)
    
    #--------------------------------------------------------------------------
    # Figure control_R
    #--------------------------------------------------------------------------
    
    p3 = ggplot(sims, mapping = aes(x=time,y=R_0,group=sim, color=factor(sim))) +
        # labs(x = "", y = TeX('Basic Reproductive No. ($R_{0}$)')) +
        labs(x = "", y = TeX('$R_{0}$')) + 
        theme(
            # axis.title.y = element_text(size=24, vjust=0.35),
            axis.title.y = element_text(size=34, vjust=0.35),
            axis.text.x = element_text(size=24, vjust=-0.85),
            axis.text.y = element_text(size=24, vjust=0.35),
            plot.margin = unit(c(0.5,0.5,0,.5), "cm")
        ) +
        geom_line() + guides(color=FALSE)
    cutoff <- data.frame(yintercept=1, cutoff=factor(1))
    p3 = p3+geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff, show.legend=FALSE) 
    # print(p3)
    
    p2 = ggplot(sims, mapping = aes(x=time,y=R_eff,group=sim, color=factor(sim))) +
        # labs(x = "", y = TeX('Eff. Reproductive No. ($R_{f}$)')) + 
        labs(x = "", y = TeX('$R_{f}$')) + 
        theme(
            # axis.title.y = element_text(size=24, vjust=0.35),
            axis.title.y = element_text(size=34, vjust=0.35),
            axis.text.x = element_text(size=24, vjust=-0.85),
            axis.text.y = element_text(size=24, vjust=0.35),
            plot.margin = unit(c(.5,.5,0,.5), "cm")
        ) +
        # scale_y_continuous(labels = scales::dollar_format("")) +
        geom_line() + guides(color=FALSE)
    
    cutoff <- data.frame(yintercept=1, cutoff=factor(1))
    p2 = p2+geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff, show.legend=FALSE) 
    # print(p2)
    
    multiplot(p3,p2,
              cols = 1
              # ,labpos=list(c(0.5,0.03), c(0.03,0.5))
    )
    

    #--------------------------------------------------------------------------
    # Internal state plots
    #--------------------------------------------------------------------------
    
    pl0 = ggplot(sims, mapping = aes(x=time,y=C,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$C$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)
        ) +
        geom_line() + guides(color=FALSE)
    pl1 = ggplot(sims, mapping = aes(x=time,y=S_H,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$S_{H}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)
        ) +
        geom_line() + guides(color=FALSE)
    pl2 = ggplot(sims, mapping = aes(x=time,y=E_H,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$E_{H}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)
        ) +
        geom_line() + guides(color=FALSE)
    pl3 = ggplot(sims, mapping = aes(x=time,y=I_H,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$I_{H}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)
        ) +
        geom_line() + guides(color=FALSE)    
    pl4 = ggplot(sims, mapping = aes(x=time,y=R_H,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$R_{H}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)  
        ) +
        geom_line() + guides(color=FALSE)
    pl5 = ggplot(sims, mapping = aes(x=time,y=N_H,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$N_{H}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)
        ) +
        geom_line() + guides(color=FALSE)
    pl6 = ggplot(sims, mapping = aes(x=time,y=J_G,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$J_{G}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)  
        ) +
        geom_line() + guides(color=FALSE)
    pl7 = ggplot(sims, mapping = aes(x=time,y=S_G,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$S_{G}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18) 
        ) +
        geom_line() + guides(color=FALSE)
    pl8 = ggplot(sims, mapping = aes(x=time,y=E_G,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$E_{G}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)  
        ) +
        geom_line() + guides(color=FALSE)
    pl9 = ggplot(sims, mapping = aes(x=time,y=I_G,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$I_{G}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)  
        ) +
        geom_line() + guides(color=FALSE)
    pl10 = ggplot(sims, mapping = aes(x=time,y=N_G,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$N_{G}$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)  
        ) +         
        geom_line() + guides(color=FALSE)
    pl11 = ggplot(sims, mapping = aes(x=time,y=noise,group=sim, color=factor(sim))) +
        labs(x = "", y = TeX('$\\zeta$')) + 
        theme(
            axis.title.y = element_text(size=28, vjust=0.35),
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18)
        ) +
        scale_y_continuous(labels = scales::dollar_format("")) +
        geom_line() + guides(color=FALSE)
    cutoff <- data.frame(yintercept=0, cutoff=factor(0))
    pl11 = pl11+geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff, show.legend=FALSE) 

#     # To print all internal states on one plot
#     multiplot(pl0,pl1,pl2,pl3,pl4,pl5,pl6,pl7,pl8,pl9,pl10,
#               cols = 3,
#               # layout = matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow=TRUE),
#               labs=list("Year", ""),
#               labpos=list(c(0.5,0.02), c(0.01,0.5))
#     )
    
    #--------------------------------------------------------------------------
    # Figure SOb0
    #--------------------------------------------------------------------------

    p3 = ggplot(bakedRDS,aes(x=S_0, y=b0)) + 
        labs(x = TeX('Initial Human Susceptible Population ($S_{0}$)'), y = TeX('Transmission average value ($b_{0})$')) + 
        theme(
            axis.title.x = element_text(size=34, vjust=0.35),
            axis.title.y = element_text(size=34, vjust=0.35),
            axis.text.x = element_text(size=30, vjust=-0.85),
            axis.text.y = element_text(size=30, vjust=0.35)  
        ) +
        scale_y_continuous(labels = scales::dollar_format("")) +
        scale_x_continuous(labels = scales::dollar_format("")) +
        geom_point()
        # print(p3)
}





