library(magrittr)
library(reshape2)
library(ggplot2)
library(plyr)
library(latex2exp)
source("multiplot.R")
# ----load-mif_run-------------------------------------------------------

plot = TRUE
pdf = FALSE # pdf: T png: F

rdsName = "Models/control2b0/control2b0" #ControlG

rdsName%>%
    paste(".rds", sep = "") %>%
    readRDS() -> bakedRDS

# # run statistics 
sum(bakedRDS$etime)/8
summary(bakedRDS)

# # ----clean-mif_run-------------------------------------------------------
# # Remove failed particles

bakedRDS %<>%
    subset(nfail.max==0) %>%
    mutate(S_0=exp(signif(log(S_0),5))) %>%
    ddply(~S_0,subset,rank(-loglik)<=20)

# order from least negative loglik to most and take top line/best particle
bakedRDS = bakedRDS[order(bakedRDS$loglik, decreasing = TRUE) ,]
ordered = bakedRDS[1,]
# remove rownames for input into pomp
rownames(ordered) = NULL


# # ----plot-mif_run-------------------------------------------------------

# won't work without previous theta -- look up in associate 3_*.r file
for (i in 1:length(theta)){ # all parameters up to epsilon
    
    if (!plot){ # to save files
        if (pdf){
            pdf( # create a file for each plot by concat strings
                paste(
                    "Plots/",
                    paste(
                        paste(
                            rdsName,
                            names(theta)[i],
                            sep="-"
                        ),
                        ".pdf",
                        sep=""
                    ),
                    sep=""
                )
            )
        } 
        else { #png
            png( # create a file for each plot by concat strings
                paste(
                    "Plots/",
                    paste(
                        # paste(
                            names(theta)[i],
                            # rdsName,
                            # names(theta)[i],
                            # sep="-"
                         # ),
                        ".png",
                        sep=""
                    ),
                    sep=""
                ),
                width = 870, height = 600, units = "px"
            ) 
        }
        
        
        
        
        switch((names(theta)[i]), 
               S_0={
                   xAxis = TeX('Initial host susceptible pop. ($S_{0}$)')
               },
               rho={
                   xAxis = TeX('Reporting rate ($\\rho$)')  
               },
               psi={
                   xAxis = TeX('Overdispersion factor ($\\psi$)')  
               },
               b0={
                   xAxis = TeX('Transmission average value ($b_{0}$)')  
               },
               b0H={
                   xAxis = TeX('Host transmission average ($b_{0}^{H}$)')  
               },
               b0G={
                   xAxis = TeX('Vector transmission average ($b_{0}^{G}$)')  
               },
               b1={
                   xAxis = TeX('Transmission amplitude ($b_{1}$)')  
               },
               c1={
                   xAxis = TeX('EIP scale parameter ($c_{1}$)')  
               },
               c2={
                   xAxis = TeX('EIP shape parameter ($c_{2}$)')  
               },
               muG={
                   xAxis = TeX('Avg. adult vector lifespan $(\\hat{\\mu}_{G})^{-1}$') 
               },
               eta={
                   xAxis = TeX('Infected host importation ($\\eta$)')  
               },
               epsilon={
                   xAxis = TeX('Gamma white noise intensity ($\\epsilon$)')  
               }
        )
        
        
        assign(
            # p,
            paste0("p",i),
            bakedRDS %>%
                ggplot(aes_string(x=(names(theta)[i]),y="loglik")) +
                labs(x=xAxis,y="Profile Log-Likelihood") +
                scale_x_continuous(labels = scales::dollar_format("")) +
                theme(
                    axis.title.x = element_text(size=34, vjust=-0.85),
                    axis.title.y = element_text(size=34, vjust=0.35),
                    axis.text.x = element_text(size=27, vjust=-0.85),
                    axis.text.y = element_text(size=27, vjust=0.35),
                    plot.margin = unit(c(1,2,1,1), "cm")
                    
                ) +
                geom_point()+
                geom_smooth(method="loess")
        ) %>%
        plot()
        dev.off()
    }
        
    else{ # show plot, no save
        
        xAxis = ""
        
        switch((names(theta)[i]), 
               S_0={
                    xAxis = TeX('Initial host susceptible pop. ($S_{0}$)')
               },
               rho={
                   xAxis = TeX('Reporting rate ($\\rho$)')  
               },
               psi={
                   xAxis = TeX('Overdispersion factor ($\\psi$)')  
               },
               b0={
                   xAxis = TeX('Transmission average value ($b_{0}$)')  
               },
               b0H={
                   xAxis = TeX('Host transmission average ($b_{0}^{H}$)')  
               },
               b0G={
                   xAxis = TeX('Vector transmission average ($b_{0}^{G}$)')  
               },
               b1={
                   xAxis = TeX('Transmission amplitude ($b_{1}$)')  
               },
               c1={
                   xAxis = TeX('EIP scale parameter ($c_{1}$)')  
               },
               c2={
                   xAxis = TeX('EIP shape parameter ($c_{2}$)')  
               },
               muG={
                   xAxis = TeX('Control weighted parameter $(\\hat{\\mu}_{G})^{-1}$') 
               },
               eta={
                   xAxis = TeX('Infected host importation ($\\eta$)')  
               },
               epsilon={
                   xAxis = TeX('Gamma white noise intensity ($\\epsilon$)')  
               }
        )
        if (substr(names(theta)[i],1,1)=='c' ) {
            labels = scale_x_continuous()
            xAxisTextSize = 10
            angle = 75
        } else {
            labels = scale_x_continuous(labels = scales::dollar_format(""))
            xAxisTextSize = 14
            angle = 0
        }
        assign(
            paste0("p",i),
            bakedRDS %>%
                ggplot(aes_string(x=(names(theta)[i]),y="loglik")) +
                labs(x=xAxis,y="") +
                labels +
                theme(
                    axis.title.x = element_text(size=18, vjust=-0.85),
                    axis.text.x = element_text(size=xAxisTextSize, angle = angle),
                    axis.text.y = element_text(size=16, vjust=0.35)
                    
                ) +
                geom_point()+
                geom_smooth(method="loess")
        ) 
        

    }

    
    

}

multiplot(p1,p2,p3,p4,p5,p6,p8,p9,p10, # control2b0
          layout = matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow=TRUE),
          labs=list("Parameter", "Profile Log-Likelihood"),
          labpos=list(c(0.5,-0.03), c(0.01,0.5))
)

#------------------------#------------------------#------------------------#------------------------
# Formatting for scientific notation
fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
#     if (parse(text=l) == '0.0e+00'){
#         l <- 0
#     }
    parse(text=l)
}