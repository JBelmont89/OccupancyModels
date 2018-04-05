#  Bayesian vs classical occupancy model
library(shiny)
#Bayesian and classical analysis
library(AHMbook)
library(jagsUI)

#For plotting
library(ggplot2)
library(ggmcmc)
library(gridExtra)
library(reshape2)
library(dplyr)


runApp(shinyApp(
  ui = fluidPage(
    titlePanel("Occupancy model simulation"),
    sidebarLayout(position = "left",
                  sidebarPanel(h4("Model parameters"),
                               withMathJax(),
                               HTML("This App fits the following two-state occupancy model with both classical and Bayesian methods: $$z_i\\sim \\text{Bernoulli}(\\psi)$$
                                                    $$y_i\\sim \\text{Binomial}(J,p)$$"),
                               HTML("Where \\( z_i \\) is the true occupancy state at M sites with occupacy probability\\( \\psi \\) and
                                    \\( z_i \\) are the observed occurences with detection probability \\( p \\) on the J surveys.
                                    For Bayesian analysis, a non-informative Uniform (0,1) prior is used for both \\( \\psi \\) and \\( p \\) parameters."),
                               sliderInput("psi", label= "Occupancy probability (\\( \\psi \\))", min = 0, max = 1, value = 0.5,step = 0.05),
                               sliderInput("p", label = "Detection probability (p)", min = 0, max = 1, value = 0.5,step = 0.05),
                               numericInput("M", label = "Sites",value = 100,min = 100,1000,step = 50),
                               numericInput("J", label = "Surveys",value = 3,min = 2,5,step = 1),
                               h4("Jags controls"),
                               radioButtons("iter",
                                            label = "Iterations",
                                            choices = list("1000" = 1000,"2000" = 2000,"5000" = 5000,"10000" = 10000),
                                            selected = 5000,inline=T),
                               radioButtons("thin",
                                            label = "Thinning",
                                            choices = list("1" = 1,"2" = 2,"5" = 5,"10" = 10),
                                            selected = 5,inline=T),
                               radioButtons("bn",
                                            label = "Burnin period",
                                            choices = list("1000" = 1000,"2000" = 2000,"3000" = 3000,"4000" = 4000),
                                            selected = 2000,inline=T)),
                  mainPanel(tabsetPanel(type = "tabs",
                                        tabPanel("Comparisson",plotOutput("main")),
                                        tabPanel("Trace plot",plotOutput("trace")),
                                        tabPanel("Autocorrelation Plot",plotOutput("autocorr")),
                                        tabPanel("Gelman Diagnostic Plot",plotOutput("gelman")),
                                        tabPanel("References",h2("References"),
                                                 p("- Kéry, M., & Royle, J. A. (2015).",em("Applied Hierarchical Modeling in Ecology: Analysis of distribution,
abundance and species richness in R and BUGS: Volume 1: Prelude and Static Models")," Academic Pres"),
                                                 p("- Royle, J. A., & Dorazio, R. M. (2008).", em("Hierarchical modeling and inference in ecology: the analysis
                                  of data from populations, metapopulations and communities"),". Elsevier."),
                                                 p("- Fernández-i-Marín, X. (2016). ggmcmc: Analysis of MCMC samples and Bayesian inference.",em("Journal of Statistical Software"),", 70(9)."),
                                                 p("- Wickham, H., & Wickham, M. H. (2007).", em("The ggplot package"),"."),
                                                 p("- Hornik, K., Leisch, F., & Zeileis, A. (2003).",em("JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling"),". In Proceedings of DSC (Vol. 2, pp. 1-1).")
                                                 ))
                  )
    )
  ),
  server = function(input, output) {
    Dataset <- reactive({
      set.seed(123)
      sim_dat1 <- simOcc(mean.occ=input$psi,beta1=0, beta2=0, beta3=0, mean.det=input$p,M = input$M,J = input$J,
                         time.effects=c(0, 0), alpha1=0, alpha2=0, alpha3=0, sd.lp=0, b=0,show.plot = F)
      sim.data <- list(y = sim_dat1$y, M = sim_dat1$M, J = sim_dat1$J)
      cat("
          model {
          # Priors
          psi ~ dunif(0, 1)
          p ~ dunif(0, 1)
          # Likelihood
          for (i in 1:M) { # Loop over sites
          z[i] ~ dbern(psi) # State model
          for (j in 1:J) { # Loop over replicate surveys
          y[i,j] ~ dbern(z[i]*p) # Observation model 
          }
          }
          }
          ", fill = TRUE,file="sim_occ.txt")
      zst <- apply(sim_dat1$y, 1, max) 
      inits <- function(){list(z = zst)}
      params <- c("psi", "p")
      withProgress(message = "Application loading", value = 0, {
        incProgress(0.2, detail = "Burn-in phase")
        Sys.sleep(0.5)
        incProgress(0.5, detail = "Sampling from joint posterior")
        OccuJAGs <- jags(sim.data, inits, params, "sim_occ.txt", n.chains = 3,
                         n.thin = as.numeric(input$thin), n.iter = as.numeric(input$iter), n.burnin = as.numeric(input$bn))
        incProgress(0.3, detail = "Done")
      })
      umf <- unmarkedFrameOccu(y = sim_dat1$y) 
      fm1 <- occu(~1 ~1, data = umf)
      list(sim_dat1 = sim_dat1,OccuJAGs=OccuJAGs,fm1=fm1)
    })
    
    output$main <- renderPlot({
      
      sim_dat1 <- Dataset()$sim_dat1
      OccuJAGs <- Dataset()$OccuJAGs
      fm1 <- Dataset()$fm1
      plot.df <- data.frame(probs=c(input$psi,OccuJAGs$mean$psi,backTransform(fm1, "state")@estimate,
                                    input$p,OccuJAGs$mean$p,backTransform(fm1, "det")@estimate),
                            method= rep(c("True","Bayesian","Classical"),2),
                            low = c(NA,OccuJAGs$q2.5$psi, confint(backTransform(fm1, "state"))[1],
                                    NA, OccuJAGs$q2.5$p,plogis(confint(fm1, type='det'))[1]),
                            high = c(NA,OccuJAGs$q97.5$psi,confint(backTransform(fm1, "state"))[2],
                                     NA, OccuJAGs$q97.5$p, plogis(confint(fm1, type='det'))[2]),
                            parameters= rep(c("psi","p"),each=3)
      )
      
      
      ggplot(plot.df,aes(x=parameters,y=probs,colour=method,group=method))+
        ylim(0,1)+
        geom_errorbar(aes(x=parameters, ymin=low, ymax=high),width = 0.25,alpha=0.45,na.rm = T,
                      position = position_dodge(width = .06))+
        geom_point(aes(shape=method),position=position_dodge(width = .06),size=2)+
        labs(y="Detection/Occupancy probabilities",x="Parameters")+
        scale_colour_manual(name="",values=c("purple","orange","tomato"),guide = FALSE)+
        scale_shape_manual(name="",values=c(16,16,17))+
        scale_x_discrete(labels=c("p"="p","psi"=expression(psi)))+
        guides(shape = guide_legend(override.aes = 
                                      list(size=2,shape=c(16,16,17), 
                                           colour=c("purple","orange","tomato") ))) +
        theme(text = element_text(size=12),
              axis.text.x = element_text(hjust=1,size = 14))
      
    },height = 600)
    output$trace <- renderPlot({
      OccuJAGs <- Dataset()$OccuJAGs
      S <- ggs(OccuJAGs$samples)
      f1 <- ggs_density(S,greek=T)
      f2 <- ggs_traceplot(S,greek=T)+theme(legend.position = 0)
      grid.arrange(f2, f1, ncol=2, nrow=1)
    })
    output$autocorr <- renderPlot({
      OccuJAGs <- Dataset()$OccuJAGs
      ggs_autocorrelation(S,greek=T)
    })
    output$gelman <- renderPlot({
      OccuJAGs <- Dataset()$OccuJAGs
      gp.dat = gelman.plot(OccuJAGs$samples)
      df = data.frame(bind_rows(as.data.frame(gp.dat[["shrink"]][,,1]), 
                                as.data.frame(gp.dat[["shrink"]][,,2])), 
                      q=rep(dimnames(gp.dat[["shrink"]])[[3]], each=nrow(gp.dat[["shrink"]][,,1])),
                      last.iter=rep(gp.dat[["last.iter"]], length(gp.dat)))
      ggplot(melt(df, c("q","last.iter")), 
             aes(last.iter, value, colour=q, linetype=q)) + 
        geom_hline(yintercept=1, colour="grey30", lwd=0.2) +
        geom_line() +
        facet_wrap(~variable, labeller= labeller(.cols=function(x) gsub("V", "Chain ", x))) +
        labs(x="Last Iteration in Chain", y="Shrink Factor",
             colour="Quantile", linetype="Quantile") +
        scale_linetype_manual(values=c(2,1))
    })
    }
    ),launch.browser = T)







