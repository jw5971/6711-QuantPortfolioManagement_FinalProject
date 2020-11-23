library(shiny)
library(readr)
library(quadprog)
library(quantmod)
library(tseries)
library(zoo)

################################1. Data Sourcing##############################
st_art <- "2007-03-23" #earliest time that all ETFs data are available
en_d <- "2020-10-01" # end date 2020-09-30
quo_te <- "Close"

# load ETFs data 
Inv_universe <- c("fxe",
                  "ewj",
                  "gld",
                  "qqq",
                  "spy",
                  "shv",
                  "dba",
                  "uso",
                  "xbi",
                  "ilf",
                  "gaf",
                  "epp",
                  "fez")

ETFs <- lapply(Inv_universe,function(x) get.hist.quote(instrument = x, start = st_art, end = en_d, quote = quo_te))
ETFs_DF <- do.call(cbind.data.frame, ETFs)
names(ETFs_DF) <- Inv_universe
ETFs_DF
#Read data from drive
setwd("C:/Users/JIWU5/OneDrive/Desktop/quantitative-portfolio-management-master/Project")
fama3 <- read_csv("./fama_3.csv")
fama3_DF <- do.call(data.frame, fama3) 
rownames(fama3_DF) <-fama3_DF$Date #set index
fama3_DF <- fama3_DF[,-1] # drop date
fama3_DF <- (fama3_DF+1)^(1/365)-1 #annulized to daily
All_Data <- merge(ETFs_DF, fama3_DF, by=0, all=TRUE) #merge by row names (by=0 or by="row.names")
rownames(All_Data) <- All_Data$Row.names
All_Data <- All_Data[,-1]
head(All_Data)
################################# 2. Data Processing ##########################
return_ETFS <- as.data.frame(diff(as.matrix(log(ETFs_DF))))
All_Data_Ret <- merge(return_ETFS, fama3_DF, by=0, all=TRUE)
All_Data_Ret <- na.omit(All_Data_Ret)
return_ETFSlessRF <- return_ETFS - All_Data_Ret$RF

All_Data_RetlessRF <- merge(return_ETFSlessRF, fama3_DF, by=0, all=TRUE)
All_Data_RetlessRF <- na.omit(All_Data_RetlessRF)

################################################################################
ETFs_fit_func <- function(dataset,betaset)
{
  
  rho_ETFs <- lapply(Inv_universe,function(x) 
  {
    betaset[x,'Mkt.RF']*dataset$Mkt.RF+
      betaset[x,'SMB']*dataset$SMB+
      betaset[x,'HML']*dataset$HML+
      betaset[x,'(Intercept)']
  })
  rho_ETFs_set <- do.call(cbind, rho_ETFs)
  colnames(rho_ETFs_set) <- Inv_universe
  return(rho_ETFs_set)
}

#tail(All_Data_RetlessRF)
###############################################################################
divide_data <- function(per_iods,crisis_day) {
  short <- per_iods
  crisis <- which(All_Data_RetlessRF$Row.names == crisis_day)
  
  bf_c1 <- crisis - 1.5*short
  bf_c2 <- crisis - 0.5*short
  med_c1 <- crisis - 0.5*short
  med_c2 <- crisis + short*0.5
  af_c1 <- crisis + short*0.5
  af_c2 <- crisis + short*1.5
  
  
  Bf_c <- All_Data_RetlessRF[bf_c1:bf_c2,]
  #head(Bf_short)
  Med_c <- All_Data_RetlessRF[med_c1:med_c2,]
  Aft_c <- All_Data_RetlessRF[af_c1:af_c2,]
  #########################OLS Regression for beta #########################
  var_x <- "Mkt.RF + SMB + HML"
  fomula_string <- paste(Inv_universe,var_x, sep = "~")
  names(fomula_string) <- Inv_universe
  
  d <- list()
  d$Bf_c <- Bf_c
  d$Med_c <- Med_c
  d$Aft_c <- Aft_c
  
  regre_all_result <- lapply(d, function(y) 
  {
    regre_ssion <- lapply(Inv_universe,
                          function(x) lm(as.formula(fomula_string[x]),data = y)$coeff)
    regre_result <- do.call(rbind, regre_ssion)
    rownames(regre_result) <- Inv_universe
    return(regre_result)
  } )
  
  rho_ETFs_set = list()
  
  rho_ETFs_set$Bf_c <- ETFs_fit_func(Bf_c,regre_all_result$Bf_c)
  rho_ETFs_set$Med_c <- ETFs_fit_func(Med_c,regre_all_result$Med_c)
  rho_ETFs_set$Aft_c <- ETFs_fit_func(Aft_c,regre_all_result$Aft_c)
  l <- list()
  l$rho_ETFs_set <- rho_ETFs_set
  l$regre_all_result <- regre_all_result
  l$asset <- d
  return(l)
}

################################################################################

opt_Tar_ret <- function(ma_t,TargetRet) 
{
  nAssets <- dim(ma_t)[2]
  upperB <- 2
  lowerB <- -2
  ub <- rep(upperB, nAssets)
  lb <- rep(lowerB, nAssets)
  cov.mat <- cov(ma_t, use = "complete.obs")
  diag(cov.mat) <- 1
  
  Amat <- as.matrix(cbind(rep(1, times=13),colMeans(ma_t),diag(nAssets), -diag(nAssets)))
  #TargetRet <- 0.15
  TargetReturn = (TargetRet+1)**(1/365)-1
  bvec <- c(1, TargetReturn,lb,-ub)
  tar_ret_result <- solve.QP(Dmat = cov.mat, dvec = rep(0,13), Amat = Amat, bvec = bvec, meq = 2)
  tar_ret_result$solution
}

#######################################################################

opt_beta <- function(rho_mat,beta_mat,Beta)
{
  nAssets <- dim(rho_mat)[2]
  upperB <- 2
  lowerB <- -2
  ub <- rep(upperB, nAssets)
  lb <- rep(lowerB, nAssets)
  cov.mat <- cov(rho_mat, use = "complete.obs")
  diag(cov.mat) <- 1
  
  Amat <- as.matrix(cbind(rep(1, times=13),beta_mat[,"Mkt.RF"],diag(nAssets), -diag(nAssets)))
  #Beta <- 0.15 #0.5/1/1.5
  TarBeta <- (Beta+1)**(1/365)-1
  bvec <- c(1, TarBeta,lb, -ub)
  result <- solve.QP(Dmat = cov.mat, dvec = colMeans(rho_mat), Amat = Amat, bvec = bvec, meq = 2)
  names(result$solution) <- Inv_universe
  result$solution
  
}

Beta_w <- function(Term,Beta,crisis_day){
  opt_beta_w <- list()
  #cal_cu <- divide_data(Term)
  cal_cu <- divide_data(Term,crisis_day)
  opt_beta_w$Bf_c <- opt_beta(cal_cu$rho_ETFs_set$Bf_c,cal_cu$regre_all_result$Bf_c,Beta = Beta)
  opt_beta_w$Med_c <- opt_beta(cal_cu$rho_ETFs_set$Med_c,cal_cu$regre_all_result$Med_c,Beta = Beta)
  opt_beta_w$Aft_c <- opt_beta(cal_cu$rho_ETFs_set$Aft_c,cal_cu$regre_all_result$Aft_c,Beta = Beta)
  opt_beta_w
}

Ret_w <- function(Term,targetRet,crisis_day){
  rho_ETFs_set <- divide_data(Term,crisis_day)$rho_ETFs_set
  opt_ret_w <- lapply(rho_ETFs_set,function(x) {opt_Tar_ret(x,targetRet)})
  opt_ret_w
}

perf_set<- function(weights,return)
{
  S <- t(matrix(weights,ncol = 13) %*% t(return[,2:14]))
  #cumprod(1+)
  S <- as.data.frame(S)
  rownames(S) <- return$Row.names
  
  cum_ret <- prod(1+S)
  #cumprod(1+S)
  Mean_Geo_ret <- prod(1+S)^(1/dim(S)[1])-1
  Min_ret <- min(S)
  MXDD_10 <- max(rollapply(S,10,PerformanceAnalytics::maxDrawdown))
  S_td <- PerformanceAnalytics::StdDev(S)
  shar_pe <- colMeans(S - return$RF)/S_td
  sk_ew <- PerformanceAnalytics:: skewness(S)
  kurt <- PerformanceAnalytics::kurtosis(S)
  V_aR <- PerformanceAnalytics:: VaR(S)
  C_VaR <- PerformanceAnalytics::ETL(S)
  perf_mat <- c(cum_ret,Mean_Geo_ret,Min_ret,MXDD_10,S_td,shar_pe,sk_ew,kurt,V_aR,C_VaR)
  names(perf_mat) <- c("Cumulative_Return","Mean_Geometric_Return","Min_Return","10_days_MaxDrawdown","Volatility","Sharpe_Ratio","skewness","kurtosis","VaR","CVaR")
  t(perf_mat)
}

Ret_performance <- function(weights,return)
{
  return_performance_set <- list()
  return_performance_set$Bf_c <- perf_set(weights$Bf_c,return$Bf_c)
  return_performance_set$Med_c <- perf_set(weights$Med_c,return$Med_c)
  return_performance_set$Aft_c <- perf_set(weights$Aft_c,return$Aft_c)
  return_performance_table <- do.call(rbind,return_performance_set)
  return_performance_table
}

Beta_performance <- function(weights,return)
{
  Beta_performance_set <- list()
  Beta_performance_set$Bf_c <- perf_set(weights$Bf_c,return$Bf_c)
  Beta_performance_set$Med_c <- perf_set(weights$Med_c,return$Med_c)
  Beta_performance_set$Aft_c <- perf_set(weights$Aft_c,return$Aft_c)
  Beta_performance_table <- do.call(rbind,Beta_performance_set)
  Beta_performance_table
}

###############################################
# Define UI for Market Risk application
inter_face <- shiny::pageWithSidebar(
  
  # Application title.
  headerPanel("Long/Short Market Risk Portfolio Performance"),
  
  sidebarPanel(
    selectInput("Method", "Choose a method:", 
                choices = c("Target Return", "Target Beta")),
    
    selectInput("Crisis", "Choose a crisis date:", 
                choices = c("2008-09-16", "2010-04-28","2014-12-01")),
    
    column(width=12, sliderInput("Term", label="Time window:",
                                 min=30, max=100, value=60, step=10)),
    column(width=12, sliderInput("TargetRet", label="Target Return:",
                                 min=0.0, max=0.5, value=0.15, step=0.05)),
    column(width=12, sliderInput("Beta", label="Target Beta:",
                                 min=-1, max=2, value=1, step=0.5)),
    
    helpText("Note: The data unit is percentage(%)."),
    
    submitButton("Update View")
  ),
  
  mainPanel(
    
    h4("Weights"),
    tableOutput("view"),
    h4("Before Crisis"),
    plotOutput("Bf_c"),
    h4("In Crisis"),
    plotOutput("Med_c"),
    h4("After Crisis"),
    plotOutput("Aft_c"),
    h4("Performance"),
    tableOutput("Performance")
    
  )
)

#######################################################################
ser_ver <-function(input, output) {
  
  # Return the requested dataset
  datasetInput <- reactive({
  
    switch(input$Method,
             "Target Return" =  list(Ret_w(input$Term,input$TargetRet,input$Crisis),
                                     Ret_performance(Ret_w(input$Term,input$TargetRet,input$Crisis),
                                                     divide_data(input$Term,input$Crisis)$asset)),
             "Target Beta" = list(Beta_w(input$Term,input$Beta,input$Crisis),
                                  Beta_performance(Beta_w(input$Term,input$Beta,input$Crisis),
                                                   divide_data(input$Term,input$Crisis)$asset)
)
           )
  })

  output$view <- renderTable({
    da_ta <- as.matrix(do.call(rbind,datasetInput()[[1]]))
    colnames(da_ta) <- Inv_universe
    rownames(da_ta) <- c("Before_crisis","In_Crisis","After_Crisis")
    da_ta
    
  },colnames = TRUE,rownames = TRUE, digits = 3)
  
  output$Aft_c <- renderPlot({
    barplot(datasetInput()[[1]]$Aft_c, names.arg = Inv_universe)
    
  })
  
  output$Med_c <- renderPlot({
    barplot(datasetInput()[[1]]$Med_c, names.arg = Inv_universe)
  })  
  
  output$Bf_c <- renderPlot({
    
    barplot(datasetInput()[[1]]$Bf_c, names.arg = Inv_universe)
  })
  
  output$Performance <- renderTable({
    
    out <- as.matrix(t(datasetInput()[[2]]))
    colnames(out) <- c("Before Crisis","In Crisis","After Crisis")
    out
  },rownames = TRUE,colnames = TRUE, digits = 3)
  
}

shiny::shinyApp(ui=inter_face, server=ser_ver)