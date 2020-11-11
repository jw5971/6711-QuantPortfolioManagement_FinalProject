#####################Load package
library(quadprog)
library(tseries)
library(readr)
library(zoo)
library(quantmod)

# install.packages(c("readr"))

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
# head(fama3)
# tail(fama3)

fama3_DF <- do.call(data.frame, fama3) 
rownames(fama3_DF) <-fama3_DF$Date #set index
fama3_DF <- fama3_DF[,-1] # drop date
fama3_DF <- (fama3_DF+1)^(1/365)-1 #annulized to daily
# head(fama3_DF)
# head(rownames(fama3_DF))
# merge two datasets
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
head(All_Data_RetlessRF)

write.csv(All_Data_RetlessRF, file = "All_Data_RetlessRF.csv")
#####################  separate the dataset  ####################
#September 15, 2008: Lehman Brothers Bankruptcy Triggered Global Panic.
short <- 100
#60
# long <- 200
# med <- 110
#before crisis - short_period/long_period
#  crisis - sho
crisis <- which(All_Data_RetlessRF$Row.names == "2008-09-16")

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
#################################### 3. Regression #########################
#OLS model
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

regre_all_result

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

rho_ETFs_set = list()

rho_ETFs_set$Bf_c <- ETFs_fit_func(Bf_c,regre_all_result$Bf_c)
rho_ETFs_set$Med_c <- ETFs_fit_func(Med_c,regre_all_result$Med_c)
rho_ETFs_set$Aft_c <- ETFs_fit_func(Aft_c,regre_all_result$Aft_c)

#rownames(rho_ETFs_set$Aft_long) <- Aft_long[,"Row.names"]
# output
#write.csv(rho_ETFs_set$Aft_long, file = "test_data.csv")
#write.csv(regre_result, file = "regre_result.csv")
########################### 4. optimization Target return #################

opt_Tar_ret <- function(ma_t) 
  {
    nAssets <- dim(ma_t)[2]
    upperB <- 2
    lowerB <- -2
    ub <- rep(upperB, nAssets)
    lb <- rep(lowerB, nAssets)
    cov.mat <- cov(ma_t, use = "complete.obs")
    diag(cov.mat) <- 1

    Amat <- as.matrix(cbind(rep(1, times=13),colMeans(ma_t),diag(nAssets), -diag(nAssets)))
    # TargetReturn <- 0.15
    TargetReturn = (0.15+1)**(1/365)-1
    bvec <- c(1, TargetReturn,lb,-ub)
    tar_ret_result <- solve.QP(Dmat = cov.mat, dvec = rep(0,13), Amat = Amat, bvec = bvec, meq = 2)
    tar_ret_result$solution
    
    #verify:
    #sum(tar_ret_result$solution* colMeans(ma_t))
    #vol:
    #t(tar_ret_result$solution) %*% cov.mat %*% tar_ret_result$solution
}

opt_ret_w <- lapply(rho_ETFs_set,opt_Tar_ret)

################################# 5. optimization Target beta ################
# beta_mat <- regre_all_result
# rho_mat <- rho_ETFs_set

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
#sum(result$solution* beta_mat$Bf_short[,"Mkt.RF"])
opt_beta_w <- list()
Be_ta <- 0.5
#
opt_beta_w$Bf_c <- opt_beta(rho_ETFs_set$Bf_c,regre_all_result$Bf_c,Beta = Be_ta)
opt_beta_w$Med_c <- opt_beta(rho_ETFs_set$Med_c,regre_all_result$Med_c,Beta = Be_ta)
opt_beta_w$Aft_c <- opt_beta(rho_ETFs_set$Aft_c,regre_all_result$Aft_c,Beta = Be_ta)

################################### 6. portfolio analysis ######################

#opt return (annual return 15%)

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

return_performance_set <- list()
return_performance_set$Bf_c <- perf_set(opt_ret_w$Bf_c,Bf_c)
return_performance_set$Med_c <- perf_set(opt_ret_w$Med_c,Med_c)
return_performance_set$Aft_c <- perf_set(opt_ret_w$Aft_c,Aft_c)

#rbind(performance_set)
return_performance_table <- do.call(rbind,return_performance_set)
rownames(return_performance_table) <- c("Bf_c","Med_c","Aft_c")

write.csv(return_performance_table, file = "return_performance_table.csv")

#*(253)^(1/2)  annualized
#opt return (annual target beta 15%)
#opt_beta_w

Beta_performance_set <- list()
Beta_performance_set$Bf_c <- perf_set(opt_beta_w$Bf_c,Bf_c)
Beta_performance_set$Med_c <- perf_set(opt_beta_w$Med_c,Med_c)
Beta_performance_set$Aft_c <- perf_set(opt_beta_w$Aft_c,Aft_c)

#rbind(performance_set)
Beta_performance_table <- do.call(rbind,Beta_performance_set)
rownames(Beta_performance_table) <- c("Bf_c","Med_c","Aft_c")


write.csv(Beta_performance_table, file = "Beta_performance_table_05.csv")

######################## 7.security analysis ###################################

perf_set_security<- function(date,return,Rf)
{

  S <- as.data.frame(return)
  rownames(S) <- date
  #cumprod(1+)
  #rownames(S) <- t(date)
  cum_ret <- apply(1+S,2,prod)
  #cumprod(1+S) prod(1+S)
  Mean_Geo_ret <- apply(1+S,2,function(x) prod(x)^(1/dim(S)[1])-1 ) #prod(1+S)^(1/dim(S)[1])-1
  
  Min_ret <- apply(S,2,min)
  
  MXDD_10 <- apply(S,2,function(x) max(rollapply(x,10,PerformanceAnalytics::maxDrawdown)))
  
  S_td <- PerformanceAnalytics::StdDev(S)
  shar_pe <- colMeans(S - Rf)/S_td
  sk_ew <- PerformanceAnalytics:: skewness(S)
  kurt <- PerformanceAnalytics::kurtosis(S)
  V_aR <- PerformanceAnalytics:: VaR(S)
  C_VaR <- PerformanceAnalytics::ETL(S)
  perf_mat <- c(cum_ret,Mean_Geo_ret,Min_ret,MXDD_10,S_td,shar_pe,sk_ew,kurt,V_aR,C_VaR)
  names(perf_mat) <- c("Cumulative_Return","Mean_Geometric_Return","Min_Return","10_days_MaxDrawdown","Volatility","Sharpe_Ratio","skewness","kurtosis","VaR","CVaR")
  t(perf_mat)
}

#apply(, 2, perf_set_security)


performance_security <- list()
performance_security$Bf_c <- apply(Bf_c[,2:14],2,function(x) perf_set_security(Bf_c[,"Row.names"],x,Bf_c$RF))
performance_security$Med_c <- apply(Med_c[,2:14],2,function(x) perf_set_security(Med_c[,"Row.names"],x,Med_c$RF))
performance_security$Aft_c <- apply(Aft_c[,2:14],2,function(x) perf_set_security(Aft_c[,"Row.names"],x,Aft_c$RF))

write.csv(performance_security, file = "performance_security.csv")

# sp_y <- All_Data_RetlessRF[bf_c1:af_c2,]
# perf_set_security(sp_y[,"Row.names"],sp_y$spy,sp_y$RF)

#write.csv(rho_ETFs_set$Aft_long, file = "test_data.csv")

#write.csv(regre_all_result$Aft_long, file = "test_data_beta.csv")


plotRet <- function(weights,return,na_me = "Portfolio")
  {
    a <- as.data.frame(t(matrix(weights,ncol = 13) %*% t(return[,2:14])))
    rownames(a) <- return[,"Row.names"]
    colnames(a) <- na_me
    x11()
    PerformanceAnalytics::charts.PerformanceSummary(a,methods = "ModifiedVaR",
                                                    wealth.index = TRUE)
}

plotRet(opt_beta_w$Bf_c,Bf_c)
plotRet(opt_beta_w$Med_c,Med_c)
plotRet(opt_beta_w$Aft_c,Aft_c)

