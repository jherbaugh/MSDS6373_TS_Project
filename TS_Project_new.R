

library(nnfor)
library(forecast)
library(vars)
library(tswge)
library(naniar)
library(DataExplorer)
library(GGally)

hk_eco <- read.csv(file.choose(), header = TRUE)

head(hk_eco)


##### Doing EDA on data

hk_eda <- hk_eco

str(hk_eda)

summary(hk_eda)

ggpairs(hk_eda[2:14])



#choosing hk_eda$HSI...volume to run Univarate Model
hk <- (hk_eco$HSI_volume)

#scale data by dividing by 1 million to normalize data
hk.scale <- hk_eco$HSI_volume/1000000

#plotting different timeframes in data to game insight
plotts.sample.wge(hk.scale)# True Acf is slowly dumping  with spectral density showing evidence of a wandering behavior at peak f=0 as seen in the realization
# with a hint of cyclic behavior gives an indication of a dominant (1-B) present hence we proceed to difference the data


# Taking First difference to remove (1-B) from full  scaled data
hk.df1= artrans.wge(hk.scale,phi.tr = 1)#Differencing the data shows that the sample autocorrelation looks white noise  and come cyclic behavior but the realization does not show white noise.
# we will proceed to model the data


# Model identification for the transformed or difference data using AIC and BIC approach
aic5.wge(hk.df1, p= 0:13, q= 0:5, type = "aic")

aic5.wge(hk.df1, p= 0:13, q= 0:5, type = "bic")

# Since AIC value is lower compared to that of BIC, we use ARMA(3,4) from the AIC results


# Use estimate p & q to get estimates of phi and thetas
est.hk <- est.arma.wge(hk.df1,p=3, q=4)

est.hk$avar

mean(hk.scale)

#Final Model
# (1-B)(1 + 0.462B - 0.060B^2 + 0.633B^3)(X_t + 1492.191) = (1 - 0.121B - 0.447B^2 - 0.705B^3 + 0.327B^4)a_t  $avar = 256933.7


########### forecast model forward 1000 volumes into future

#ARMA(3,4)
fore.arma.wge(hk.scale, phi = est.hk$phi,theta = est.hk$theta, n.ahead = 180, lastn = FALSE,limits = FALSE)


#ARIMA(3,1,4)
fore.aruma.wge(hk.scale,d=1,est.hk$phi,theta = est.hk$theta, n.ahead = 180, lastn = FALSE,limits = FALSE)



# Calculating ASE for Last 180 observations
n = 180
x = hk.scale
x.pred.scale = fore.arma.wge(hk.scale, phi = est.hk$phi,theta = est.hk$theta, n.ahead = 180, lastn = FALSE,limits = FALSE)
ASE.scale = mean((x[(length(x)-n+1):(length(x))]-x.pred.scale$f)^2)
print(paste0('ASE: ',round(ASE.scale,2)))

#ARMA ASE: 252237.46



# Calculating ASE for Last 180 observations of ARIMA model
n = 180
x = hk.scale
x.pred.scale1 = fore.aruma.wge(hk.scale,d=1,est.hk$phi,theta = est.hk$theta, n.ahead = 180, lastn = FALSE,limits = FALSE)
ASE.scale1 = mean((x[(length(x)-n+1):(length(x))]-x.pred.scale1$f)^2)
print(paste0('ASE: ',round(ASE.scale1,2)))

#ARIMA ASE: 233115.11


#Rolling Window ASE

Rolling_Window_ASE = function(series, trainingSize, horizon = 1, s = 0, d = 0, phis = 0, thetas = 0)
{
  ASEHolder = numeric()

  for( i in 1:(length(series)-(trainingSize + horizon) + 1))
  {

    forecasts = fore.aruma.wge(series[i:(i+(trainingSize-1))],phi = phis, theta = thetas, s = s, d = d,n.ahead = horizon)

    ASE = mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts$f)^2)

    ASEHolder[i] = ASE

  }

  ASEHolder
  hist(ASEHolder)
  WindowedASE = mean(ASEHolder)

  print("The Summary Statistics for the Rolling Window ASE Are:")
  print(summary(ASEHolder))
  print(paste("The Rolling Window ASE is: ",WindowedASE))
  return(WindowedASE)
}


Rolling_Window_ASE(hk.scale, s = 12, phis = est.hk$phi, horizon = 15, trainingSize = 100)


###################### using Private.Domestic..Price.Index ##########################################################

hk_ind <- hk_eda$Private_Domestic.Price_Index.


#plotting different time frames in data to game insight
plotts.sample.wge(hk_ind)# realization is shows wandering behavior with no sign of cyclic behavior
# it also has its True sample autocorrelation showing identical rho_k = 1
#spectral density shows a dominant (1-B) behavior at peak f=0 hence an indication for an ARIMA model we then proceed to take the first different

# Taking First difference to remove (1-B) from full  scaled data
hk.ind.df1= artrans.wge(hk_ind,phi.tr = 1) # Taking the first difference of the data (1-B) to stationarize the data

acf(hk.ind.df1) # Acfs are slowing dumping with a slight hint of cyclic behavior

hk.ind.df2= artrans.wge(hk.ind.df1,phi.tr = 1) #taking the second difference of the first difference (1-B)^2 shows significance at lag 1 with hint of cyclic behavior

hk.ind.df3 = artrans.wge(hk.ind.df2, phi.tr = c(rep(0,11),1))# taking the 12th difference of the data (1-B^12) to stationarize the data but shows seasonality at lag 12
#hence we remove the seasonality in the data with this differencing like in the airline model

plotts.sample.wge(hk.ind.df3) # sample Autocorrelations shows seasonal behavior at lag 12
#spectral density of the realization has no peak at f=o but has 8 peaks between f=0 and f=0.5

parzen.wge(hk.ind.df3) # wander with a strong indication of cyclic behavior

# Find a model for the transformed or difference data using AIC and BIC with the assumption data is stationary
aic5.wge(hk.ind.df3, p= 0:13, q= 0:5, type = "aic")

aic5.wge(hk.ind.df3, p= 0:13, q= 0:5, type = "bic")


# Use estimate p & q to get estimates of phi and thetas
est.hk.ind <- est.arma.wge(hk.ind.df3,p=12, q=5)

est.hk.ind$avar

mean(hk_ind)



#Final Model
# (1-B)^2 (1-B^12)(1 -0.0065B - 0.0053B^2 + 0.0581B^3 + 0.2819B^4 + 0.2321B^5 + 0.2023B^6 + 0.1594B^7 + 0.1309B^8 + 0.1402B^9 + 0.0337^10 - 0.0020B^11 + 0.4638B^12)(X_t + 197.065) = (1 - 1.7602B + 0.7959B^2 + 0.0444B^3 - 0.2856B^4)a_t  $avar = 0.0955


# forecast model forward 180 observations into future

#ARMA(12,5)
fore.arma.wge(hk_ind, phi = est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 180, lastn = FALSE,limits = FALSE)


#ARIMA(12,2,5)
fore.aruma.wge(hk_ind,d=1,s =12,est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 180, lastn = FALSE,limits = FALSE)



# Calculating ASE for Last 180 observations
n = 180
x = hk_ind
x.pred.ind = fore.arma.wge(hk_ind, phi = est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 180, lastn = FALSE,limits = FALSE)
ASE.ind = mean((x[(length(x)-n+1):(length(x))]-x.pred.ind$f)^2)
print(paste0('ASE: ',round(ASE.ind,2)))

#ARMA ASE: 313876.36



# Calculating ASE for Last 180 observations using ARIMA Model
n = 180
x = hk_ind
x.pred.ind1 = fore.aruma.wge(hk_ind,d=1,phi= est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 180, lastn = TRUE,limits = FALSE)
ASE.ind1 = mean((x[(length(x)-n+1):(length(x))]-x.pred.ind1$f)^2)
print(paste0('ASE: ',round(ASE.ind1,2)))

#ARIMA ASE: 1738977.48


Rolling_Window_ASE(hk_ind, s = 12, phis = est.hk.ind$phi,theta = est.hk.ind$theta, horizon = 15, trainingSize = 1000)




# end of code

#########################################################################################################################
























#choosing hk_eda$First_hand_sales_quantity to run Univarate Model

hk_close <- (hk_eda$HSI...close)


#plotting different timeframes in data to game insight
plotts.sample.wge(hk_close)# realization shows wandering behavior with an indication of cyclic behavior
# it also has its True sample autocorrelation showing identical rho_k = 1
#spectral density shows a dominant (1-B) behavior at peak f=0 hence an indication for an ARIMA model suited for this data hence we then proceed to take the first different
# an indication of nonstationary process

#plotts.sample.wge(hk_close[1:925]) #perhaps before 2008 economic or financial crisis

#plotts.sample.wge(hk_close[926:4233]) #perhaps after 2008 economic or financial crisis



# Taking First difference to remove (1-B) from full  scaled data
hk_close.df1= artrans.wge(hk_close,phi.tr = 1)# first differencing shows sample autocorrelation depicting white noise
#Although the dominant (1-B) is removed, it still has a seasonal treand to it hence will require removing the (1-B)^12

# Taking First difference to remove (1-B^12) from full  scaled data
hk_close.df12= artrans.wge(hk_close.df1, phi.tr = c(rep(0,11),1))


# Find a model for the transformed or difference data using AIC and BIC using the Second difference
aic5.wge(hk_close.df12, p= 0:13, q= 0:5, type = "aic")

aic5.wge(hk_close.df12, p= 0:13, q= 0:5, type = "bic")



# Since AIC value is lower compared to that of BIC, we use ARMA(3,4) from the AIC results


# Use estimate p & q to get estimates of phi for  original data
est.hk.close <- est.ar.wge(hk_close.df12,p=12, type = "burg")

est.hk.close$avar

mean(hk_close)

#Final Model
# (1-0.9989B)(X_t + 20781.49) = a_t  $avar = 80411.43


# Use estimate p & q to get estimates of phi for DF12
#est.hk.close12 <- est.ar.wge(hk_close.df12,p=1)

#est.hk.close12$avar

#mean(hk_close)

#Final Model for df12
# (1-0.914857B)(X_t + 20781.49) = a_t  $avar = 147220.8


########### forecast model forward 180 volumes into future

# Forecast with AR(1) original data

fore.arma.wge(hk_close, phi = est.hk.close$phi, n.ahead = 180, lastn = FALSE,limits = FALSE)

#AR(1) forecast for original data d=1
#fore.aruma.wge(hk_close,est.hk.close$phi, n.ahead = 1000, lastn = FALSE,limits = FALSE)


# Forecast with AR(1)  forDF12
#fore.arma.wge(hk_close, phi = est.hk.close12$phi, n.ahead = 180, lastn = FALSE,limits = FALSE)


#AR(1) forecast for df12 s=12
fore.aruma.wge(hk_close,s=12,est.hk.close12$phi, n.ahead = 180, lastn = FALSE,limits = FALSE)


# Calculating ASE for Last 1000 observations
n = 180
x = hk_close
x.pred.close = fore.arma.wge(hk_close, phi = est.hk.close$phi, n.ahead = 1000, lastn = FALSE,limits = FALSE)
ASE.close = mean((x[(length(x)-n+1):(length(x))]-x.pred.close$f)^2)
print(paste0('ASE: ',round(ASE.close,2)))

# ASE: 12812559.32" # AR(1)



# Calculating ASE for Last 1000 observations
n = 180
x = hk_close
x.pred.close12 = fore.arma.wge(hk_close, phi = est.hk.close12$phi, n.ahead = 180, lastn = FALSE,limits = FALSE)
ASE.close12 = mean((x[(length(x)-n+1):(length(x))]-x.pred.close12$f)^2)
print(paste0('ASE: ',round(ASE.close12,2)))


#ASE: 43388280.54"


# Calculating ASE for Last 1000 observations for DF 12
n = 1000
x = hk_close
x.pred.scale = fore.arma.wge(hk_close, phi = est.hk.salesquant$phi, n.ahead = 1000, lastn = FALSE,limits = FALSE)
ASE.scale = mean((x[(length(x)-n+1):(length(x))]-x.pred.scale$f)^2)
print(paste0('ASE: ',round(ASE.scale,2)))




# Calculating ASE for Last 1000 observations for DF12 s=12
n = 1000
x = hk_close
x.pred.close12 = fore.aruma.wge(hk_close,s=12,est.hk.close12$phi, n.ahead = 180, lastn = FALSE,limits = FALSE)
ASE.close12 = mean((x[(length(x)-n+1):(length(x))]-x.pred.close12$f)^2)
print(paste0('ASE: ',round(ASE.close12,2)))

















################ perhaps before 2008 economic or financial crisis
hk.pre.2008 <- hk.scale[1:925]

plotts.sample.wge(hk.pre.2008)


# Taking First difference to remove (1-B)
hk.pre2008.df1= artrans.wge(hk.pre.2008,phi.tr = 1)
 acf(hk.pre2008.df1)

hk.pre2008.df2 = artrans.wge(hk.pre2008.df1, phi.tr = c(rep(0,11),1))

# Find a model for the transformed or difference data using AIC and BIC
aic5.wge(hk.pre2008.df1, p= 0:13, q= 0:5, type = "aic")

aic5.wge(hk.pre2008.df1, p= 0:13, q= 0:5, type = "bic")


# Find a model for the transformed or difference data using AIC and BIC
aic5.wge(hk.pre2008.df2, p= 0:13, q= 0:5, type = "aic")

aic5.wge(hk.pre2008.df2, p= 0:13, q= 0:5, type = "bic")



# Since BIC value is lower compared to that of AIC, we use ARMA(1,1) from the BIC results



# Use estimate p & q to get estimates of phi and thetas
est.hk.pre2008 <- est.arma.wge(hk.pre2008.df1,p=1, q=1)

est.hk.pre2008$avar

mean(hk.pre.2008)



#Final Model
# (1-B)(1 + 0.215B)(X_t + 316.9835) = (1 - 0.900B)a_t  $avar = 10849.39



# Use estimate p & q to get estimates of phi and thetas
est12.hk.pre2008 <- est.arma.wge(hk.pre2008.df2,p=12, q=1)

est12.hk.pre2008$avar

mean(hk.pre.2008)
########### forecast model forward 1000 volumes into future

#ARMA(1,1)
fore.arma.wge(hk.pre2008.df2, phi = est12.hk.pre2008$phi,theta = est12.hk.pre2008$theta, n.ahead = 1000, lastn = FALSE,limits = FALSE)



#ARIMA(1,1,1)
fore.aruma.wge(hk.pre2008.df2,d=1,est.hk.pre2008$phi,theta = est.hk.pre2008$theta, n.ahead = 1000, lastn = FALSE,limits = FALSE)




# Calculating ASE for Last 1000 observations
x.2008 = hk.pre.2008
n = 1000
x.pred08 = fore.arma.wge(hk.pre2008.df2, phi = est12.hk.pre2008$phi,theta = est12.hk.pre2008$theta, n.ahead = 1000, lastn = FALSE,limits = FALSE)
ASE = mean((x.2008[(length(x.2008)-n+1):(length(x.2008))]-x.pred08$f)^2)
print(paste0("ARMA_ASE_:",round(ASE,2)))


x1.2008 <- hk.pre.2008
n = 1000
x1.pred = fore.aruma.wge(x1.2008,d=1,est.hk.pre2008$phi,theta = est.hk.pre2008$theta, n.ahead = 1000, lastn = FALSE,limits = FALSE)
ASE.d1 = mean((x1.2008[(length(x1.2008)-n+1):(length(x1.2008))]-x1.pred$f)^2)
print(paste0("ARIMA_ASE_:",round(ASE.d1,2)))






plotts.sample.wge(hk_eda$Private.Domestic..Price.Index) #Price index of private housing in HK - Y


