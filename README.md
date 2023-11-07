---
title: "Financial Mathematics Case Study"
author: "Massan Sarsenbayev"
output:
  bookdown::pdf_document2
documentclass: report
geometry: "left=3cm,right=3cm,top=2cm,bottom=2cm"
classoption: a4paper
always_allow_html: true
---



## Introduction
This report aims to analyze a portfolio of seven selected stocks, including their historical prices, expected returns, volatility, correlations, and risk. The portfolio will be optimized for a target return using an appropriate solver function, and the efficient frontier curve will be plotted for different target returns. The Sharpe ratios will be calculated using a risk-free investment, and the equation of the Capital Market Line will be determined. Linear regression analysis will be used to calculate the beta of each asset in the portfolio and to estimate the Value at Risk (VaR) of the portfolio. The volatility of a single asset in the portfolio will be estimated using ARCH/GARCH models, and the best model will be identified. Finally, the findings will be summarized in a non-technical language for a potential investor, and other relevant performance measurements will be discussed.


```{r, include=FALSE}
library(xts)
library(zoo)
library(readr)
library(foreach)
library(DEoptim)
library(iterators)
library(fGarch)
library(Rglpk)
library(ROI)
library(ROI.plugin.glpk)
library(ROI.plugin.quadprog)
library(plotly)
library(dplyr)
library(quantmod)
library(PerformanceAnalytics)
library(imputeTS)
library(corrplot)
library(ggplot2)
library(ggcorrplot)  
library(PortfolioAnalytics)
library(tseries)
library(rugarch)
```


###########################################
II
 *Use an appropriate Solver function to determine the portfolio risk and the percentage of investment in each asset in your portfolio for a target return of your choice.
 *draw the efficient 	frontier curve

https://www.codingfinance.com/post/2018-05-31-portfolio-opt-in-r/


```{r}


stock_tickers <- c("FANG", "ATCO", "SHEL", "AAPL", "NGLOY", "AMRK", "TJX");

portfolioPrices <- NULL
for(ticker in stock_tickers) {
  portfolioPrices <- cbind(portfolioPrices,
  getSymbols.yahoo(ticker, from='2015-12-14', periodicity = 'daily', 
                   auto.assign=FALSE)[,4])
};


portfolioReturns <- na.omit(diff(log(portfolioPrices)));
print("Expected returns of stocks in portfolio:")
daily_returns <- colMeans(portfolioReturns);
daily_returns ;

annual_returns <- daily_returns*250

print("Daily volatilities of stocks:")
volatilities <- StdDev(portfolioReturns);
volatilities;

print("Correlations of stocks in portfolio:")
correlation <- cor(portfolioReturns);
correlation

print("Variance Covariance Matrix:")
Var_Cov_Matrix<-cov(portfolioReturns)*250;
Var_Cov_Matrix;


```

```{r}
simpleReturns <- NULL
for(i in 1:ncol(portfolioPrices)) {
  simpleReturns <- cbind(portfolioReturns, diff(portfolioPrices[,i]) / portfolioPrices[-1,i])
};
```


```{r}
# plotting the correlation heatmap
ggcorrplot(correlation);
```


```{r}
portf <- portfolio.spec(colnames(portfolioReturns))

portf <- add.constraint(portf, type="weight_sum", min_sum=1, max_sum=1)
portf <- add.constraint(portf, type="box", min=0.05, max=0.40)
portf <- add.objective(portf, type="return", name="mean")
portf <- add.objective(portf, type="risk", name="StdDev")
portf <- add.constraint(portf, type = "mean", mean_target = 0.20)


optPort <- optimize.portfolio(portfolioReturns, portf, 
                              optimize_method = "ROI", trace=TRUE)


print("Optimum weights of stocks in portfolio:")
weights(optPort)

weights <- weights(optPort)

portf_risk <- sqrt(t(weights) %*% (Var_Cov_Matrix %*% weights))
print(paste0("Portfolio risk: ", portf_risk))

chart.Weights(optPort)
```

```{r}

w1 <- c(0.0500000,   0.0500000,   0.0500000,   0.3480800,   0.2851958,   0.1667242, 0.0500000) 
names (w1) <- stock_tickers
print(w1)

## To construct the portfolio weights (equally weighted) 
v <- c(1,1,1,1,1,1,1)
names (v) <- stock_tickers
w2 <- c(v/sum(v))
print(w2)

w3 <- c(0.09508637, -0.30551382, -0.83734765,  0.82794644,  0.72427894,  0.30812314,  0.18742657)
names (w3) <- stock_tickers
print(w3)

w4 <- c(0.1000000,   0.1000000,   0.1000000,   0.2429739,   0.2187530,   0.1382731,   0.1000000 )
names (w4) <- stock_tickers
print(w4)

w5 <- c(1.275774e-17,  3.285755e-17, -6.658931e-17,  4.531861e-01,  3.516386e-01,  1.951754e-01, -2.775558e-17)
names (w5) <- stock_tickers
print(w5)

```
```{r}
w1_expected_returns <- sum(annual_returns*w1)
w1_expected_returns
w2_expected_returns <- sum(annual_returns*w2)
w2_expected_returns
w3_expected_returns <- sum(annual_returns*w3)
w3_expected_returns
w4_expected_returns <- sum(annual_returns*w4)
w4_expected_returns
w5_expected_returns <- sum(annual_returns*w5)
w5_expected_returns
expected_returns <- c(0.2189405, 0.1464668, 0.4770716, 0.1811587, 0.2567224)
expected_returns
```
```{r}
expected_returns123 <- c(w1_expected_returns, w2_expected_returns, w3_expected_returns, w4_expected_returns, w5_expected_returns)
expected_returns123
```


```{r}
#w1 StdDev
utc <- upper.tri(Var_Cov_Matrix)
wt.var <- sum(diag(Var_Cov_Matrix) * w1^2)
wt.cov <- sum(w1[row(Var_Cov_Matrix)[utc]] *
                w1[col(Var_Cov_Matrix)[utc]] *
                Var_Cov_Matrix[utc])
variance <- wt.var + 2 * wt.cov
w1_standard_dev <- sqrt(variance)
print(w1_standard_dev)
#w2 StdDev
utc <- upper.tri(Var_Cov_Matrix)
wt.var <- sum(diag(Var_Cov_Matrix) * w2^2)
wt.cov <- sum(w2[row(Var_Cov_Matrix)[utc]] *
                w2[col(Var_Cov_Matrix)[utc]] *
                Var_Cov_Matrix[utc])
variance <- wt.var + 2 * wt.cov
w2_standard_dev <- sqrt(variance)
print(w2_standard_dev)
#w3 StdDev
utc <- upper.tri(Var_Cov_Matrix)
wt.var <- sum(diag(Var_Cov_Matrix) * w3^2)
wt.cov <- sum(w3[row(Var_Cov_Matrix)[utc]] *
                w3[col(Var_Cov_Matrix)[utc]] *
                Var_Cov_Matrix[utc])
variance <- wt.var + 2 * wt.cov
w3_standard_dev <- sqrt(variance)
print(w3_standard_dev)
#w4 StdDev
utc <- upper.tri(Var_Cov_Matrix)
wt.var <- sum(diag(Var_Cov_Matrix) * w4^2)
wt.cov <- sum(w4[row(Var_Cov_Matrix)[utc]] *
                w4[col(Var_Cov_Matrix)[utc]] *
                Var_Cov_Matrix[utc])
variance <- wt.var + 2 * wt.cov
w4_standard_dev <- sqrt(variance)
print(w4_standard_dev)
#w5 StdDev
utc <- upper.tri(Var_Cov_Matrix)
wt.var <- sum(diag(Var_Cov_Matrix) * w5^2)
wt.cov <- sum(w5[row(Var_Cov_Matrix)[utc]] *
                w5[col(Var_Cov_Matrix)[utc]] *
                Var_Cov_Matrix[utc])
variance <- wt.var + 2 * wt.cov
w5_standard_dev <- sqrt(variance)
print(w5_standard_dev)

volatilities <- c(0.2707854, 0.2654307, 0.4739584, 0.2648562, 0.2878772)
volatilities
```


The optimal portfolio risk is the level of risk that is most appropriate for an investor's risk tolerance and investment objectives. It is the level of risk that maximizes the expected return of the portfolio for a given level of risk, or minimizes the risk for a given level of expected return.

The optimal portfolio risk can be determined using the efficient frontier curve, which is a graphical representation of the trade-off between expected return and risk for a given set of investments. The efficient frontier curve is used to identify the optimal portfolios that provide the highest expected return for a given level of risk, or the minimum risk for a given level of expected return.

```{r}

ef <- extractEfficientFrontier(optPort, match.col = "StdDev", n.portfolios = 20,
                         risk_aversion = NULL)

chart.EfficientFrontier(ef,
            match.col = "StdDev", n.portfolios = 20, xlim = NULL, 
            ylim = NULL, cex.axis = 0.8, element.color = "darkgray", 
            main = "Efficient Frontier",
            RAR.text = "SR", rf = 0, tangent.line = TRUE, cex.legend = 0.8,
            chart.assets = TRUE, labels.assets = TRUE, pch.assets = 21,
            cex.assets = 0.8)
```


```{r}
# Generate vector of returns
#minret <- 0.02/100
expected_returns <- expected_returns

minret <- min(expected_returns)
maxret <- max(expected_returns)
#maxret <- maxret.opt$weights %*% meanReturns
 
vec <- seq(minret, maxret, length.out = 100)



eff.frontier <- data.frame(Risk =vector("numeric", length(vec)) ,
                           Return = vector("numeric", length(vec)))
 
wts <- mat.or.vec(nr = length(vec), nc = ncol(portfolioReturns))
colnames(wts) <- colnames(portfolioReturns)
 
for(i in 1:length(vec)){
  
  # Creates a new portfolio object using p and adds mean as an objective
  
  portf <- add.constraint(portf, type = "return", name = "mean", return_target = vec[i])
  
  # Creates a new portfolio object using p and adds var as an objective
  portf <- add.objective(portf, type = "risk", name = "var")
  
# Creates a new portfolio object using p and adds a weight_concentration
# objective. The conc_aversion parameter controls how much concentration is
# penalized. The portfolio concentration is defined as the Herfindahl Hirschman
# Index of the weights.
  portf <- add.objective(portf, type = "weight_concentration", name = "HHI",
                             conc_aversion = 0.01)
 
  eff.opt <- optimize.portfolio(portfolioReturns, portf, optimize_method = "ROI")
  
  eff.frontier$Risk[i] <- sqrt(t(eff.opt$weights) %*% Var_Cov_Matrix %*% eff.opt$weights)
  
  eff.frontier$Return[i] <- eff.opt$weights %*% expected_returns
  
  
  
  wts[i,] = eff.opt$weights
  
 # print(paste(round(i/length(vec) * 100, 0), "% done..."))
  
}
eff.frontier$Sharperatio <- eff.frontier$Return / eff.frontier$Risk

```


```{r}


# Generate random portfolios
randomport<- random_portfolios(portf, permutations = 50, rp_method = "sample")

portf <- add.constraint(portfolio = portf, type = "full_investment")
portf <- add.constraint(portf, type="long_only")
# Get minimum variance portfolio
minvar.port <- add.objective(portf, type = "risk", name = "var")
 
# Optimize
minvar.opt <- optimize.portfolio(portfolioReturns, minvar.port, optimize_method = "random", 
                                 rp = randomport)
 
# Generate maximum return portfolio
maxret.port <- add.objective(portf, type = "return", name = "mean")
 
# Optimize
maxret.opt <- optimize.portfolio(portfolioReturns, maxret.port, optimize_method = "random", 
                                 rp = randomport)



```


```{r}

feasible.sd <- apply(randomport, 1, function(x){
  return(sqrt(matrix(x, nrow = 1) %*% Var_Cov_Matrix %*% matrix(x, ncol = 1)))
})
 
feasible.means <- apply(randomport, 1, function(x){
  return(x %*% expected_returns)
})
 
feasible.sr <- feasible.means / feasible.sd
portf <- plot_ly() %>%
  add_trace(x = feasible.sd, y = feasible.means, color = feasible.sr, 
        mode = "markers", type = "scattergl", showlegend = F,
        
        marker = list(size = 3, opacity = 0.9, 
                      colorbar = list(title = "Sharpe Ratio"))) %>%
  add_trace(data = eff.frontier, x = ~Risk, y = ~Return, mode = "markers", type = "scattergl")%>% 
  layout(title = "Efficient Frontier",
         yaxis = list(title = "Mean Returns", tickformat = ".2%"),
         xaxis = list(title = "Standard Deviation", tickformat = ".2%"))

portf


#http://rstudio-pubs-static.s3.amazonaws.com/333559_274650a583524ac1b086f79a7133fc28.html

```



###########################################
III (may be not correct)
  *Calculate the Sharpe ratios 
  *Determine the equation of the Capital Market Line by using a risk free investment       guaranteeing a return of 1.5%

The beta of a portfolio is a measure of the portfolio's volatility relative to the market. The market is usually represented by a benchmark such as the S&P 500. A beta of 1 means that the portfolio has the same level of volatility as the market, while a beta greater than 1 indicates higher volatility and a beta less than 1 indicates lower volatility.

To calculate the beta of a portfolio using the CAPM module, you will need to provide the following arguments:
* portfolioReturn: This is the return of the portfolio over a period of time. This can be a series of daily, weekly, or monthly returns.
* benchmarkReturns: This is the return of the benchmark over the same period of time as the portfolio returns.
* riskFreeRate: This is the risk-free rate of return. This is usually a very low rate such as the return on a short-term government bond.

The CAPM.beta.bull and CAPM.beta.bear functions calculate the bull and bear betas of a portfolio, which are measures of the portfolio's sensitivity to market movements during periods of rising and falling markets, respectively.

The CAPM (Capital Asset Pricing Model) module provides two functions for calculating the alpha of a portfolio: CAPM.alpha and CAPM.jensenAlpha.

The alpha of a portfolio is a measure of the excess return of the portfolio relative to the expected return based on the market and the portfolio's beta. A positive alpha indicates that the portfolio has outperformed the market, while a negative alpha indicates underperformance.

To calculate the alpha of a portfolio using the CAPM module, you will need to provide the following arguments:

portfolioReturn: This is the return of the portfolio over a period of time. This can be a series of daily, weekly, or monthly returns.

benchmarkReturns: This is the return of the benchmark over the same period of time as the portfolio returns.

riskFreeRate: This is the risk-free rate of return. This is usually a very low rate such as the return on a short-term government bond.

The CAPM.jensenAlpha function calculates the Jensen alpha of a portfolio, which is a variant of the alpha that adjusts for the portfolio's beta.

```{r}

tickers <- stock_tickers
weights <- weights

benchmarkPrices <- getSymbols.yahoo("SPY", from="2015-12-14", 
                          periodicity = "daily", auto.assign=FALSE)[,4]
colSums(is.na(benchmarkPrices))
benchmarkReturns <- na.omit(diff(log(benchmarkPrices)))


#Rename Columns
colnames(portfolioPrices) <- tickers

#Get sum of NA per column
colSums(is.na(portfolioPrices))

#Plot
plot(portfolioPrices, legend = tickers)


#Calculate Returns For DF
dailyReturns <- na.omit(ROC(portfolioPrices, type="discrete"))

#Calculate Portfolio Returns
portfolioReturn <- Return.portfolio(portfolioReturns, weights=w5)

#Plot Performance
chart.CumReturns(portfolioReturn)
charts.PerformanceSummary(portfolioReturn)

#Calculate Metrics 
print("beta")
CAPM.beta(portfolioReturn, benchmarkReturns, .015/250)
print("beta bull")
CAPM.beta.bull(portfolioReturn, benchmarkReturns, .015/250)
print("beta bear")
CAPM.beta.bear(portfolioReturn, benchmarkReturns, .015/250)

print("jensenAlph")
#CAPM.alpha(portfolioReturn, benchmarkReturns, .035/252)
CAPM.jensenAlpha(portfolioReturn, benchmarkReturns, .015/250)

SharpeRatio(portfolioReturn, Rf = .015/250, p = 0.95, FUN = "StdDev",
            weights = NULL, annualize = FALSE)

table.AnnualizedReturns(portfolioReturn, Rf=.015/250, geometric=TRUE)



```


###########################################
III
  *Calculate the Sharpe ratios for a range of expected portfolio returns and volatilities 
  *Determine the equation of the Capital Market Line by using a risk free investment guaranteeing a return of 1.5%
  
  
You can also use the SharpeRatio function from the PerformanceAnalytics package to calculate the Sharpe ratio of a portfolio. This function allows you to specify the returns of the portfolio, the risk-free rate, and the confidence level to use for the volatility calculation.

```{r}
# Set the risk-free rate
Rf <- 0.0

# Calculate the Sharpe ratio of the portfolio
sharpe_ratio <- SharpeRatio(portfolioReturn, Rf = Rf, p = 0.95, FUN = "StdDev",
                            weights = NULL, annualize = TRUE)

# Print the Sharpe ratio
print(paste("Sharpe ratio of the portfolio:", round(sharpe_ratio, 3)))
```

The Sharpe ratio is a measure of the risk-adjusted return of an asset or a portfolio, calculated as the excess return over the risk-free rate divided by the volatility of the returns. It is used to compare the performance of different investments and to evaluate the trade-off between risk and return.

To calculate the Sharpe ratio for a range of expected portfolio returns and volatilities, you can use a loop to iterate over the different values of returns and volatilities, and calculate the Sharpe ratio for each combination.

Here is an example of how you can calculate the Sharpe ratios for a range of expected portfolio returns and volatilities in R:

In this example, the expected_returns and volatilities sequences define the range of expected returns and volatilities, respectively. The sharpe_ratios matrix is initialized with NA values and will be used to store the Sharpe ratios calculated in the loop. The loop iterates over the expected_returns and volatilities sequences, and calculates the Sharpe ratio for each combination using the formula (expected_return - Rf)/volatility. The resulting Sharpe ratios are stored in the sharpe_ratios matrix, which can then be printed or plotted to visualize the results.

```{r}

# Set the risk-free rate
Rf <- 0.0

# Define the range of expected returns and volatilities
expected_returns <- as.vector(expected_returns)
volatilities <- as.vector(volatilities)

stock_tickers_list <- list(stock_tickers, 
                      stock_tickers)

# Create a matrix to store the Sharpe ratios
sharpe_ratios_Matrix <- matrix(NA, nrow = length(expected_returns), 
                     ncol = length(volatilities), dimnames = stock_tickers_list)

# Loop over the expected returns and volatilities
for (i in 1:length(expected_returns)) {
  for (j in 1:length(volatilities)) {
    # Calculate the Sharpe ratio
    sharpe_ratios_Matrix[i, j] <- (expected_returns[i] - Rf)/volatilities[j]
  }
}

# Print the Sharpe ratios
print(sharpe_ratios_Matrix)
```

Stock Sharpe ratios can be negative for a variety of reasons. One possible reason is that the stock has performed poorly and has not provided a sufficient return to compensate for the risk taken on by the investor. Another possibility is that the stock market as a whole has experienced negative returns, which would also result in negative Sharpe ratios for individual stocks.

It's also worth noting that the Sharpe ratio is a relative measure that compares the risk and return of an investment to a benchmark, such as a risk-free asset or a market index. If the benchmark has a higher return than the investment, the Sharpe ratio will be negative, even if the investment has a positive return on its own.

```{r}
# Set the risk-free rate
Rf <- 0.0
portfolios <- c("w1", "w2", "w3", "w4", "w5")
volatilities <- as.vector(volatilities)
sharpe_ratios <- (expected_returns - Rf)/volatilities

sharpe_ratios_df <- data.frame(sharpe_ratios, colnames = portfolios)

colnames(sharpe_ratios_df) <- c("Sharpe Ratio", "Portfolios")

# Print the Sharpe ratios
print(sharpe_ratios_df)
```



```{r}
# Set the risk-free rate
Rf <- 0.0

# Set the confidence level
p <- 0.95

# Set the annualization factor
annualize <- 250

# Convert the prices to returns
asset_returns <- portfolioReturns

# Calculate the volatilities of the assets in the portfolio
volatilities <- as.vector(volatilities)

# Calculate the Sharpe ratios for each expected return
sharpe_ratios <- sapply(asset_returns, function(return, volatility) {
  SharpeRatio(R = return, Rf = Rf, p = p, annualize = annualize, 
              sigma = volatility)
}, volatility = volatilities)

# Plot the Sharpe ratios
plot(asset_returns, sharpe_ratios, xlab = "Expected Return", 
     ylab = "Sharpe Ratio",
     main = "Sharpe Ratios for Different Expected Returns")



```


The equation of the CML is of the form y = mx + b, where y is the Sharpe ratio, m is the slope of the CML (risk premium), x is the expected return, and b is the y-intercept (risk-free rate). We can use the slope and y-intercept of the CML to determine its equation.


```{r}

annual_expected_returns <- expected_returns*250

sharpe_ratios_list <- (expected_returns - Rf)/volatilities

# Plot the Sharpe ratios on a graph with risk on the x-axis and return on 
#the y-axis
plot(volatilities, annual_expected_returns, xlim = c(0, 0.6), ylim = c(0, 0.4),
 xlab = "Risk (volatility)", ylab = "Return", main = "Capital Market Line")

# Add the risk-free rate to the plot
abline(h = Rf, col = "red", lty = 2)

# Add the CML to the plot
abline(a = 0, b = max(sharpe_ratios_list), col = "blue", lty = 2)



# Add the efficient frontier to the plot
abline(eff.frontier, col = "green")

# Add a legend to the plot
legend("topright", c("Risk-free rate", "CML", "Efficient Frontier"), 
       lty = c(2, 2, 2), col = c("red", "blue", "green"))

```
```{r}
#Calculate Portfolio Returns
portfolioReturn <- Return.portfolio(portfolioReturns, weights=w5)
table.AnnualizedReturns(portfolioReturn, Rf=.015/250, geometric=TRUE)
```


###########################################
IV
Beta of assets in portfolio

https://www.investopedia.com/ask/answers/070615/what-formula-calculating-beta.asp



```{r}

# Join the data
PreturnBenchmJoined <- merge(portfolioReturn, benchmarkReturns)
# Join the data
DailyReturnsBenchmJoined <- merge(dailyReturns, benchmarkReturns)
# Save the xts object as a CSV file
write.csv(PreturnBenchmJoined, "PreturnBenchmJoined.csv")
# To read in the dataset in R
PreturnBenchmJoinedDate<-read.csv("PreturnBenchmJoined.csv", header=TRUE)
DailyReturnsBenchmJoinedDF <- as.data.frame(DailyReturnsBenchmJoined)

covarince <- cov(DailyReturnsBenchmJoinedDF)
covarince

stock_varinces <- covarince[,c("SPY.Close")]

varianceSPY <- var(DailyReturnsBenchmJoinedDF$SPY.Close)
varianceSPY

stock_varinces <- as.vector(stock_varinces)

stock_betas <- stock_varinces/varianceSPY
stock_betas

```


```{r}
library(tidyr)

beta.tab <- DailyReturnsBenchmJoinedDF %>% 
  gather(key = "asset", value = "return", -SPY.Close) %>%
  group_by(asset) %>%
  do(ols.model = lm(data = ., formula = return ~ SPY.Close)) %>% #estimate model
  mutate(beta = coef(ols.model)[2]) # get coefficients
print(beta.tab)


#https://www.r-bloggers.com/2017/01/how-to-calculate-betas-systematic-risk-for-a-large-number-of-stocks/
  
```
```{r}
mean(beta.tab$beta)
```


VAR(5%)

```{r}
PorVaR.Gaus<-VaR(portfolioReturns, p=0.95, 
           weights=w1,portfolio_method = "component",method="gaussian")
PorVaR.Gaus
```


###########################################
V
  *Estimate the volatility of a single asset in your portfolio using ARCH/GARCH and its extensions using R. 
  *Identify the best model and provide an explanation of your chosen model. 

```{r}
# create a zoo object for plotting
NGLOYPrices<- zoo(portfolioPrices$NGLOY)
RNGLOY<-ts(diff(log(NGLOYPrices)))

plot(RNGLOY, col='BLUE')

M10<-garch(RNGLOY, order=c(1,0), trace=FALSE)
M11<-garch(RNGLOY, order=c(1,1), trace=TRUE)
M20<-garch(RNGLOY, order=c(2,0), trace=FALSE)
M21<-garch(RNGLOY, order=c(2,1), trace=FALSE)
M22<-garch(RNGLOY, order=c(2,2), trace=FALSE)
M02<-garch(RNGLOY, order=c(0,2), trace=FALSE)
M12<-garch(RNGLOY, order=c(1,2), trace=FALSE)

AIC (M10, M11, M20, M21, M22, M02, M12)

#plot(M11)

summary(M11)
```

```{r}
summary(M11)
```



```{r}
library(fGarch)

# Fit a GARCH model to the asset returns
fit <- garchFit(formula = ~garch(1, 0), data = portfolioReturns$TJX.Close, 
                trace = TRUE)

# Extract the estimated volatility from the fitted model
volatilitygarch <- fitted(fit)
plot(volatilitygarch)
```


After fitting the models, the AIC function is used to compute the Akaike Information Criterion (AIC) for each model. The AIC is a measure of the goodness-of-fit of a statistical model, with a lower AIC indicating a better fit. In this case the lowest AIC value has a model g12, that indicates this is the best model.







```{r}
utc <- upper.tri(Var_Cov_Matrix)
wt.var <- sum(diag(Var_Cov_Matrix) * weights^2)
wt.cov <- sum(weights[row(Var_Cov_Matrix)[utc]] *
                weights[col(Var_Cov_Matrix)[utc]] *
                Var_Cov_Matrix[utc])
variance <- wt.var + 2 * wt.cov
standard_dev <- sqrt(variance)
print(standard_dev)
```


```{r}
annual_volatility <- sqrt(250)*volatilities
annual_volatility
```




```{r}
annual_expected_returns <- expected_returns*250
annual_expected_returns
```

```{r}
# create a zoo object for plotting
AMRKPrices<- zoo(portfolioPrices$AAPL)
RAMRK<-ts(diff(log(AMRKPrices)))
```


```{r}
g1<-ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1)),
               mean.model=list(armaOrder=c(1,0)),distribution.model="std")

```

```{r}
garch11<-ugarchfit(g1,data = NGLOYPrices)
garch11
```

```{r}
vole <- ts(garch11@fit$sigma^2, start = c(2015,12), end = c(2022,12), frequency = 72)
plot(vole, xlab="", ylab="",main="Volatility NGLOY")
```



