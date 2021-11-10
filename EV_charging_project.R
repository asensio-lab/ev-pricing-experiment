##### SECTION 0: LOAD REQUIRED PACKAGES, FUNCTIONS, AND DATA

## LOAD DATA

EV_data <- read.csv(file.choose()) #open sessions.csv

## LOAD PACKAGES AND DEFINE FUNCTIONS

# General display items.
library(ggplot2)
library(extrafont)
font_import()
loadfonts(device = 'pdf')
library(reshape)

# McCrary tests and display items.
library(rddensity)
library(lpdensity)
rdplotdensity2<-function (rdd, X, plot_title, plotRange = NULL, plotN = 10, 
                          plotGrid = c("es", "qs"), alpha = 0.01, type = NULL, CItype = NULL, 
                          title = "", xlabel = "", ylabel = "", lty = NULL, lwd = NULL, lcol = NULL, 
                          pty = NULL, pwd = NULL, pcol = NULL, CIshade = NULL, CIcol = NULL, 
                          legendTitle = NULL, legendGroups = NULL) 
{
  c <- rdd$opt$c
  p <- rdd$opt$p
  q <- rdd$opt$q
  hl <- rdd$h$left
  hr <- rdd$h$right
  kernel <- rdd$opt$kernel
  if (length(plotRange) == 0) {
    plotRange <- c(max(min(X), c - 3 * hl), min(max(X), c +  3 * hr))
  }
  else if (length(plotRange) != 2) {
    stop("Plot range incorrectly specified.\\n")
  }
  else if (plotRange[1] >= c | plotRange[2] <= c) {
    stop("Plot range incorrectly specified.\\n")
  }
  if (length(plotN) == 0) {
    plotN <- c(10, 10)
  }
  else if (length(plotN) == 1) {
    plotN <- c(plotN, plotN)
  }
  else if (length(plotN) > 2) {
    stop("Number of grid points incorrectly specified.\\n")
  }
  if (plotN[1] <= 1 | plotN[2] <= 1) {
    stop("Number of grid points incorrectly specified.\\n")
  }
  if (length(plotGrid) == 0) {
    plotGrid <- "es"
  }
  else {
    plotGrid <- plotGrid[1]
  }
  if (!plotGrid %in% c("es", "qs")) {
    stop("Grid specification invalid.\\n")
  }
  scalel <- (sum(X <= c) - 1)/(length(X) - 1)
  scaler <- (sum(X >= c) - 1)/(length(X) - 1)
  if (plotGrid == "es") {
    gridl <- seq(plotRange[1], c, length.out = plotN[1])
    gridl[plotN[1]] <- c
    gridr <- seq(c, plotRange[2], length.out = plotN[2])
    gridr[1] <- c
  }
  else {
    gridl <- seq(mean(X <= plotRange[1]), mean(X <= c), length.out = plotN[1])
    gridl <- quantile(X, gridl)
    gridr <- seq(mean(X <= c), mean(X <= plotRange[2]), length.out = plotN[2])
    gridr <- quantile(X, gridr)
    gridl[plotN[1]] <- c
    gridr[1] <- c
  }
  Estl <- lpdensity(data = X[X <= c], grid = gridl, bw = hl, 
                    p = p, q = q, v = 1, kernel = kernel, scale = scalel)
  Estr <- lpdensity(data = X[X >= c], grid = gridr, bw = hr, 
                    p = p, q = q, v = 1, kernel = kernel, scale = scaler)
  Estplot <- lpdensity.plot(Estl, Estr, alpha = alpha, type = type, 
                            CItype = CItype, title = title, xlabel = xlabel, ylabel = ylabel, 
                            lty = lty, lwd = lwd, lcol = lcol, pty = pty, pwd = pwd, 
                            pcol = pcol, CIshade = CIshade, CIcol = CIcol, legendTitle = legendTitle, 
                            legendGroups = legendGroups) +  geom_vline(xintercept = ctpt,linetype="dotted")+ 
    labs(title = plot_title, x="Duration of Charge (hours)", y=  expression("Density of Sessions"))+
    scale_x_continuous(limits = c(1, 7), breaks= c(0:7))+
    scale_y_continuous(limits = c(0.0,0.45), breaks = c(0.0,0.1,0.2,0.3,0.4))+
    theme_bw()+ theme(legend.position = "none")+geom_point() +
    theme(axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          plot.title = element_text(size=14 ,hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank() ) +
    theme(text=element_text(family="Helvetica"))
  print(Estplot)
  return(list(Estl = Estl, Estr = Estr, Estplot = Estplot))
}

# Regression discontinuity for announcement effect.
library(rdd)
library(splines)

# Quantile regression for pricing effect.
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(modelsummary)
library(stargazer)
library(lubridate)
library(estimatr)

##### SECTION 1: DATA CLEANING & PRE-PROCESSING

# Post-pricing starts on 26 June, 2018.
# Divide pre-pricing and post-pricing periods.
EV_data$pricing = as.Date(EV_data$startTime,"%d/%m/%y") >= '2018-06-26'

# The 17255th session in the dataset has a duration of more than 18 days.
# Remove session 17255.
EV_data <- subset(EV_data, sessionID != 17255)

# Large spike at shortest 5% of data (less than 6 minutes).
# These sessions represent errors in station operation - remove sessions.
hours <- EV_data$durationHours
hours_low <- quantile(hours, probs = seq(0,1,0.05))
lowest <- hours_low[2]
EV_data <- subset(EV_data, durationHours > lowest)
# Generate days of week column
EV_data$dayOfWeek <- weekdays(as.Date(EV_data$startTime, '%d/%m/%y'))

##### SECTION 2: DISTRIBUTION OF SESSIONS AND SORTING BEHAVIOR

# Subset Announcement period.
EV_announ <- subset(EV_data, as.Date(startTime,"%d/%m/%y") >= '2015-04-01' &
                      as.Date(startTime,"%d/%m/%y") < '2018-06-26')

# Subset Post-Pricing period.
EV_postpr <- subset(EV_data, as.Date(startTime,"%d/%m/%y") >= '2018-06-26')

## GENERATE HISTOGRAMS SHOWING DISTRIBUTION OF SESSION DURATION (Figure 1).

# Announcement period histogram (Figure 1 left).
ggplot(EV_announ, aes(x=durationHours)) + 
  geom_histogram(binwidth=0.2,
                 color = 'black',
                 fill = 'light grey') +
  scale_x_continuous(breaks = seq(0, 14, by = 2), limits = c(0,14)) +
  scale_y_continuous(limits = c(0,800)) +
  geom_vline(xintercept = 4, linetype = 'dashed') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Helvetica", size=28)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="a) Announcement Period",x="Duration of Charge (hours)", y = "Count of Sessions")

# Post-pricing period histogram (Figure 1 right).
ggplot(EV_postpr, aes(x=durationHours)) + 
  geom_histogram(binwidth=0.2,
                 color = 'black',
                 fill = 'light grey') +
  scale_x_continuous(breaks = seq(0, 14, by = 2), limits = c(0,14)) +
  scale_y_continuous(limits = c(0,800)) +
  geom_vline(xintercept = 4, linetype = 'dashed') +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Helvetica", size=28)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="b) Pricing Period",x="Duration of Charge (hours)", y = "Count of Sessions")

## MCCRARY TEST FOR SORTING BEHAVIOR AT 4 HOURS (c = 4) (Supp. Figure 6).

# Announcement period (c = 4 hours) (Supp. Figure 6 left).
# Change p to test different order equations. Change c to test placebos.
RDD_announ <- rddensity(EV_announ$durationHours, c = 4, p = 3)
summary(RDD_announ)
ctpt = 4
rdplotdensity2(RDD_announ, EV_announ$durationHours, 'McCrary Density Test for Announcement Period (April 2015 - June 2018)',
              xlabel='Charging Duration (hours)', 
              ylabel='Density')

# Post-pricing period (c = 4 hours) (Supp. Figure 6 right).
# Change p to test different order equations. Change c to test placebos.
RDD_postpr <- rddensity(EV_postpr$durationHours, c = 4, p = 3)
summary(RDD_postpr)
ctpt = 4
rdplotdensity2(RDD_postpr, EV_postpr$durationHours, 'McCrary Density Test for Post-Pricing Period (June 2018 - Dec 2019)',
              xlabel='Charging Duration (hours)', 
              ylabel='Density')

##### SECTION 3: REGRESSION DISCONTINUITY

# Announcement period only. Post-pricing period not suitable for RD.

## PRE-PROCESSING

# Remove all one-time users because differences cannot be calculated with them.
users_freq <- data.frame(table(EV_data$userID))
freq_users <- subset(users_freq, Freq != 1)
freq_data <- subset(EV_data, userID %in% freq_users$Var1)

# Subset all sessions in Announcement period.
freq_announ <- subset(freq_data, as.Date(startTime,"%d/%m/%y") >= '2015-04-01' &
                        as.Date(startTime,"%d/%m/%y") < '2018-06-26')
users_announ <- data.frame(table(freq_announ$userID))
freq_users_announ <- subset(users_announ, Freq != 1)
freq_announ <- subset(freq_announ, userID %in% freq_users_announ$Var1)

# Calculate differences between prior and current session energy consumption.
# Include all day of week and month of year dummies of current session.
users <- unique(freq_announ$userID)
diffs <- matrix(ncol = 7)
for (user in users) {
  user_data <- subset(freq_announ, userID == user)[,c('startTime', 'Energy','durationHours','parkingLot','stationID',
                                                      'startMonth','dayOfWeek')]
  cons_diff <- diff(user_data$Energy) # Dependent var = difference in energy use between previous and current session
  start_tim <- user_data$startTime[-nrow(user_data)]
  durations <- user_data$durationHours[-nrow(user_data)]
  parkinglt <- user_data$parkingLot[-nrow(user_data)]
  stationID <- user_data$stationID[-nrow(user_data)]
  doW <- user_data$dayOfWeek[-nrow(user_data)]
  moY <- user_data$startMonth[-nrow(user_data)]
  differenc <- cbind(durations, cons_diff, parkinglt, stationID, 
                     doW, moY, start_tim) #include extra columns for more variables
  diffs <- rbind(diffs, differenc)
}
diffs <- diffs[-1,]

# Take log of differences.
min_diff <- ceiling(abs(min(as.numeric(diffs[,2]))))
diffs[,2] <- as.numeric(diffs[,2]) + min_diff
diffs[,2] <- log(as.numeric(diffs[,2]))
diff_data <- as.data.frame(diffs)

# Create time dummy variables
wk <- diff_data$doW
mo <- diff_data$moY
wk_dummy <- model.matrix(~ wk - 1, data = diff_data)
mo_dummy <- model.matrix(~ mo - 1, data = diff_data)
class(wk_dummy) <- "numeric"
class(mo_dummy) <- "numeric"
time_dummies <- cbind(wk_dummy,mo_dummy)
diff_data <- cbind(diff_data, time_dummies)

## STATIC REGRESSION DISCONTINUITY.

# Set up variables for equation.
y <- as.numeric(diff_data$cons_diff) # Dependent variable
x <- as.numeric(diff_data$durations) # Running variable
c1 <- diff_data$parkinglt # Clustering variable
c2 <- diff_data$stationID # Alternate clustering variable
# Month dummies.
m1 <- diff_data$mo1; m2 <- diff_data$mo2; m3 <- diff_data$mo3; m4 <- diff_data$mo4
m5 <- diff_data$mo5; m6 <- diff_data$mo6; m7 <- diff_data$mo7; m8 <- diff_data$mo8
m9 <- diff_data$mo9; m10 <- diff_data$mo10; m11 <- diff_data$mo11; m12 <- diff_data$mo12
# Day of week dummies.
d1 <- diff_data$wkMonday; d2 <- diff_data$wkTuesday; d3 <- diff_data$wkWednesday; d4 <- diff_data$wkThursday
d5 <- diff_data$wkFriday; d6 <- diff_data$wkSaturday; d7 <- diff_data$wkSunday

# Generate polynomial spline.
b <- bs(x, knots = c(4))
b1 <- b[,1]
b2 <- b[,2]
b3 <- b[,3]
b4 <- b[,4]

# Cutpoint = 4 hours (change for placebo test).
ctpt = 4

# Generate Imbens-Kalyanaraman bandwidths.
# kernel = 'epanechnikov' or 'triangular'. Obtain results for both to check robustness.
bw <- with(diff_data, IKbandwidth(x, y, cutpoint = ctpt, kernel = 'epanechnikov'))

# RD clustered by parking lot (c1). To cluster by station ID, change to cluster = c2
rd_diff_c1 <- RDestimate(y ~ x, data = diff_data, cutpoint = ctpt, bw = bw, cluster = c1)
summary(rd_diff_c1)

# RD clustered by parking lot (c1) with spline.
rd_diff_c1_sp <- RDestimate(y ~ x | b1+b2+b3+b4, data = diff_data, cutpoint = ctpt,
                            bw = bw, cluster = c1)
summary(rd_diff_c1_sp)

# RD clustered by parking lot (c1) with spline and time (day of week, month) covariates (Supp. Table 2).
rd_diff_c1_t <- RDestimate(y ~ x | b1+b2+b3+b4 + 
                             m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+
                             d1+d2+d3+d4+d5+d6+d7, 
                           data = diff_data, cutpoint = ctpt, bw = bw, cluster = c1)
summary(rd_diff_c1_t)

# Visualize results.
diff_data$over4 = diff_data$durations > 4
ggplot(diff_data, aes(x = durations, y = cons_diff, color = over4)) +
  geom_point(alpha = 0.7, stroke = 0.4, size=0.01) +
  geom_vline(xintercept=4, linetype="dashed", size = 0.5) + 
  geom_hline(yintercept=log(min_diff), linetype="dashed", size = 0.3, colour = "grey50") +
  coord_cartesian(ylim=c(3.8,4.1),xlim=c(0, 8)) +
  stat_smooth(method = "loess", formula = y~x, fill= "grey15", size = 0.65) +
  scale_colour_manual(values = c("grey30", "firebrick1")) +
  labs(title = "", x="Duration of Charge (hours)", y=  expression("Log Change in Energy Consumption")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none", axis.text=element_text(size=14),
        axis.title=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Regression Discontinuity Analysis, Announcement Period') +
  theme(text=element_text(family="Times New Roman", size=16))

## DYNAMIC RD: EXPANDING WINDOW OVER ANNOUNCEMENT PERIOD

# Reorder sessions based on startTime.
diff_data <- diff_data[order(as.Date(diff_data$start_tim, '%d/%m/%Y')),]
row.names(diff_data) <- NULL
dynamic_diff <- matrix(ncol=4) # Set up matrix to store dynamic RD results
numbers <- (1:nrow(diff_data))[(1:nrow(diff_data) %% 100) == 0] # Define window (increasing at 100 increments)
bw_dyn <- with(diff_data, IKbandwidth(x, y, cutpoint = 4, kernel = 'epanechnikov'))
for (i in numbers) {
  dynamic_data <- diff_data[1:i,] # Loop through windows, expanding by 100
  y_dyn <- as.numeric(dynamic_data$cons_diff) # Dependent variable
  x_dyn <- as.numeric(dynamic_data$durations) # Running variable
  c1_dyn <- dynamic_data$parkinglt # Cluster by parking lot
  b1_dyn <- b1[1:i]; b2_dyn <- b2[1:i]; b3_dyn <- b3[1:i]; b4_dyn <- b4[1:i] # spline
  # Month dummies for sessions within the time window
  m1 <- dynamic_data$mo1; m2 <- dynamic_data$mo2; m3 <- dynamic_data$mo3; m4 <- dynamic_data$mo4
  m5 <- dynamic_data$mo5; m6 <- dynamic_data$mo6; m7 <- dynamic_data$mo7; m8 <- dynamic_data$mo8
  m9 <- dynamic_data$mo9; m10 <- dynamic_data$mo10; m11 <- dynamic_data$mo11; m12 <- dynamic_data$mo12
  # Day of week dummies for sessions the time window
  d1 <- dynamic_data$wkMonday; d2 <- dynamic_data$wkTuesday; d3 <- dynamic_data$wkWednesday; d4 <- dynamic_data$wkThursday
  d5 <- dynamic_data$wkFriday; d6 <- dynamic_data$wkSaturday; d7 <- dynamic_data$wkSunday
  # Use Imbens-Kalyanaraman bandwidth
  bw_dyn <- with(dynamic_data, IKbandwidth(x_dyn, y_dyn, cutpoint = 4, kernel = 'epanechnikov'))
  # Obtain dynamic RD results for the time window
  dynamic_rd <- RDestimate(y_dyn ~ x_dyn | b1_dyn+b2_dyn+b3_dyn+b4_dyn + 
                             m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+
                             d1+d2+d3+d4+d5+d6+d7,
                           data = dynamic_data, cutpoint = 4, bw = bw_dyn,
                           cluster = c1_dyn)
  # Record RD estimate + upper and lower bound
  dyn_row <- c(i,dynamic_rd$est[1], dynamic_rd$est[1]-dynamic_rd$se[1], dynamic_rd$est[1]+dynamic_rd$se[1])
  # Append results as row to dynamic RD results table
  dynamic_diff <- rbind(dynamic_diff, dyn_row)
  print(i) # Record progress
}
dynamic_diff <- dynamic_diff[-1,]
colnames(dynamic_diff) <- c('no obs','est','lower cf','upper cf')
df_dynamicdiff <- as.data.frame(dynamic_diff)

# Visualize results (Supp. Figure 8).
ggplot(df_dynamicdiff, aes(df_dynamicdiff[,1], df_dynamicdiff[,2])) +
  #Title and axes labels
  labs(title = 'Dynamic Regression-Discontinuity Analysis, Announcement',
       y = expression('Estimate of Treatment Effect (%)'),
       x = expression('Cumulative Usage (Number of Sessions)')) +
  #Axes scale
  scale_x_continuous(breaks = c(0,2000,4000,6000,8000,10000)) + 
  scale_y_continuous(limits = c(-0.05,0.21), breaks = c(-0.05,0.00,0.05,0.10,0.15,0.20)) +
  #Confidence interval
  geom_ribbon(aes(ymin=df_dynamicdiff[,3], ymax=df_dynamicdiff[,4]), alpha=0.2) +
  #Horizontal 0 line
  geom_hline(yintercept = 0, colour = 'grey70') +
  #RD estimate line
  geom_line(aes(y = df_dynamicdiff[,2], col = "estimate"),colour="black")+
  #Vertical lines indicating time boundaries of Announcement period.
  geom_vline(xintercept = 0, linetype= 'dotted', colour = 'grey50') +
  geom_vline(xintercept = 11700, linetype= 'dotted', colour = 'grey50') +
  #Theme of plot
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #Center plot title
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(family="Helvetica", size=14))

##### SECTION 4: QUANTILE REGRESSION

## PRE-PROCESSING

# Subset recurring users (155 users).
# Recurring users span both pre- and post-pricing periods.
EV_users <- aggregate(x = EV_data$pricing, by = list(EV_data$userID), FUN = mean)
both_sessions <- subset(EV_users, EV_users$x > 0 & EV_users$x < 1)
EV_recurring <- subset(EV_data, EV_data$userID %in% both_sessions$Group.1)

# Generate day of week categorical variables for recurring users (if not done already).
EV_recurring$dayOfWeek <- weekdays(as.Date(EV_recurring$startTime, '%d/%m/%y'))

## OVERALL AVERAGE EFFECTS FOR ALL RECURRING USERS

pooled_all <- plm(formula = durationHours ~ pricing + factor(parkingLot) + 
                   factor(dayOfWeek) + factor(startMonth), 
                 data=EV_recurring, model="pooling", index = c('parkingLot'))
coeftest(pooled_all, vcov=vcovHC(pooled_all, type="HC1", cluster="group"))
summary(pooled_all)

## CONDITIONAL QUANTILE REGRESSION, QUINTILES BASED ON BASELINE (PRE-PRICING) USAGE

# Define recurring user quintiles based on pre-pricing usage.
pre_pricing <- subset(EV_recurring, pricing == FALSE)
users <- unique(pre_pricing$userID)
user_quantiles <- matrix(ncol = 2)
for (user in users) {
  user_dat <- subset(pre_pricing, userID == user)
  avg_dura <- mean(user_dat$durationHours)
  user_avg_dura <- cbind(user, avg_dura)
  user_quantiles <- rbind(user_quantiles, user_avg_dura)
}
user_quantiles <- user_quantiles[-1,]
user_quantiles <- as.data.frame(user_quantiles)
# Assign users to quintiles.
user_quantiles$quantile_dura <- ntile(user_quantiles$avg_dura, 5)
# Add quintile column to EV recurring users dataset.
EV_recurring$QuantileDuration <- user_quantiles$quantile_dura[match(EV_recurring$userID, user_quantiles$user)]

# Subset quintiles.
quan_01 <- subset(EV_recurring, QuantileDuration == 1)
quan_02 <- subset(EV_recurring, QuantileDuration == 2)
quan_03 <- subset(EV_recurring, QuantileDuration == 3)
quan_04 <- subset(EV_recurring, QuantileDuration == 4)
quan_05 <- subset(EV_recurring, QuantileDuration == 5)

# Subset low-demand and high-demand users.
low_demand <- subset(EV_recurring, QuantileDuration < 3)
high_demand <- subset(EV_recurring, QuantileDuration >= 3)

# Run quantile regression

dat = quan_05 # Select quintile (01 to 05, low_demand or high_demand) (Supp. Table 3).
pooled_quan <- plm(formula=durationHours ~ pricing + factor(parkingLot) + 
                   factor(dayOfWeek) + factor(startMonth), 
                 data=dat, model="pooling", index = c('parkingLot'))
coeftest(pooled_quan, vcov=vcovHC(pooled_quan, type="HC1", cluster="group"))
summary(pooled_quan)

# Use alternate packages to check result robustness and standard errors.
modelsummary(pooled_quan, statistic_override = vcovHC)
se = vcovHC(pooled_quan, type = "HC1") %>% diag() %>% sqrt()
stargazer(pooled_quan, pooled_quan, 
          se = list(NULL, se), type = 'text')

## CONDITIONAL QUANTILE REGRESSION, DYNAMIC EFFECTS

dat = low_demand # Select quintile (01 to 05, low_demand or high_demand)

# Create matrix of dynamic results: cumulative effects of policy in 1-month intervals.
month_elapse = 1 # Start at 1 month after pricing implementation
filter_date <- strftime(ymd(as.Date('2018-06-26')) %m+% months(month_elapse)) # Define window
dyn_price <- matrix(ncol=5) # Create matrix to store dynamic results
while (filter_date <= as.Date('2020-01-26')) {
  filter_dat = subset(dat, as.Date(startTime,"%d/%m/%y") <= filter_date)
  pooled_dat <- plm(formula=durationHours ~ pricing+ factor(parkingLot) + factor(dayOfWeek) + factor(startMonth), data=filter_dat, model="pooling", index = c('parkingLot'))
  test_dat <- coeftest(pooled_dat, vcov=vcovHC(pooled_dat, type="HC1", cluster="group"))
  row_est <- test_dat[2,1]
  row_CI_low <- test_dat[2,1] - 1.96*test_dat[2,2]
  row_CI_high <- test_dat[2,1] + 1.96*test_dat[2,2]
  dyn_row <- c(month_elapse, row_est, row_CI_low, row_CI_high, nrow(subset(filter_dat, pricing == TRUE)))
  dyn_price <- rbind(dyn_price, dyn_row)
  month_elapse = month_elapse + 1
  filter_date <- strftime(ymd(as.Date('2018-06-26')) %m+% months(month_elapse))
  print(filter_date)
}
dyn_price <- dyn_price[-1,]
colnames(dyn_price) <- c('months_elapsed','est','lower cf','upper cf','post_sessions')
df_dyn_price <- as.data.frame(dyn_price)

# Visualize results (Supp. Figure 9).
ggplot(df_dyn_price, aes(df_dyn_price[,1], df_dyn_price[,2])) +
  #Title and axes labels
  labs(title = 'Price effects (expanding window), All Recurring Users',
       y = expression('Change in charging duration (hr)'),
       x = expression('Months elapsed since pricing implementation')) +
  #Axes scale
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) + 
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) +
  #Confidence interval
  geom_ribbon(aes(ymin=df_dyn_price[,3], ymax=df_dyn_price[,4]), alpha=0.2) +
  #Horizontal 0 line
  geom_hline(yintercept = 0, colour = 'grey70') +
  #RD estimate line
  geom_line(aes(y = df_dyn_price[,2], col = "estimate"),colour="black")+
  #Theme of plot
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #Center plot title
  theme(plot.title = element_text(hjust = 0.5)) +
  #Font Times New Roman (import extrafont library and load fonts)
  theme(text=element_text(family="Helvetica", size=12))

# Robustness check: gamma function GLM

dat = low_demand # Select quintile (01 to 05, low_demand or high_demand)

gamma <- glm(durationHours ~ pricing + factor(parkingLot) + 
               factor(dayOfWeek) + factor(startMonth), 
             family=Gamma(link="log"),dat)

summary(gamma, cluster="parkingLot")

## PLACEBO TESTS FOR QUANTILE REGRESSION

# Alternate date as placebo.

dat_placebo <- subset(high_demand, pricing == FALSE) # Subset high-demand users (most impacted by pricing)
# Assign any date as placebo
dat_placebo$placebo <- as.Date(dat_placebo$startTime,"%d/%m/%y") >= '2016-07-01'

# Placebo test formula (Supp. Table 4).
pooled_placebo <- plm(formula=durationHours ~ placebo+ factor(parkingLot) + 
                        factor(dayOfWeek) + factor(startMonth), 
                      data=dat_placebo, model="pooling", index = c('parkingLot'))
coeftest(pooled_placebo, vcov=vcovHC(pooled_placebo, type="HC1", cluster="group"))
mean(subset(dat_placebo, placebo==FALSE)$durationHours)

modelsummary(pooled_placebo, statistic_override = vcovHC)
se = vcovHC(pooled_placebo, type = "HC1") %>% diag() %>% sqrt()
stargazer(pooled_placebo, pooled_placebo, 
          se = list(NULL, se), type = 'text')

# Average Effects of Policy Treatment as Continuous Variable

EV_recurring_pre <- subset(EV_recurring, pricing == FALSE)
EV_recurring_post <- subset(EV_recurring, pricing == TRUE)
EV_avgs_pre <- aggregate(durationHours~userID, EV_recurring_pre, FUN=mean)
EV_avgs_post <- aggregate(durationHours~userID, EV_recurring_post, FUN=mean)
EV_analysis <- merge(EV_avgs_pre, EV_avgs_post, by.x = "userID", by.y = "userID")
colnames(EV_analysis) <- c('userID','duration_pre','duration_post')
diff_pol <- EV_analysis$duration_post - EV_analysis$duration_pre
EV_analysis$diff <- diff_pol

# Obtain count of sessions for each user as covariate
EV_count <- as.data.frame(table(EV_both_sessions$userID))
EV_count_pre <- as.data.frame(table(EV_both_sessions_pre$userID))
EV_count_post <- as.data.frame(table(EV_both_sessions_post$userID))
EV_analysis <- merge(EV_analysis, EV_count, by.x = "userID", by.y = "Var1")

# Obtain dummy variable of whether session surpassed 4 hour price jump
EV_analysis$surpass <- EV_analysis$duration_pre > 4

# Obtain dependent variable (log of difference)
for_ln <- EV_analysis$diff + ceiling(abs(min(EV_analysis$diff)))
EV_analysis$ln_diff <- log(for_ln)

# Analysis (Supp. Table 5)
model_1 <- lm(formula = ln_diff ~ duration_pre + Freq, data=EV_analysis)
coeftest(model_1, vcov = vcovHC(model_1, type="HC1"))
model_1 <- lm_robust(formula = ln_diff ~ duration_pre + Freq, data=EV_analysis, se_type='stata')
summary(model_1)

model_2 <- lm(formula = ln_diff ~ factor(surpass) + Freq, data=EV_analysis)
coeftest(model_2, vcov = vcovHC(model_2, type="HC1"))
model_2 <- lm_robust(formula = ln_diff ~ factor(surpass) + Freq, data=EV_analysis, se_type='stata')
summary(model_2)

model_3 <- lm(formula = ln_diff ~ duration_pre + factor(surpass) + Freq, data=EV_analysis)
coeftest(model_3, vcov = vcovHC(m, type="HC1"))
model_3 <- lm_robust(formula = ln_diff ~ duration_pre + factor(surpass) + Freq, data=EV_analysis, se_type='stata')
summary(model_3)

##### SECTION 5: SUMMARY STATISTICS

# All users (Supp. Table 1, "All Sessions" column 1).

mean(EV_data$durationHours)
sd(EV_data$durationHours)
min(EV_data$durationHours)
max(EV_data$durationHours)
nrow(EV_data)

mean(EV_data$Energy)
sd(EV_data$Energy)
min(EV_data$Energy)
max(EV_data$Energy)

mean(table(EV_data$userID))
sd(table(EV_data$userID))
min(table(EV_data$userID))
max(table(EV_data$userID))
length(table(EV_data$userID))

# Recurring users (Supp. Table 1, "All Sessions" column 2).

mean(EV_recurring$durationHours)
sd(EV_recurring$durationHours)
min(EV_recurring$durationHours)
max(EV_recurring$durationHours)
nrow(EV_recurring)

mean(EV_recurring$Energy)
sd(EV_recurring$Energy)
min(EV_recurring$Energy)
max(EV_recurring$Energy)

mean(table(EV_recurring$userID))
sd(table(EV_recurring$userID))
min(table(EV_recurring$userID))
max(table(EV_recurring$userID))
length(table(EV_recurring$userID))

# Quintiles (Supp. Table 1, "All Sessions" columns 3-7).

dat = quan_05 # change quintile (01 to 05, low_demand, or high_demand)

mean(dat$durationHours)
sd(dat$durationHours)
min(dat$durationHours)
max(dat$durationHours)
nrow(dat)

mean(dat$Energy)
sd(dat$Energy)
min(dat$Energy)
max(dat$Energy)

mean(table(dat$userID))
sd(table(dat$userID))
min(table(dat$userID))
max(table(dat$userID))
length(table(dat$userID))

# Check differences in usage between user types

user_freq <- as.data.frame(table(EV_data$userID))
user_freq_one <- subset(user_freq, Freq == 1)
user_freq_freq <- subset(user_freq, Freq > 1)
one_time <- subset(EV_data, userID %in% user_freq_one$Var1) # users with only one session
more_time <- subset(EV_data, userID %in% user_freq_freq$Var1) # users with more than one session
mean(one_time$durationHours)
mean(one_time$Energy)
mean(more_time$durationHours)
mean(more_time$Energy)
mean(table(more_time$userID))

# t-test if one-time + multiple-time users have same usage.
t.test(one_time$Energy, more_time$Energy)

# Pre-pricing summary statistics, all users (Supp. Table 1, "Pre-Pricing" column 1).

mean(subset(EV_data, pricing==FALSE)$durationHours)
sd(subset(EV_data, pricing==FALSE)$durationHours)
min(subset(EV_data, pricing==FALSE)$durationHours)
max(subset(EV_data, pricing==FALSE)$durationHours)
nrow(subset(EV_data, pricing==FALSE))

mean(subset(EV_data, pricing==FALSE)$Energy)
sd(subset(EV_data, pricing==FALSE)$Energy)
min(subset(EV_data, pricing==FALSE)$Energy)
max(subset(EV_data, pricing==FALSE)$Energy)

mean(table(subset(EV_data,pricing==FALSE)$userID))
sd(table(subset(EV_data,pricing==FALSE)$userID))
min(table(subset(EV_data,pricing==FALSE)$userID))
max(table(subset(EV_data,pricing==FALSE)$userID))
length(table(subset(EV_data,pricing==FALSE)$userID))

# Pre-pricing summary statistics, recurring users (Supp. Table 1, "Pre-Pricing" column 2).

mean(subset(EV_recurring, pricing==FALSE)$durationHours)
sd(subset(EV_recurring, pricing==FALSE)$durationHours)
min(subset(EV_recurring, pricing==FALSE)$durationHours)
max(subset(EV_recurring, pricing==FALSE)$durationHours)
nrow(subset(EV_recurring, pricing==FALSE))

mean(subset(EV_recurring, pricing==FALSE)$Energy)
sd(subset(EV_recurring, pricing==FALSE)$Energy)
min(subset(EV_recurring, pricing==FALSE)$Energy)
max(subset(EV_recurring, pricing==FALSE)$Energy)

mean(table(subset(EV_recurring,pricing==FALSE)$userID))
sd(table(subset(EV_recurring,pricing==FALSE)$userID))
min(table(subset(EV_recurring,pricing==FALSE)$userID))
max(table(subset(EV_recurring,pricing==FALSE)$userID))
length(table(subset(EV_recurring,pricing==FALSE)$userID))

# Pre-pricing summary statistics, quintiles of recurring users (Supp. Table 1, "Pre-Pricing" columns 3-7).

dat = quan_05 # change quintile (01 to 05, low_demand, or high_demand)

mean(subset(dat, pricing==FALSE)$durationHours)
sd(subset(dat, pricing==FALSE)$durationHours)
min(subset(dat, pricing==FALSE)$durationHours)
max(subset(dat, pricing==FALSE)$durationHours)
nrow(subset(dat, pricing==FALSE))

mean(subset(dat, pricing==FALSE)$Energy)
sd(subset(dat, pricing==FALSE)$Energy)
min(subset(dat, pricing==FALSE)$Energy)
max(subset(dat, pricing==FALSE)$Energy)

mean(table(subset(dat,pricing==FALSE)$userID))
sd(table(subset(dat,pricing==FALSE)$userID))
min(table(subset(dat,pricing==FALSE)$userID))
max(table(subset(dat,pricing==FALSE)$userID))
length(table(subset(dat,pricing==FALSE)$userID))

# Post-pricing summary statistics, all users (Supp. Table 1, "Post-Pricing" column 1).

mean(subset(EV_data, pricing==TRUE)$durationHours)
sd(subset(EV_data, pricing==TRUE)$durationHours)
min(subset(EV_data, pricing==TRUE)$durationHours)
max(subset(EV_data, pricing==TRUE)$durationHours)
nrow(subset(EV_data, pricing==TRUE))

mean(subset(EV_data, pricing==TRUE)$Energy)
sd(subset(EV_data, pricing==TRUE)$Energy)
min(subset(EV_data, pricing==TRUE)$Energy)
max(subset(EV_data, pricing==TRUE)$Energy)

mean(table(subset(EV_data,pricing==TRUE)$userID))
sd(table(subset(EV_data,pricing==TRUE)$userID))
min(table(subset(EV_data,pricing==TRUE)$userID))
max(table(subset(EV_data,pricing==TRUE)$userID))
length(table(subset(EV_data,pricing==TRUE)$userID))

# Post-pricing summary statistics, recurring users (Supp. Table 1, "Post-Pricing" column 2).

mean(subset(EV_recurring, pricing==TRUE)$durationHours)
sd(subset(EV_recurring, pricing==TRUE)$durationHours)
min(subset(EV_recurring, pricing==TRUE)$durationHours)
max(subset(EV_recurring, pricing==TRUE)$durationHours)
nrow(subset(EV_recurring, pricing==TRUE))

mean(subset(EV_recurring, pricing==TRUE)$Energy)
sd(subset(EV_recurring, pricing==TRUE)$Energy)
min(subset(EV_recurring, pricing==TRUE)$Energy)
max(subset(EV_recurring, pricing==TRUE)$Energy)

mean(table(subset(EV_recurring,pricing==TRUE)$userID))
sd(table(subset(EV_recurring,pricing==TRUE)$userID))
min(table(subset(EV_recurring,pricing==TRUE)$userID))
max(table(subset(EV_recurring,pricing==TRUE)$userID))
length(table(subset(EV_recurring,pricing==TRUE)$userID))

# Post-pricing summary statistics, quintiles of recurring users (Supp. Table 1, "Post-Pricing" columns 3-7).

dat = quan_05 # change quintile (01 to 05, low_demand, or high_demand)

mean(subset(dat, pricing==TRUE)$durationHours)
sd(subset(dat, pricing==TRUE)$durationHours)
min(subset(dat, pricing==TRUE)$durationHours)
max(subset(dat, pricing==TRUE)$durationHours)
nrow(subset(dat, pricing==TRUE))

mean(subset(dat, pricing==TRUE)$Energy)
sd(subset(dat, pricing==TRUE)$Energy)
min(subset(dat, pricing==TRUE)$Energy)
max(subset(dat, pricing==TRUE)$Energy)

mean(table(subset(dat,pricing==TRUE)$userID))
sd(table(subset(dat,pricing==TRUE)$userID))
min(table(subset(dat,pricing==TRUE)$userID))
max(table(subset(dat,pricing==TRUE)$userID))
length(table(subset(dat,pricing==TRUE)$userID))

# Test if recurring and non-recurring users differ in usage.
mean(EV_recurring$durationHours)
mean(EV_recurring$Energy)
mean(table(EV_recurring$userID))
not_recurring <- subset(more_time, !(userID %in% EV_recurring$userID))
not_recurring_2 <- subset(EV_data, !(userID %in% EV_recurring$userID))
t.test(not_recurring$durationHours, EV_recurring$durationHours)
t.test(not_recurring$Energy, EV_recurring$Energy)
t.test(table(not_recurring$userID), table(EV_recurring$userID))

# Check if low_demand and high_demand users differ in usage (mean and distribution).
mean(subset(low_demand, pricing == FALSE)$durationHours)
mean(subset(high_demand, pricing == FALSE)$durationHours)
mean(subset(low_demand, pricing == TRUE)$durationHours)
mean(subset(high_demand, pricing == TRUE)$durationHours)
t.test(subset(high_demand, pricing == FALSE)$durationHours, subset(high_demand, pricing == TRUE)$durationHours)
t.test(subset(low_demand, pricing == FALSE)$durationHours, subset(high_demand, pricing == FALSE)$durationHours)
t.test(subset(low_demand, pricing == TRUE)$durationHours, subset(high_demand, pricing == TRUE)$durationHours)
t.test(subset(low_demand, pricing == FALSE)$durationHours, subset(low_demand, pricing == TRUE)$durationHours)
ks.test(subset(low_demand, pricing == FALSE)$durationHours, subset(high_demand, pricing == FALSE)$durationHours)

## Histograms based on user type (low_demand vs. high_demand) (Supp. Figure 5).

low_demand_graph <- low_demand
low_demand_graph$pricing <-replace(low_demand_graph$pricing, low_demand_graph$pricing == FALSE, 'Pre-Pricing')
low_demand_graph$pricing <-replace(low_demand_graph$pricing, low_demand_graph$pricing == TRUE, 'Post-Pricing')
low_demand_graph$pricing_f = factor(low_demand_graph$pricing, levels=c('Pre-Pricing','Post-Pricing'))

low_demand_graph$graph <- low_demand_graph$pricing
low_demand_graph$graph <- replace(low_demand_graph$graph, 
                              ((as.Date(low_demand_graph$startTime,"%d/%m/%y") < '2018-06-26') &
                                 (as.Date(low_demand_graph$startTime,"%d/%m/%y") >= '2015-04-01')),
                              'Announcement')
low_demand_graph2 <- subset(low_demand_graph, (graph == 'Announcement') | (graph == 'Post-Pricing') )

ggplot(as.data.frame(low_demand_graph), aes(x=durationHours)) + 
  geom_histogram(color = 'white', fill='lightcoral', binwidth = 0.1, size = 0.1)+
  facet_grid(pricing_f~.) +
  #scale_x_continuous(breaks = seq(0,10 , by = 1))+
  xlim(4,11)+
  ylim(0,20) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Low-Demand Users') +
  ylab('Frequency') +
  xlab('Duration of Charge (hrs)') + 
  theme(text=element_text(family='Helvetica', size=30)) +
  theme(axis.text = element_text(size = 18)) +
  theme(plot.title = element_text(hjust = 0.5))

high_demand_graph <- high_demand
high_demand_graph$pricing <-replace(high_demand_graph$pricing, high_demand_graph$pricing == FALSE, 'Pre-Pricing')
high_demand_graph$pricing <-replace(high_demand_graph$pricing, high_demand_graph$pricing == TRUE, 'Post-Pricing')
high_demand_graph$pricing_f = factor(high_demand_graph$pricing, levels=c('Pre-Pricing','Post-Pricing'))

high_demand_graph$graph <- high_demand_graph$pricing
high_demand_graph$graph <- replace(high_demand_graph$graph, 
                                ((as.Date(high_demand_graph$startTime,"%d/%m/%y") < '2018-06-26') &
                                   (as.Date(high_demand_graph$startTime,"%d/%m/%y") >= '2015-04-01')),
                                'Announcement')
high_demand_graph2 <- subset(high_demand_graph, (graph == 'Announcement') | (graph == 'Post-Pricing') )

ggplot(as.data.frame(high_demand_graph), aes(x=durationHours)) + 
  geom_histogram(color = 'white', fill='gray70', binwidth = 0.2, size = 0.1)+
  facet_grid(graph~.) +
  scale_x_continuous(breaks = seq(0,25 , by = 2))+
  ylim(0,300) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('High-Demand Users') +
  ylab('Frequency') +
  xlab('Duration of Charge (hrs)') + 
  theme(text=element_text(family='Helvetica', size=30)) +
  theme(axis.text = element_text(size = 18)) +
  theme(plot.title = element_text(hjust = 0.5))

## Histogram of plug-in and plug-out hours (Supp. Figure 4).

hours <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
hours <- as.data.frame(hours)

# Pre-pricing all (Supp. Figure 4 upper left)
EV_data_pre <- subset(EV_data, pricing == FALSE)
pre_all_starts <- as.data.frame(table(EV_data_pre$startHour))
pre_all_ends <- as.data.frame(table(EV_data_pre$endHour))
pre_all_dist = merge(x=hours, y=pre_all_starts, by.x = 'hours', by.y = 'Var1')
pre_all_dist = merge(x=pre_all_dist, y=pre_all_ends, by.x = 'hours', by.y = 'Var1')
names(pre_all_dist) <- c('Hour','plug-in','plug-out')
pre_all_melt <- melt(pre_all_dist, id='Hour')
ggplot(pre_all_melt,aes(x=Hour,y=value,fill=variable)) + 
  geom_bar(stat='identity', position='identity', alpha=.5) +
  scale_fill_manual(values=c("gray70","lightcoral")) +
  scale_x_continuous(breaks = seq(0,23 , by = 1))+
  ylim(0,2000) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('All Users, Pre-Pricing') +
  ylab('Frequency') +
  xlab('') + 
  theme(text=element_text(family='Helvetica', size=18)) +
  theme(axis.text = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5))

# Post-pricing all (Supp. Figure 4 upper right)
EV_data_post <- subset(EV_data, pricing == TRUE)
post_all_starts <- as.data.frame(table(EV_data_post$startHour))
post_all_ends <- as.data.frame(table(EV_data_post$endHour))
post_all_dist = merge(x=hours, y=post_all_starts, by.x = 'hours', by.y = 'Var1')
post_all_dist = merge(x=post_all_dist, y=post_all_ends, by.x = 'hours', by.y = 'Var1')
names(post_all_dist) <- c('Hour','plug-in','plug-out')
post_all_melt <- melt(post_all_dist, id='Hour')
ggplot(post_all_melt,aes(x=Hour,y=value,fill=variable)) + 
  geom_bar(stat='identity', position='identity', alpha=.5) +
  scale_fill_manual(values=c("gray70","lightcoral")) +
  scale_x_continuous(breaks = seq(0,23 , by = 1))+
  ylim(0,2000) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('All Users, Post-Pricing') +
  ylab('') +
  xlab('') + 
  theme(text=element_text(family='Helvetica', size=18)) +
  theme(axis.text = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5))

# Pre-pricing recurring only (Supp. Figure 4 lower left)
EV_recurring_pre <- subset(EV_recurring, pricing == FALSE)
pre_starts <- as.data.frame(table(EV_recurring_pre$startHour))
pre_ends <- as.data.frame(table(EV_recurring_pre$endHour))
pre_dist = merge(x=hours, y=pre_starts, by.x = 'hours', by.y = 'Var1')
pre_dist = merge(x=pre_dist, y=pre_ends, by.x = 'hours', by.y = 'Var1')
names(pre_dist) <- c('Hour','plug-in','plug-out')
pre_melt <- melt(pre_dist, id='Hour')
ggplot(pre_melt,aes(x=Hour,y=value,fill=variable)) + 
  geom_bar(stat='identity', position='identity', alpha=.5) +
  scale_fill_manual(values=c("gray70","lightcoral")) +
  scale_x_continuous(breaks = seq(0,23 , by = 1))+
  ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Long-Term Users, Pre-Pricing') +
  ylab('Frequency') +
  xlab('Time of Day (Hour)') + 
  theme(text=element_text(family='Helvetica', size=18)) +
  theme(axis.text = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5))

# Post-pricing recurring only (Supp. Figure 4 lower right)
EV_recurring_post <- subset(EV_recurring, pricing == TRUE)
post_starts <- table(EV_recurring_post$startHour)
post_ends <- table(EV_recurring_post$endHour)
post_dist = merge(x=hours, y=post_starts, by.x = 'hours', by.y = 'Var1')
post_dist = merge(x=post_dist, y=post_ends, by.x = 'hours', by.y = 'Var1')
names(post_dist) <- c('Hour','plug-in','plug-out')
post_melt <- melt(post_dist, id='Hour')
ggplot(post_melt,aes(x=Hour,y=value,fill=variable)) + 
  geom_bar(stat='identity', position='identity', alpha=.5) +
  scale_fill_manual(values=c("gray70","lightcoral")) +
  scale_x_continuous(breaks = seq(0,23 , by = 1))+
  ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Long-Term Users, Post-Pricing') +
  ylab('') +
  xlab('Time of Day (Hour)') + 
  theme(text=element_text(family='Helvetica', size=18)) +
  theme(axis.text = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5))

## Imbens-Kalyanaraman optimal bandwidth charts (Supp. Figure 7)

bndwdth<-c(0.25, 0.5,  0.75,  1,  1.5, 2, 2.5, 3)
mat <- data.frame(ncol=4, nrow=8)
for( i in 1:8){
  bw <- bndwdth[i]*with(diff_data, IKbandwidth(x, y, cutpoint = 4, kernel = 'epanechnikov'))
  RRDD <- RDestimate(y ~ x, data = diff_data, cutpoint = 4, bw = bw, cluster = c1)
  mat[i,1]<-bndwdth[i]
  mat[i,2]<-RRDD$ci[1]
  mat[i,3]<-RRDD$ci[4]
  mat[i,4]<-RRDD$est[1]
}

ggplot(mat, aes(mat[,1])) + 
  labs(title = "Epanechnikov Kernel", x="% of Imbens and Kalyanaraman Optimal Bandwidth", y=  expression("Estimate of Coefficient"))+
  geom_ribbon(aes(ymin=mat[,2], ymax=mat[,3]), alpha=0.2) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey70') +
  geom_line(aes(y = mat[,4], col = "estimate"),colour="black")+
  scale_x_continuous(breaks = seq(0,3 , by = 0.5))+
  scale_y_continuous(limits=c(-0.05,0.05))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=14,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(text=element_text(family="Helvetica", size=12))
