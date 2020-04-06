# WBAT/CTD profile synchronization and plotting

## Developed by Julek Chawarski [ julek.chawarski@mi.mun.ca ]
## Initial commit: April-05-2020

### This script was developed to facilitate the synchronization of WBAT profiles with depth 
### and other environmental parameters collected by the CTD-rosette. WBAT profiles are 
### processed and cleaned in Echoview 10. 

# load packages
library(ggplot2)
library(dplyr)
library(oce)
library(reshape2)
library(lubridate)

setwd("C:\\Users\\jchawars\\OneDrive - Memorial University of Newfoundland\\NEG")

NEG_CTD <- read.ctd.sbe("CTD\\UpDown\\26da.2017.9.8_110618.cnv")  #read a single CTD files that coincides with WBAT profile
NEG_CTD.data <- data.frame(NEG_CTD[["data"]])  #select just the data portion

# # # NEG CTD data only records timeS (time in seconds) in the time slot, therfore we use startime as a reference

time <- NEG_CTD[["metadata"]]$startTime   # collect the start time of the cast
time_seconds <- period_to_seconds(lubridate::seconds(time))  # convert start time to seconds (timeY)

NEG_CTD.data$Time <- NEG_CTD.data$timeS + time_seconds # convert timeS to seconds by adding to start time

# Read in the WBAT data

Stn3_wbat <- read.csv("WBAT\\Export\\Stn3_wbat_200kHz_sv.csv")
Stn3_wbat <- read.csv("WBAT\\Export\\Stn3_wbat_70kHz_sv.csv")

Stn3_wbat$Date_M <- as.Date(as.character(Stn3_wbat$Date_M), "%Y%m%d")  # convert Date_M to class "Date"
Stn3_wbat$datetime <- as.POSIXct(paste(Stn3_wbat$Date_M, Stn3_wbat$Time_M), format="%Y-%m-%d %H:%M:%S", tz="UTC")

Stn3_wbat$Time <- period_to_seconds(lubridate::seconds(Stn3_wbat$datetime))

Stn3_wbat$Time <- Stn3_wbat$Time - 116    # <--- Include the correction factor for time difference between CTD and WBAT 
# For example - Stn3 (GearNo 8) has WBAT ahead of CTD by 1:54 == 116 seconds

Stn3_wbat$Time <- as.character(Stn3_wbat$Time) #convert to character vector for joining function
NEG_CTD.data$Time <- as.character(round(NEG_CTD.data$Time, digits=0))

wbat.ctd <- Stn3_wbat %>% left_join(., NEG_CTD.data, by="Time")


# --- In Progress --- 
# calculate absorption coefficient to adjust Sv and TS values based on sound speed

# Following Doonan et al 2003 (ICES) -- ignored boric acid component of equation

Temp <- 0.35
Sal <- 34.54
Depth <- 50

# Following Doonan et al 2003 (ICES) -- ignored boric acid component of equation
c  <- 1412 + 3.21*Temp +1.19*Sal + 0.0167*Depth   # soundSpeed
f  <- 38                                          # Acoustic Frequency (kHz)
A2 <- 22.19*Sal*(1+0.0017*Temp)
f2 <- 1.8e7*exp(-1818/(Temp+273.1))
P2 <- exp(-1.76e4*Depth)
A3 <- 4.937e-4 - 2.59e-5*Temp + 9.11e-7*Temp^2 - 1.5e-8*Temp^3
P3 <- 1-3.83e-5*Depth + 4.9e-10*Depth^2

alpha <- A2*P2*f2*f^2/(f2^2+f)/c + A3*P3*f^2

# Following Francois & Garrison 1982

c <_ 1448.96
pH <- 8

A1 <- 8.86*10^(0.78*pH-5)
A2 <- (21.44*Sal*(1+ (0.025*Temp)))/c
A3 <- 4.937e-4 - 2.59e-5*Temp + 9.11e-7*Temp^2 - (1.5e-8*Temp^3)
f1 <- 2.8*sqrt(Sal/35)*10^(4-(1245/(Temp+273)))
f2 <- (8.17*10^(8-(1990/Temp+273)))/1+(0.0018*(Sal-35))
P2 <- 1-1.37e-4*Depth + 6.2e-9*Depth^2
P3 <- 1-3.83e-5*Depth + 4.9e-10*Depth^2

alpha <- f^2 * ((A1*f1/f1^2+f^2) + (A2*P2*f2/f2^2+f^2) + (A3*P3))


# --- In Progress --- 


# create profile by custom binwidth

d <- wbat.ctd$depth
N <- wbat.ctd$Sv_mean # calibrated Sv_values
db <- seq(0, max(d, na.rm=T), 0.5) #creates a depth range vector at x resolution (start range, end range, bin size[x])
wbat.binned <- binMean1D(d, N, db) # supplies the mean value for each depth bin
wbat.binned<- data.frame(c(wbat.binned["xmids"],wbat.binned["result"]))

d_scale <- seq(0, max(wbat.ctd$depth, na.rm=T), 200)

ggplot(wbat.binned, aes(x=xmids, y=result)) + geom_line(size=1) +
  coord_flip() + scale_x_reverse(breaks=d_scale) +
  labs(x="Depth (m)", y="Sv mean") +       # 
  scale_y_continuous(breaks=c(-85,-80,-75,-70,-65,-60)) +                   #standardized scale for turbidty
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.key=element_blank(),
        legend.key.size = unit(1, "line"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size =1), 
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), size=12), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), size=12),
        axis.ticks.length=unit(-0.25, "cm"))






