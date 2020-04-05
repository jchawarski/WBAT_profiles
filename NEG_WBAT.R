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

NEG_CTD <- read.ctd.sbe("CTD\\Dat\\26da.2017.9.8.cnv")  #read a single CTD files that coincides with WBAT profile
NEG_CTD.data <- data.frame(NEG_CTD[["data"]])  #select just the data portion

# # # NEG CTD data only records timeS (time in seconds) in the time slot, therfore we use startime as a reference

time <- NEG_CTD[["metadata"]]$startTime   # collect the start time of the cast
time_seconds <- period_to_seconds(lubridate::seconds(time))  # convert start time to seconds (timeY)

NEG_CTD.data$Time <- NEG_CTD.data$timeS + time_seconds # convert timeS to seconds by adding to start time

# Read in the WBAT data

Stn3_wbat <- read.csv("WBAT\\Export\\Stn3_wbat_200kHz_sv.csv")

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






