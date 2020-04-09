# WBAT/CTD profile synchronization and plotting

## Developed by Julek Chawarski [ julek.chawarski@mi.mun.ca ]
## Initial commit: April-05-2020

### This script was developed to facilitate the synchronization of WBAT profiles with depth 
### and other environmental parameters collected by the CTD-rosette. WBAT profiles are 
### processed and cleaned using impulse noise removal in Echoview 10. WBAT profiles are exported
### as .csv by individual pings (Cell export 1 ping X 100 m bins).

# # # NEG2017 CTD data only records timeS (time in seconds) in the time slot, therfore we use startime as a reference
# # # Process varies slightly for CTDs with TimeY and TimeJ variables

library(ggplot2)
library(dplyr)
library(oce)
library(reshape2)
library(lubridate)

setwd("C:\\Users\\jchawars\\OneDrive - Memorial University of Newfoundland\\NEG")

CTD <- read.ctd.sbe("CTD\\UpDown\\26da.2017.9.8_110618.cnv")  #read a single CTD files that coincides with WBAT profile
CTD.data <- data.frame(CTD[["data"]])  #select just the data portion

time <- CTD[["metadata"]]$startTime   # collect the start time of the cast
time_seconds <- period_to_seconds(lubridate::seconds(time))  # convert start time to seconds (timeY)

CTD.data$Time <- CTD.data$timeS + time_seconds # convert timeS to seconds by adding to start time

# Read in the WBAT data

wbat <- read.csv("WBAT\\Export\\Stn3_wbat_200kHz_sv.csv")
wbat <- read.csv("WBAT\\Export\\Stn3_wbat_70kHz_sv.csv")

wbat$Date_M <- as.Date(as.character(wbat$Date_M), "%Y%m%d")  # convert Date_M to class "Date"
wbat$datetime <- as.POSIXct(paste(wbat$Date_M, wbat$Time_M), format="%Y-%m-%d %H:%M:%S", tz="UTC")
wbat$Time <- period_to_seconds(lubridate::seconds(wbat$datetime))

wbat$Time <- wbat$Time - 116                                                 # <--- Include the correction factor for time difference between CTD and WBAT 
                                                                              # For example - Stn3 (GearNo 8) has WBAT ahead of CTD by 1:54 == 116 seconds

wbat$Time <- as.character(wbat$Time) #convert to character vector for joining function
CTD.data$Time <- as.character(round(CTD.data$Time, digits=0))
wbat.ctd <- wbat %>% left_join(., CTD.data, by="Time") # join ctd and wbat data.frames by matching time cases

# -- Correction factor for Sv equation to apply absorption coefficient -- #

# echoview default absorption coefficient - below values are derived from ecs file
# swSoundsAbsorption(f[Hz], Salinity, Temperature, Pressure, pH, formualtion= "francois-garrison")

abs.default <-  swSoundAbsorption(200000, 35, 8, 0, 8, formulation= "francois-garrison") #200 kHz
abs.default <-  swSoundAbsorption(70000, 35, 8, 0, 8, formulation= "francois-garrison") # 70 kHz

wbat.ctd$Corr_factor <- abs.default - wbat.ctd$Sv_mean 

wbat.ctd$abs_coeff <- swSoundAbsorption(frequency= 70000, # adjust based on f
                                        salinity = wbat.ctd$salinity,
                                        temperature = wbat.ctd$temperature,
                                        pressure = wbat.ctd$depth,
                                        pH =8,
                                        formulation = "francois-garrison")

wbat.ctd$Sv_corr <- wbat.ctd$abs_coeff - wbat.ctd$Corr_factor 

# create profile by custom bin width for depth

d <- wbat.ctd$depth   #
N <- wbat.ctd$Sv_corr # calibrated Sv_values
db <- seq(0, max(d, na.rm=T), 0.5) #creates a depth range vector at x resolution (start range, end range, bin size[x])
wbat.binned <- binMean1D(d, N, db) # supplies the mean value for each depth bin
wbat.binned<- data.frame(c(wbat.binned["xmids"],wbat.binned["result"]))

d_scale <- seq(0, max(wbat.ctd$depth, na.rm=T), 200)

ggplot(wbat.binned, aes(x=xmids, y=result)) + geom_line(size=1) +
  coord_flip() + scale_x_reverse(breaks=d_scale) +
  labs(x="Depth [m])", y="Sv mean [db re 1m] ") +       # 
#  scale_y_continuous(breaks=c(-90, -85,-80,-75,-70,-65,-60)) +    #200 kHz
  scale_y_continuous(breaks=c(-100,-95,-90,-85,-80)) +    # 70 kHz
  
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






