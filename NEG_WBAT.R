# WBAT/CTD profile synchronization and plotting

## Developed by Julek Chawarski [ julek.chawarski@mi.mun.ca ]
## Initial commit: April-05-2020

### This script was developed to facilitate the synchronization of WBAT profiles with depth 
### and other environmental parameters collected by the CTD-rosette. WBAT profiles are 
### processed and cleaned using impulse noise removal in Echoview 11. WBAT profiles are exported
### as .csv by individual pings (Cell export 1 ping X 100 m bins).

# # # NEG2017 CTD data only records timeS (time in seconds) in the time slot, therfore we use startime as a reference
# # # Process varies slightly for CTDs with TimeY and TimeJ variables

library(ggplot2)
library(dplyr)
library(oce)
library(reshape2)
library(lubridate)
library(imputeTS)

setwd("C:\\Users\\jchawars\\OneDrive - Memorial University of Newfoundland\\NEG")

CTD <- read.ctd.sbe("CTD\\UpDown\\26da.2017.9.8_110618.cnv")  #read a single CTD files that coincides with WBAT profile
CTD.data <- data.frame(CTD[["data"]])  #select just the data portion

time <- CTD[["metadata"]]$startTime   # collect the start time of the cast
time_seconds <- period_to_seconds(lubridate::seconds(time))  # convert start time to seconds (timeY)

CTD.data$Interval <- CTD.data$timeS + time_seconds # convert timeS to seconds by adding to start time

# Read in the WBAT data
wbat <- read.csv2("WBAT\\Export\\St.3_200kHz_Sv.txt", sep = ",", skip=6,header=FALSE, fill = TRUE)
#wbat <- read.csv2("St.3_70kHz_Sv.txt", sep = ",", skip=6,header=FALSE, fill = TRUE)

# Give new names to columns:
nms <- c("PingNumber","Frequency","Date","Time","Latitude","Longitude","RangeStart","RangeStop","DepthStart","DepthStop","SampleCount", sprintf("Sv%02d", seq(1,3696)))

# Assign columns to data.frame
colnames(wbat) <- nms
wbat <- wbat[!is.na(wbat$SampleCount),] # Removing NAs (for some reason one extra column (after column #3707) with no column name but only Sv values were exported from LSSS. This messed up the data.frame, so I chose to remove these values)

## Absorption coefficient profile correction

wbat.meta <- wbat[1:dim(wbat)[1],1:11]      # select the metadata portion 

# convert the probe time to the same format as the CTD -- NOTE: ECHOVIEW 11 has fixed their export headers...

wbat.meta$Time <- gsub("(..)(..)(..)(..)", "\\1:\\2:\\3:\\4", wbat.meta$Time)
wbat.meta$Date<- as.Date(as.character(wbat.meta$Date), "%Y%m%d")
wbat.meta$datetime <- as.POSIXct(paste(wbat.meta$Date, wbat.meta$Time), format="%Y-%m-%d %H:%M:%S", tz="UTC")

wbat.meta$Interval <- period_to_seconds(lubridate::seconds(wbat.meta$datetime))

wbat.ctd <- wbat.meta %>% 
  left_join(., CTD.data, by="Interval") %>%     # join ctd and wbat data.frames by matching time cases
  mutate_all(., as.numeric) %>%           # convert all values to numeric
  group_by(PingNumber) %>%                # there are multiple values per ping so they need to be average to wbat sampling resolution
  summarise_all(., mean) %>%               # calculate mean value for each parameter per ping
  na_interpolation(.)

sal <- wbat.ctd$salinity
temp <- wbat.ctd$temperature
p <- wbat.ctd$pressure

require(oce)
c_new <- swSoundSpeed(salinity = sal,temperature = temp, pressure = p) # new sound speed

# set constant variables from orginial calibration file -- These are all dependent on transducer type and default data colllection settings
c <-1482.41             # original sound speed
coeff_abs <- 0.05    # original coefficient of absorption
t <-2.14 *10^-4         # PulseCompressedEffectivePulseLength (sec)
y <- -20.7                # two-way beam angle  --- Adjust based on transducer 

f_nominal <-200000 # nominal frequency in Hz   -- Adjust based on transducer
f <- 200000        # center frequency         NOTE: for wideband data, these are different values

equi_beam_angle <-10^((y+20*log(f_nominal/f))/10)          # calculate equivalent beam angle
coeff_abs_new <-swSoundAbsorption(frequency= f_nominal,    # new absorption coefficient
                                  salinity = sal,
                                  temperature = temp,
                                  pressure = p, 
                                  pH =8,
                                  formulation = "francois-garrison") 

wbat.sv <- wbat[1:dim(wbat)[1], 12:dim(wbat)[2]]  # selects only the Sv matrix from the exported file
wbat.sv[wbat.sv==""] <- -999                    # replaces the blank cells of matrix
wbat.sv <- matrix(as.numeric(unlist(wbat.sv)),nrow=nrow(wbat.sv))

# Set variables for range matrix calculation
Range_start <- as.numeric(wbat.meta$RangeStart)[1] # all values in the df are the same.
Range_stop <- as.numeric(wbat.meta$RangeStop)[1]
Sample_count <- as.numeric(wbat.meta$SampleCount)[1]

# calculate range for each horizontal sample
wbat.range <- sapply(1:dim(wbat.sv)[2], function(i) Range_start + ((Range_stop-Range_start)/Sample_count)*(i+0.5))    

wbat.range <- matrix(wbat.range, nrow=dim(wbat.sv)[1], ncol = dim(wbat.sv)[2], byrow = T) # creates the full range matrix

# calculate time matrix
time.mat <- apply(wbat.range,1:2, function(i)(c*t/4+i)*2/c) #time for the signal to return 

# calculate a new range matrix
wbat.range_new <- apply(time.mat, 2, function(i) i*c_new/2-c_new*t/4)   #new range matrix

# calculate power matrix
power<- wbat.sv-20*log10(wbat.range)-2*(coeff_abs)*(wbat.range)+10*log10(c*t*equi_beam_angle/2) 

## Calculate exact Sv for each sample using the equation Sv = Power + 20logR * alpha * R_new - 10log10(c*tau*psi/2)
Sv_new <- power+20*log10(wbat.range_new)+2*coeff_abs_new*wbat.range_new  # -10*log10(c_new*t*equi_beam_angle/2) 
x <- 10*log10(c_new*t*equi_beam_angle/2)  # 10log10(c*tau*psi/2) is a vector - subtraction of vector from matrices requires special formulas  
Sv_new <- apply(Sv_new, 2, function(i) i-x) # subtract the final term from the rest of the calculation

## Complete the dataframes 

#replace -999 values with NA
wbat.sv[wbat.sv < -990] <- NA
Sv_new[Sv_new < -990] <- NA

require(gmodels)
# calculate confidence intervals
Sv_new_ci <- t(as.data.frame(apply(Sv_new, 1, function(i) gmodels::ci(i, na.rm=T))))
colnames(Sv_new_ci) <- c("Estimate_new", "lowCI_new", "highCI_new", "stderr_new")

Sv_old_ci <- t(data.frame(apply(wbat.sv, 1, function(i) gmodels::ci(i, na.rm=T))))
colnames(Sv_old_ci) <- c("Estimate_old", "lowCI_old", "highCI_old", "stderr_old")

# bind the mean and 95% condfidence intervals with the CTD data
wbat.ctd <- cbind(wbat.ctd, Sv_old_ci, Sv_new_ci)
wbat.ctd$wbat_site <- "Stn.3"
wbat.ctd$ctd_site <- "26da.2017.9.8_110618"


write.csv(wbat.ctd, "NEG2017_200khz_wbat_ctd.csv")
wbat.ctd <- read.csv("NEG2017_200khz_wbat_ctd.csv")

## Plotting the corrected mean Sv value as a profile ##

Sv_label <- expression(paste("Mean S"["v"]," [dB re 1 m" ^-1,"]"))
d_scale <- seq(0, max(wbat.ctd$pressure), 200) # set depth scale for plotting

# comparison plot -- this plot shows both the original and the adjusted mean Sv value for the profile
sv.plot <-  ggplot(wbat.ctd, aes(x=pressure, y=Estimate_new)) + 
  geom_line(colour="blue", size=0.5) +  geom_ribbon(aes(ymin=lowCI_new, ymax=highCI_new), linetype=1, alpha=0.2, fill="blue") + 
  geom_line(aes(x=pressure, y=Estimate_old), inherit.aes = F, colour="red") + 
  geom_ribbon(aes(ymin=lowCI_old, ymax=highCI_old), linetype=1, alpha=0.2, fill="red", inherit.aes = T) +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y=Sv_label) +
  # xlim(510, 490) + ylim(-90,-85) +
  theme_bw()

# Corrected Sv profile
sv.plot <-  ggplot(wbat.ctd, aes(x=pressure, y=Estimate_new)) + 
  geom_line(colour="black", size=0.5) +  geom_ribbon(aes(ymin=lowCI_new, ymax=highCI_new), linetype=1, alpha=0.2, fill="blue") + 
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y=Sv_label) +
  theme_bw()


## Plotting a depth-referenced echogram from acoustic probe ##
require(imager)
echogram <- as.cimg(Sv_new)  # converts Sv matrix into an image file
map.shift <- function(x,y) list(x=wbat.ctd$pressure,y=y)  
echo.shift <- imwarp(echogram,map=map.shift, direction="forward", coordinates = "absolute", boundary="dirichlet", interpolation = "nearest") 
echo.shift[echo.shift > -10] <- NA
echo.small <- as.data.frame(resize(echo.shift,round(width(echo.shift)/5),round(height(echo.shift)/5),  # reduces the Sv matrix to roughly 4% size for easy plotting
                                   interpolation_type = 2)) # interpolates using a moving average interpolation

require(pals)
echoplot <- ggplot(echo.small, aes(x*5, y/4)) +  
  geom_raster(aes(fill=value)) +
  scale_fill_gradientn(colours=ocean.haline(100), guide = "colourbar", name="MVBS")  + 
  coord_flip() + scale_x_reverse(breaks=d_scale) + 
  xlim(980,0) + ylim(5,190) +
  labs(x="Pressure [dBar]", y= "Range [m]") +
  theme(legend.position="left") +
  theme_bw()

## Building custom CTD plots

T.plot <- ggplot(wbat.ctd, aes(x=pressure, y=temperature)) + geom_line(colour="blue") +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y="Temperature [C]") +
  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank())

S.plot <- ggplot(wbat.ctd, aes(x=pressure, y=salinity)) + geom_line(colour="red") +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y="Salinity [psu]") +
  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank())

F.plot <- ggplot(wbat.ctd, aes(x=pressure, y=fluorometer)) + geom_line(colour="green") +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y="Flourescence [mg/m**3]") +
  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank())

SP.plot <- ggplot(wbat.ctd, aes(x=pressure, y=par)) + geom_line(colour="orange") +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y="PAR [µeinsteins/s/m^2]") +
  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank())

O.plot <- ggplot(wbat.ctd, aes(x=pressure, y=oxygen)) + geom_line(colour="purple") +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y="Oxygen [µM]") +
  theme_bw() + theme(axis.title.y=element_blank(), axis.text.y=element_blank())

Ham.plot <- ggplot(wbat.ctd, aes(x=pressure, y=Ham)) + geom_line(colour="purple") +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y="HAM") +
  theme_bw()# + theme(axis.title.y=element_blank(), axis.text.y=element_blank())


require(cowplot)
plot_grid(echoplot, sv.plot, T.plot, S.plot, O.plot,  nrow=1, ncol=5,rel_widths = c(2, 1, 1, 1, 1))








