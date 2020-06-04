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

### Sv - MEAN VOLUME BACKSCATTER - ###

setwd("C:\\Users\\jchawars\\OneDrive - Memorial University of Newfoundland\\NEG")

CTD <- read.ctd.sbe("CTD\\UpDown\\26da.2017.9.395_110618.cnv")  #read a single CTD files that coincides with WBAT profile
CTD.data <- data.frame(CTD[["data"]])  #select just the data portion

time <- CTD[["metadata"]]$startTime   # collect the start time of the cast
time_seconds <- period_to_seconds(lubridate::seconds(time))  # convert start time to seconds (timeY)

CTD.data$Interval <- CTD.data$timeS + time_seconds - 5 # convert timeS to seconds by adding to start time
                             # number of seconds from adjustments in CTD_WBAT-Time-at-depth.csv
# Read in the WBAT data
wbat <- read.csv2("WBAT\\All_Sv_70_200kHz\\St.80_70kHz_Sv_T20170908_21582701-20170908_22142300.txt", sep = ",", skip=6,header=FALSE, fill = TRUE)
#wbat <- read.csv2("St.3_70kHz_Sv.txt", sep = ",", skip=6,header=FALSE, fill = TRUE)

# Give new names to columns:
nms <- c("PingNumber","Frequency","Date","Time","Latitude","Longitude","RangeStart","RangeStop","DepthStart","DepthStop","SampleCount", sprintf("Sv%02d", seq(1,(dim(wbat)[2])-11)))

# Assign columns to data.frame
colnames(wbat) <- nms
wbat <- wbat[!is.na(wbat$SampleCount),] # Removing NAs (for some reason one extra column (after column #3707) with no column name but only Sv values were exported from LSSS. This messed up the data.frame, so I chose to remove these values)

## Absorption coefficient profile correction

wbat.meta <- wbat[1:dim(wbat)[1],1:11]      # select the metadata portion 
head(wbat.meta)
# convert the probe time to the same format as the CTD -- NOTE: ECHOVIEW 11 has fixed their export headers...


wbat.meta$Time <- gsub("(..)(..)(..)(..)", "\\1:\\2:\\3:\\4", wbat.meta$Time) # adjust the first (..) to (.) for single digit hours... smh
wbat.meta$Date<- as.Date(as.character(wbat.meta$Date), "%Y%m%d")
wbat.meta$datetime <- as.POSIXct(paste(wbat.meta$Date, wbat.meta$Time), format="%Y-%m-%d %H:%M:%S", tz="UTC")

wbat.meta$Interval <- period_to_seconds(lubridate::seconds(wbat.meta$datetime))   # add or subtract seconds to adjust the time difference between instruments

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
c <-1440.94             # original sound speed
coeff_abs <- 0.01678    # original coefficient of absorption
t <-1.95 *10^-4         # PulseCompressedEffectivePulseLength (sec)
y <- -13                # two-way beam angle  --- Adjust based on transducer 

f_nominal <-70000 # nominal frequency in Hz   -- Adjust based on transducer
f <- 70000        # center frequency         NOTE: for wideband data, these are different values

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

## Plotting the corrected mean Sv value as a profile ##

Sv_label <- expression(paste("Mean S"["v"]," [dB re 1 m" ^-1,"]"))
d_scale <- seq(0, max(wbat.ctd$pressure), 100) # set depth scale for plotting

# comparison plot -- this plot shows both the original and the adjusted mean Sv value for the profile
 ggplot(wbat.ctd, aes(x=pressure, y=Estimate_new)) + 
  geom_line(colour="blue", size=0.5) +  geom_ribbon(aes(ymin=lowCI_new, ymax=highCI_new), linetype=1, alpha=0.2, fill="blue") + 
  geom_line(aes(x=pressure, y=Estimate_old), inherit.aes = F, colour="red") + 
  geom_ribbon(aes(ymin=lowCI_old, ymax=highCI_old), linetype=1, alpha=0.2, fill="red", inherit.aes = T) +
  coord_flip()+ scale_x_reverse(breaks=d_scale) + labs(x= "Pressure [dBar]", y=Sv_label) +
  # xlim(510, 490) + ylim(-90,-85) +
  theme_bw()

       wbat.ctd$wbat_site <- "Stn.84"
       wbat.ctd$ctd_site <- "26da.2017.9.395_110618"
       write.csv(wbat.ctd, "NEG2017_Stn.80_Gr.395_70kHz_WBAT-CTD-SV.csv")
       
 
## SV plotting ## 
 
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
  xlim(max(wbat.ctd$pressure),0) + ylim(5,190) +
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
comp.plot <- plot_grid(sv.plot, T.plot, S.plot, O.plot,  nrow=1, ncol=5,rel_widths = c(2, 1, 1, 1))

title <- ggdraw() + 
  draw_label(
    "NEG 2017 Stn 3 Gear 8 - 200 kHz WBAT/CTD",
    fontface = 'bold',
    x = 0,
    hjust = -0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
plot_grid(
  title, comp.plot,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)


# .csv filename NEG2017_Stn.3_Gr.8_200kHz_WBAT-CTD
# .png filename NEG_St.3_Gr.8_200kHz_wbatctd

files <- list.files(path="C:\\Users\\jchawars\\OneDrive - Memorial University of Newfoundland\\NEG\\Processed WBAT\\SV", full.names = T, pattern= "*.csv")  # load files from CTD folder

probe.all <- lapply(files, function(i) read.csv(i)) %>% 
  lapply(., mutate_if, is.integer, as.character) %>% 
  lapply(., mutate_if, is.numeric, as.character) %>% 
  bind_rows()

write.csv(probe.all, "NEG2017_WBAT_AllSv.csv")

%>%
  dplyr::select(wbat_site,
                Frequency,
                temperature, 
                salinity, 
                pressure,
                Estimate_new, 
                lowCI_new,
                highCI_new) %>%
  # Latitude,
  # Longitude
  
  mutate_at(1:2, factor) %>% 
  mutate_at(3:7, as.numeric)

Sv_label <- expression(paste("Mean S"["v"]," [dB re 1 m" ^-1,"]"))

ggplot(probe.all, aes(x=pressure, y=Estimate_new)) + geom_point(aes(colour=Frequency)) +
  stat_smooth(formula = y ~ s(x, k = 45), method = "gam", se = T, linetype="dashed", color="black", aes(colour=Frequency)) + 
  coord_flip()+ scale_x_reverse() + labs(x= "Pressure [dBar]", y=Sv_label) +
  theme_bw()










## TS - Target depth alignment and Target Strength Analysis ##

setwd("C:\\Users\\jchawars\\OneDrive - Memorial University of Newfoundland\\NEG")

CTD <- read.ctd.sbe("CTD\\UpDown\\26da.2017.9.395_110618.cnv")  #read a single CTD files that coincides with WBAT profile
CTD.data <- data.frame(CTD[["data"]])  #select just the data portion

time <- CTD[["metadata"]]$startTime   # collect the start time of the cast
time_seconds <- period_to_seconds(lubridate::seconds(time))  # convert start time to seconds (timeY)

CTD.data$Interval <- CTD.data$timeS + time_seconds - 5  # convert timeS to seconds by adding to start time and  add or subtract seconds to adjust the time difference between instruments
# number of seconds from adjustments in CTD_WBAT-Time-at-depth.csv

wbat.ts <- read.csv2("WBAT\\70kHz TS\\St.80_70kHz_minTS-90dB_TS_T20170908_21570201-20170908_22142300.txt", sep = ",", skip=19,header=FALSE, fill = TRUE)
# Give new names to columns:
nms <- c("Date","Time","Latitude","Longitude","Range","TSC","TSU","AlongshipAngle","AthwartshipAngle","sV_of_peak")

# Assign columns to data.frame
colnames(wbat.ts) <- nms
head(wbat.ts)

wbat.ts$Time <- gsub("(..)(..)(..)(..)", "\\1:\\2:\\3:\\4", wbat.ts$Time) # adjust the first (..) to (.) for single digit hours... smh
wbat.ts$Date<- as.Date(as.character(wbat.ts$Date), "%Y%m%d")
wbat.ts$datetime <- as.POSIXct(paste(wbat.ts$Date, wbat.ts$Time), format="%Y-%m-%d %H:%M:%S", tz="UTC")
wbat.ts$Interval <- period_to_seconds(lubridate::seconds(wbat.ts$datetime))   

wbat.ctd <- wbat.ts %>% 
  left_join(., CTD.data, by="Interval") %>%     # join ctd and wbat data.frames by matching time cases
  mutate_all(., as.numeric)          # convert all values to numeric

sal <- wbat.ctd$salinity
temp <- wbat.ctd$temperature
p <- wbat.ctd$pressure

require(oce)
c_new <- swSoundSpeed(salinity = sal,temperature = temp, pressure = p) # new sound speed

# set constant variables from orginial calibration file -- These are all dependent on transducer type and default data colllection settings
c <-1440.94             # original sound speed
coeff_abs <- 0.01678    # original coefficient of absorption
t <-1.95 *10^-4         # PulseCompressedEffectivePulseLength (sec)
y <- -13                # two-way beam angle  --- Adjust based on transducer 

f_nominal <-70000 # nominal frequency in Hz   -- Adjust based on transducer
f <- 70000        # center frequency         NOTE: for wideband data, these are different values

equi_beam_angle <-10^((y+20*log(f_nominal/f))/10)          # calculate equivalent beam angle
coeff_abs_new <-swSoundAbsorption(frequency= f_nominal,    # new absorption coefficient
                                  salinity = sal,
                                  temperature = temp,
                                  pressure = p, 
                                  pH =8,
                                  formulation = "francois-garrison") 

wbat.ctd$signal_time <- (c*t/4+wbat.ctd$Range)*(2/c) #time for the signal to return 


# calculate new range based on sound speed
wbat.ctd$Range_new <- wbat.ctd$signal_time*c_new/2-c_new*t/4   #new range matrix

# calculate power matrix
wbat.ctd$power<- wbat.ctd$TSC-40*log10(wbat.ctd$Range)-2*(coeff_abs)*(wbat.ctd$Range) 

## Calculate exact TS for each target using the equation TS = Power + 40logR + 2* alpha * R_new 
wbat.ctd$TSC_new <- wbat.ctd$power+40*log10(wbat.ctd$Range_new)+2*coeff_abs_new*wbat.ctd$Range_new  

ggplot(wbat.ctd, aes(x=pressure, y=TSC_new)) + 
  geom_point(colour="blue", size=0.5) +  
  coord_flip()+ scale_x_reverse() + labs(x= "Pressure [dBar]", y= "Target Strength [dB]") +
  theme_bw()


wbat.ctd$wbat_site <- "Stn.80"
wbat.ctd$threshold <- "-90dB"
write.csv(wbat.ctd, "NEG2017_Stn80._Gr.395_70kHZ_-90db_WBAT-CTD-TS.csv")

wbat.ctd %>% 
  mutate(bin_dist = factor(pressure%/%bin_size*50)) %>% 
  ggplot(aes(x = bin_dist, y = TSC_new)) +  geom_boxplot() + coord_flip() + theme_bw()





