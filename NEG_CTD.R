library(oce)
library(lattice)
library(plyr)
library(dplyr)
library(data.table)
library(ggalt)
library(ggplot2)
library(akima)
library(rgdal)
library(RColorBrewer)
library(marmap)
library(gplots)
library(raster)
library(rgeos)
library(rgdal)
library(openxlsx)
library(scatterpie)

#read bathymetry (marmap)
bath <- readGEBCO.bathy("C:/Users/jchawarski/OneDrive - Memorial University of Newfoundland/GIS/RN-2586_1525974242773/GEBCO_2014_2D_-21.3471_67.0995_-4.1141_81.7112.nc")
bath <- readGEBCO.bathy("C:/Users/jchawarski/OneDrive - Memorial University of Newfoundland/GIS/RN-3031_1541163235004/GEBCO_2014_2D_-22.6048_62.5151_5.3436_83.3503.nc")

#convert to df for ggplotting
bath.f <- fortify.bathy(bath)

#read coast (rgdal)
coast <- readOGR("C:/Users/jchawarski/OneDrive - Memorial University of Newfoundland/GIS/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp")
#define extent (raster)
coast.trim <- crop(coast, extent(-22.6048, 5.3436, 62.5151, 83.3503))
coast.trim <- crop(coast, extent(-25, 10, 62.5151, 83.3503))


####   SUMMARIZE SEABIRD CTD DATA   ####
{
files <- list.files(path="D:/Greenland Sea/NEG 2017/All CTD Files", full.names = TRUE, pattern= "*.cnv")  # load files from CTD folder
    CTD.cnv <- lapply(files, function(i){ctdTrim(read.ctd.sbe(i))})                           # concatenates and trims upcast from all CTD files into large list
        meta.tbl <- setNames(data.frame(matrix(ncol = 3, nrow = 157)), c("Site", "Date", "ID")) # create empty df for summary data
          meta.tbl[,1] <-  paste(substr(sapply(CTD.cnv, '[[', "filename"), 51, 72))             # selects and trims file name, puttin in blank matrix
          meta.tbl[,2] <- sapply(CTD.cnv, function(i){paste(unique(i[["date"]]))})                # pulls date from each @data in data.in
            uniqueID <- c(1:157)                                                                   # create unique ID by cast
              meta.tbl[,3] <-  uniqueID                                                             # assigns a unique ID

CTD.all <- lapply(CTD.cnv, function(i){ rbind(unique(data.frame(i[["data"]]))) }) %>%  # creates a list of dataframes (of each cast)
               mapply(cbind, ., ID = uniqueID, SIMPLIFY = F) %>%                         # adds ID column with uniqueID to each df in list
                    bind_rows() %>%                                                        # binds all dfs
                       left_join(meta.tbl, ., by="ID")                                      # joins 

CTD.sum.tbl <-  CTD.all %>% filter(depth > 25 & depth <200) %>%                                  # filters all data below 100 m depth
                   group_by(Site, Date) %>%                                             # group fnx, allows headers to carry over
                      dplyr::summarise(                                                 #need to specify dplyr:: fxn, because plyr has a summarise fxn as well
                        avgT = mean(temperature), 
                        avgS = mean(salinity), 
                        lat = mean(latitude, na.rm = T),
                        lon = mean(longitude, na.rm = T),
                        depth.max = max(depth)) 

CTD.sum.tbl <- subset(CTD.sum.tbl, depth.max >190)
CTD.sum.tbl <- subset(CTD.sum.tbl, lat >65)                                  #filter spatial outliers
write.csv(CTD.sum.tbl, "NEG.CTD.100m.TS.Summary.csv")                               #create file


}



#acoustic samples

#SURFACE
files <- list.files(path="F:/NEG/Surface", full.names = TRUE, pattern= "*.csv")
tables <- lapply(files, read.csv, header = TRUE)
surf.df <- do.call(plyr::rbind.fill , tables)
      
    #trim for outliers
    surf.df <- subset(surf.df, NASC < 50000)

#EPIPELAGIC
epi <- read.csv("D:/Greenland Sea/NEG 2017/Leg 1/Leg1D_38kHz_5_25m.csv")
files <- list.files(path="F:/NEG/Epipelagic", full.names = TRUE, pattern= "*.csv")
tables <- lapply(files, read.csv, header = TRUE)
epi.df <- do.call(plyr::rbind.fill , tables)

    #trim for extreme outliers
    epi.df <- subset(epi.df, NASC < 2000)
    epi.df <- subset(epi.df, Depth_mean < 200)
    
#MESOPELAGIC
files <- list.files(path="F:/NEG/Mesopelagic", full.names = TRUE, pattern= "*.csv")
tables <- lapply(files, read.csv, header = TRUE)
meso.df <- do.call(rbind.fill , tables)

    #trim for mesopelagic depth
    meso.df <- subset(meso.df, Depth_mean == 400)

#DEMERSAL
files <- list.files(path="F:/NEG/Demersal", full.names = TRUE, pattern= "*.csv")
tables <- lapply(files, read.csv, header = TRUE)
dem.df <- do.call(rbind.fill , tables)
    
    dem.df <- subset(dem.df, Depth_mean < 800)

#READ CATCH DATA
larv.df <- read.xlsx("D:/Greenland Sea/NEG 2017/Fish larvae NEG 2017.xlsx", sheet = 9)
names(larv.df)[4] <- "lat"
names(larv.df)[5] <- "lon"
larv.df[is.na(larv.df)] <- 0
ggplot() +
geom_scatterpie(data = larv.df, 
                aes(lon, lat, r = sqrt(Total)/15),
                cols = c("n.Bosa", "n.Gymno", "n.Liparis", "n.Triglops"), 
                alpha = 0.5) +
  scale_fill_manual(
    breaks = c("n.Bosa", "n.Gymno", "n.Liparis", "n.Triglops"),
    labels = c("Bosa", "Gymno", "Liparis", "Triglops"),
    values = c("n.Bosa" = "orange",
               "n.Gymno" = "palegreen3",
               "n.Liparis" = "lightslateblue",
               "n.Triglops" = "yellow"))

#akima spline interpolation - create layer

    #acoustic
interpdf <- interp2xyz(interp(x=surf.df$Lon_M, 
                              y=surf.df$Lat_M,
                              z=log(1+(surf.df$NASC)/10.13)), 
                              data.frame = T)
    #CTD
interpdf <- interp2xyz(interp(x=CTD.sum.tbl$lon, 
                              y=CTD.sum.tbl$lat,
                              z=CTD.sum.tbl$avgS), 
                              data.frame = T)


brk <- c(-50,-100,-250,-500,-1000,-1500,-2000)   #define bathymetry breaks
leg.lab <- expression("Mean T" ( degree*C))     # expression for degree
lab <- expression(paste("Log + 1/",mu, "NASC"))


interpdf %>%
  tbl_df() %>%
    ggplot(aes(x = x, y = y, z = z, fill = z)) +               #interpolation
       geom_tile(aes(fill=z)) + 
           scale_fill_distiller(name= "Log NASC", #change 
               type="seq", #change type (div, seq) 
                  palette="Spectral", #change palette (RdBu, Spectral)
                       direction =-1, 
                          na.value = "transparent") +
  
        #      ggplot() +
         #        geom_point(data= dem.df,                         # acoustic points
          #           aes(x=Lon_M,
           #              y=Lat_M,
            #             colour = log(1+NASC/10.13)),  #log(1+NASC/171)
             #               size=2,
              #              inherit.aes = F) +
               #                 scale_colour_distiller(palette = "Spectral") +
                             
             
       #       ggplot() +                     geom_point(data= CTD.sum.tbl,                         # CTD points
        #                                     aes(x=lon, 
         #                                        y=lat,
          #                                      colour = avgS),
           #                                  size=2,
            #                          inherit.aes = F) +
             #                  scale_colour_distiller(palette = "Spectral") +  #adjust palette
                                                     # direction = 1) + 
  
      


  
                geom_contour(data = bath.f,                       #bathymetry
                  aes(x=x, y=y, z=z),
                    breaks=brk,
                      size=c(0.3),
                       colour="lightsteelblue2") +
  
               geom_polygon(data= coast.trim,                      #coast
                            aes(x=long, 
                                y=lat,
                                group=group),
                                  fill= "navajowhite1", 
                                    colour ="dodgerblue4",
                                          inherit.aes = F) +
  
  
  
  
                theme(panel.grid.major =element_line(colour = "snow3",size=0.5),         #theme
                      panel.grid.minor = element_line(colour = "snow3",size=0.5), 
                      plot.title = element_text(hjust=0.5),
                      text = element_text(size=15), 
                      legend.direction = "horizontal",
                      legend.position = "bottom",
                      legend.background = element_blank(),
                      legend.text = element_text(size=10),
                      legend.title = element_text(face="bold"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "white"),
                      axis.text.y = element_text(size=10), 
                      axis.text.x = element_text(size=10),
                      axis.title.x = element_text(size=10),
                      axis.title.y = element_text(size=10),
                      panel.ontop = T) + 
                      theme_linedraw() +
                      labs(title="Interpolated Surface Acoustic Backscatter (5-25 m)", 
                           x = "Longitude",
                           y = "Latitude",
                           colour = "log NASC")  +

                      coord_map("gilbert", xlim=c(-22, 4), ylim=c(73, 80.2)) +
                    #  coord_quickmap(xlim=c(-22, 4), ylim=c(73, 80.2)) +
  
                      p +      geom_scatterpie(data = larv.df, 
                                            aes(Lon, lat, r = sqrt(Total)/10),
                                            cols = c("n.Bosa", "n.Gymno", "n.Liparis", "n.Triglops"), 
                                            alpha = 0.5) +
                            scale_fill_manual(
                              breaks = c("n.Bosa", "n.Gymno", "n.Liparis", "n.Triglops"),
                              labels = c("Bosa", "Gymno", "Liparis", "Triglops"),
                              values = c("n.Bosa" = "orange",
                                         "n.Gymno" = "palegreen3",
                                         "n.Liparis" = "lightslateblue",
                                         "n.Triglops" = "yellow"))  + coord_fixed(1.4) +
                                  ylim(72.5,82.5) + xlim(-22.5,-2.5)

 # coord_map("lambert", lat0=30, lat1=65, xlim=c(-22, 4), ylim=c(73.5, 80)) # need to work on slicing the borders to match projection
  

          
             
  

#extra plots ----  bathymetry
plot(bath, n = 1, image = TRUE, rich.colors(100), main =
       "Greenland Sea Coastal Shelf Bathymetry")

#this map is nice --- need to add sampling points etc.
p <-   autoplot(bath, geom=c("r", "c")) +
scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen") + coord_fixed()
           


#single site
{
st1 <- read.ctd.sbe("26DA.2017.10.101.cnv") #creates object CTD
st1.df <- data.frame(st.1@data) 

st1.df$site <- paste(substr(st.1@metadata[["filename"]], 22, 37)) # selects the filename from metadata, trims, and inserts into Df
st1.df$date <- paste(st.1@metadata[["date"]]) # selects the date from metadata and inserts in Df

st1.sum <- st1.df %>% filter(depth <= 100) %>%
  group_by(site, date) %>%
  summarise(avgT = mean(temperature), 
            avgS = mean(salinity), 
            lat = mean(latitude, na.rm = T),
            lon = mean(longitude, na.rm = T)) 
}

           
           