# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---\
# Case study for the COVID-19 Anthropause project
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Store a pdf for every singke individual animal 


library(DBI)
library(RSQLite)
library(lubridate)
library(raster)
require(spData)
require(tidyverse)
require(sp)
require(sf)
require(stars)
require(mapview)
require(patchwork)
require(ggmap)
require(patchwork)
library(ggsn)
data("us_states")

.dbPF <- "processed_data/mosey_mod.db" 
  # .dbPF = '/Users/diegoellis/projects/Anthropause/analysis/data/mosey_mod_anno.db'
  
  
  
load('raw_data/gHM/ghm_usa.Rdata') # Load GHM
# ghm <- raster("raw_data/gHM/gHM.tif")
ghm_sf = st_as_stars(ghm)

#---- Initialize database ----#
message("Initializing database connection...")
# invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
# invisible(assert_that(length(dbListTables(db))>0))

# Load the entire event table:
evt0 <- tbl(db, "event_clean") %>%
  collect()

dbDisconnect(db)
# indtb <- tbl(db,'individual') %>% collect() # Load a tibble with all individual animals

taxonomic_name <- 'Odocoileus hemionus'
outdir_path <- 'figs/esa_case_study'
path_to_niche <- 'out/niche_determinant_anthropause.csv'
path_to_dbbm <- 'out/dbbmm_size.csv'
# mosey_db # Is the moisy database

make_example_species_figure <- function(mosey_db, ghm, taxonomic_name, outdir_path,
                                        path_to_niche, path_to_dbbm){

  
  pal_col = c('2019'= 'steelblue4', '2020' = 'firebrick')
  
  # Load niches
  niche_det = read.csv(path_to_niche) %>%
    drop_na(total) %>%
    filter(scientificname == taxonomic_name) %>%
    mutate(fyear = factor(year))  %>%
    arrange(fyear, week)  %>% 
    group_by(fyear, week)
  
  
  dbbm_file = read.csv(path_to_dbbm)  %>% 
    filter(species == taxonomic_name) %>% 
    mutate(fyear = factor(year)) %>%
    arrange(fyear, wk)   %>%
    group_by(fyear, wk)  %>%
    select(area, sg, year, wk, ind_id) %>% 
    filter(ind_id %in% niche_det$individual)
  
  # Subset by individuals that have data for 2019 and 2020
  ind_tracking_periods = plyr::ddply(dbbm_file, 'ind_id', function(x){
    range(x$year)
  }) %>% filter(V1 == 2019 & V2 == 2020) %>% dplyr::select(ind_id)
  
  niche_det = niche_det %>% dplyr::filter(individual %in% ind_tracking_periods$ind_id)
  dbbm_file = dbbm_file %>% dplyr::filter(ind_id %in% ind_tracking_periods$ind_id)
  
  # Subset by the moose we are interested in:
  message('Subsetting only to individuals that have data for both 2019-2020')
  evt0_sub = evt0 %>% dplyr::filter(individual_id %in% ind_tracking_periods$ind_id)
  
  message(paste0('Total of ', length(unique(evt0_sub$individual_id))), ' individual animals across ', length(unique(evt0_sub$study_id)), ' studies')
  
  # message('Thinning movement data to 1 point per day')
  # 
  # evt0_sub = data.frame(evt0_sub) %>% 
  #   dplyr::select(individual_id, lon, lat, timestamp) %>%
  #   dplyr::mutate( day = day(timestamp),
  #                  year = year(timestamp)) %>% 
  #   group_by(individual_id, year, day) %>% 
  #   slice(1)
  # 
  
  # Loop through all individuals
  # 
  ind <- unique(evt0_sub$individual_id)[500]
  for(ind in unique(evt0_sub$individual_id)){
    print(ind)
    
    # Plot for a single animal? 
    tmp_ind = evt0_sub %>% dplyr::filter(individual_id == ind)
    
    tmp_ind_sp = SpatialPointsDataFrame(
      coords= tmp_ind[,c('lon', 'lat')], data = tmp_ind,
      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      # proj4string = crs(ghm)
      )
    # tmp_ind_sp_tr <- spTransform(tmp_ind_sp, crs(ghm))
    
    cbbox <- make_bbox(lon = tmp_ind_sp$lon, lat = tmp_ind_sp$lat, f = .4) #from ggmap
    sq_map <- get_map(location = cbbox, maptype = "satellite", source = "google")
    
    
    ghm_sub = raster::crop( ghm,( extent(cbbox[1], cbbox[3], cbbox[2], cbbox[4]) *2) )
    # ghm_sub = raster::crop(ghm, extent(cbbox[1], cbbox[3], cbbox[2], cbbox[4]))
    rtp <- rasterToPolygons(ghm_sub)
    # bm <- bm + geom_raster(...) # insert your raster here
    # 
    # bm <- bm + coord_cartesian()
    # 
    # 
    # geom_tile(data = mean_price_per_bedroom_per_tile, aes(x = Longitude, y = Latitude, alpha = mean_price), fill="blue") + 
    #   scale_alpha(range = c(0, 0.8))
    # 
    tmp_ind$year = factor(tmp_ind$yr) 
    
    
    case_map = ggmap(sq_map) + 
      geom_polygon(data = rtp, 
                   aes(x = long, y = lat, group = group,
                       fill = rep(rtp$gHM, each = 5)), 
                   size = 0, 
                   alpha = 0.5)  +
      scale_fill_viridis_c(option = 'plasma', na.value="NA") + 
      geom_point(data = tmp_ind, aes(x = lon, y = lat, color = year)) +
      geom_line(data = tmp_ind, aes(x = lon, y = lat, color = year)) +
      labs(x = " ", y = " ", title = paste0(taxonomic_name, " -Individual ", ind)) +
      theme_minimal() +
      scale_colour_manual(values = pal_col)+
      guides(fill=guide_legend(title="Human modification")) +
      #       coord_cartesian() +
      #    geom_stars(data = st_as_stars(ghm_sub), color = alpha(0.8))  +
      NULL
    
    niche_plot =
      (
        niche_det %>% 
          dplyr::filter(individual == ind) %>% 
          ggplot(aes(x = week, y = total, group = year), fill = year) + 
          # geom_line(aes(color = holc_grade), size = 1) +
          geom_smooth(aes(color = fyear)) + 
          labs(y = 'Total niche breadth'
               , x = '') +
          ggtitle('Niche breadth') +
          theme_minimal() + 
          # scale_colour_viridis_b(option = 'plasma')
          scale_colour_manual(values = pal_col) +
          theme(legend.position = "none") +
          
          #  facet_grid(~individual) +
          # scale_colour_brewer()
          NULL
      )
    # niche_plot
    
    dBBM_plot = 
      (
        dbbm_file %>% 
          dplyr::filter(ind_id == ind) %>% 
          ggplot(aes(x = wk, y = area, group = fyear), fill = fyear) + 
          # geom_line(aes(color = holc_grade), size = 1) +
          geom_smooth(aes(color = fyear)) + 
          labs(y = 'Space use in km^2'
               , x = '') +
          ggtitle('Space use') +
          theme_minimal() + 
          # scale_colour_viridis_b(option = 'plasma')
          scale_colour_manual(values = pal_col) +
          theme(legend.position = "none") +
          # facet_grid(~ind_id) +
          # scale_colour_brewer()
          NULL
      )
    # dBBM_plot
    
    sg_plot = 
      (
        dbbm_file %>% 
          dplyr::filter(ind_id == ind) %>% 
          ggplot(aes(x = wk, y = sg, group = fyear), fill = fyear) + 
          # geom_line(aes(color = holc_grade), size = 1) +
          geom_smooth(aes(color = fyear)) + 
          labs(y = 'Safegraph mobility count'
               , x = 'Week of the year') +
          ggtitle('Human mobility') +
          theme_minimal() + 
          scale_colour_manual(values = pal_col) +
          theme(legend.position = "none") +
          #  facet_grid(~ind_id) +
          NULL
      )
    
    
    # Inset
    lonlat_to_state <- function(pointsDF,
                                states = spData::us_states,
                                name_col = "NAME") {
      ## Convert points data.frame to an sf POINTS object
      pts <- st_as_sf(pointsDF, coords = 1:2, crs = 4326)
      
      ## Transform spatial data to some planar coordinate system
      ## (e.g. Web Mercator) as required for geometric operations
      states <- st_transform(states, crs = 3857)
      pts <- st_transform(pts, crs = 3857)
      
      ## Find names of state (if any) intersected by each point
      state_names <- states[[name_col]]
      ii <- as.integer(st_intersects(pts, states))
      state_names[ii]
    }
    
    
    testPoints = data.frame(
      x = c(bbox(tmp_ind_sp)[1,1], bbox(tmp_ind_sp)[1,2]),
      y = c(bbox(tmp_ind_sp)[2,1], bbox(tmp_ind_sp)[2,2])
    )
    
    unique(lonlat_to_state(testPoints))
    
    us_tmp = us_states %>% dplyr::filter(NAME %in% unique(lonlat_to_state(testPoints))    )
    
    
    # Make USA inset
    inset = ggplot() + 
      geom_sf(data = us_states, fill = "grey49", color='transparent') + 
      geom_sf(data = us_tmp, fill = NA, color = "red", size = 1.2) +
      theme_void()
    inset
    
    
    
    caseplot = (case_map + inset_element(inset, left = 1, right = 1.3, top = 1, bottom = 0.8)) + niche_plot / dBBM_plot / sg_plot
    
    # caseplot = case_map + niche_plot / dBBM_plot / sg_plot
    
    ggsave(paste0(outdir_path, taxonomic_name, ' ',ind, '.png'), caseplot)
    
    
  } # Finish looping through individuals
  
  
  
  
} # End of function
