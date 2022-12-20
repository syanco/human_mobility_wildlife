# --- --- --- --- --- --- --- --- --- ---
# To do: 
# --- --- --- --- --- --- --- --- --- ---

# Add space use dbmm to visualize: # Add DBMM 95% or 99% 
# Password

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Shiny App for 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

conflict_prefer("box", "shinydashboard")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("summarize", "dplyr")
require(plyr)
library(tidyverse)
library(ggthemes)
require(gridExtra)
library(patchwork)
library(brms)
library(grid)
library(cowplot)
require(lubridate)
library(dplyr)
library(lubridate)
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(ggplot2)
library(leaflet)
library(viridis)
require(shinythemes)
require(geohashTools)
library(googledrive)
library(googlesheets4)
library(rsconnect)
require(plyr)
col <- "#D3D3D3"

# Authentification details:
# gs4_user()
# AUTH_KEY_PATH = './keys/mol-wi-shiny-gdrv-map-of-life-bd9849bf59f3.json'
# DRIVE_PATH = 'https://drive.google.com/drive/folders/1iFgIDGZZ_T5Gh7bKC4Rvz3uqg_b8P22y?usp=sharing' # Link to folder
# SHEET_PATH = 'https://docs.google.com/spreadsheets/d/1g5N27qWgE6ZzIdlgwOA0lAxXudbtTObrc8hdhjTTW2Q/edit?usp=sharing' # Link to .csv
# drive_auth(path = AUTH_KEY_PATH)
# drive_token = drive_auth()
# shiny_token <- gs4_auth() #
# shiny_token = gs4_auth()
# gs4_auth(token = drive_token())
# saveRDS(shiny_token, "shiny_app_token.rds")
# Load model summaries
area_mod_summary = read.csv('/Users/diegoellis/projects/Anthropause/out/area_mod_summary_2022-10-05.csv')
area_mod_summary %>% filter(inter_sig == 'Y') # Add interaction code in shiny app
area_mod_summary %>% filter(!inter_sig == 'Y') # Add interaction code in shiny app

niche_mod_summary = read.csv('/Users/diegoellis/Downloads/niche_mod_summary_2022-10-24.csv')

# Load path to models
indir_space_models = '/Users/diegoellis/projects/Anthropause/area_additive/'
indir_niche_models = '/Users/diegoellis/projects/Anthropause/niche_controlled/'
indir_space_interaction_models = '/Users/diegoellis/projects/Anthropause/area_interaction/'
# vector of species 
species = unique(area_mod_summary$species)
# Need to fix niche plot: Maybe if else skip statement render plot of niche figure 
species = unique(niche_mod_summary$species)

# Load one point pe day mosie DB # Load and format data
# locs = read.csv('/Users/diegoellis/projects/Anthropause/one_point_per_day_20221118.csv') %>%
locs = read.csv('/Users/diegoellis/projects/Anthropause/one_point_per_day_20221118_anno.csv') %>%
  mutate(
    # month = month(timestamp),
    # year = year(timestamp),
    # date = ymd(as.Date(timestamp)),
    col_level = as.factor(year)
  )

folder_dbmm <- list.files('/Users/diegoellis/projects/Anthropause/dbbmms/', full.names = T, pattern = 'dbbmm_')

dbbmm_size = read.csv('/Users/diegoellis/Downloads/dbbmm_size.csv') %>% mutate(    col_level = as.factor(year))
niche_det = read.csv('/Users/diegoellis/Downloads/niche_determinant_anthropause.csv') %>% mutate(    col_level = as.factor(year))

ui <- fluidPage(
  theme = shinytheme("cerulean"),
  # to enable boxes from shinyDashboard
  useShinydashboard(),
  
  # Row with selectors ----
  fluidRow(
    column(3,
           # App title ----
           titlePanel("Data owner App"),
           sidebarLayout(
             sidebarPanel(HTML('<b>Sharing a crowded planet: COVID-19 reveals  the footprint of human activity on wildlife<b>'),
                          h6(HTML('Team Leads: Ruth Oliver 1,2, Scott Yanco 1, Diego Ellis Soto 1, Walter Jetz 1 <br/> 1. Yale University, 2. University of California Santa Barbara <br>
                                  Shiny App developed by Diego Ellis Soto')),
                          width = 12),
             mainPanel(
               # 'Please select your study species and studyID. You will be able to display individual tagged animals, and the amount of human infrastructure and movement of humans they were exposed to. Bellow, you will also find species level predictions of space and habitat use '
               img(id ='bls_logo', src = "bls_logo.png", height = 100, width = 300, position = "right"),
               # tags$img(src = 'bls_logo.png'),
               # includeHTML("https://www.bio-logging.net"),
               position = "right")
             , position = "right"
           )
    )
    
    # ,column(10,
    #          # App title ----
    #          titlePanel(h6(HTML('Team Leads: <br/> Ruth Oliver 1,2, Scott Yanco 1, Diego Ellis Soto 1, Walter Jetz 1 <br/> 1. Yale University, 2. University of California Santa Barbara')))# ,
    # )
  ),
  fluidRow(
    
    column(12,
           # App title ----
           titlePanel(h4("Modeling wildlife responses")),
           sidebarLayout(
             
             sidebarPanel(      (HTML(
               '</b>We modelled the response of 40 species of terrestrial birds and mammals to human mobility and infrastructure across the United States.</b> <br> 
                <br> 
            All our analysis are based on individual animal trajectories from 2019-2020 which are displayed bellow for each individual separately. For each individual animal, we estimated weekly geographic space use and niche breadth and captured these responses in the context of </b>human modification (static) and mobility (dynamic)</b>.  <br> 
             <br> 
            A diagram of our methodology is detailed bellow. We aggregated individual level space use and niche breadth to make species level inference.<br/>
             <br> 
               Our species level models allow making inference on the effect of the amount of area used, human modification, and mobility on the breadth of environmental niches. <br/>'))),
             mainPanel(
               img(src = "species_responses.png", height = 500, width = 800, position = "right")
             )
           )
    )
  ),
  fluidRow(
    
    titlePanel(h4("Requested feedback:")),
    sidebarLayout(
      sidebarPanel(      HTML(
        '<b>-We kindly ask you to fill out our Google Survey embedded in the link bellow.<b><br>
      -<b>Please add your contact information, animal handling permits, desired data sharing options, and funding acquisition <b><br>
      -<b>Please visually inspect the tracking data of individual animals</b><br>
      -<b>Please provide feedback on the modeled space use and niche breadth respones of your study species. Do these effects seem plausible?</b><br>
      -<b>We recognize that these data are often aggregated across multiple study sites and may represent species level averages.</b><br>'
      )),
      mainPanel()
    )
  ), # End of fluid row
  
  fluidRow(
    column(2, 
           # Input: Dropbox for the species ----
           selectInput(inputId = "Species", label = strong("Select a species:"), choices = "", selected = "Alces alces", multiple = FALSE)),
    column(2, 
           # Input: Dropbox for the individual ----
           selectInput(inputId = "Individual", label = strong("Select an individual:"), choices = "", selected = "", multiple = FALSE)),
    column(2, 
           # Input: Dropbox for the individual ----
           selectInput(inputId = "StudyID", label = strong("Select a studyID:"), choices = "", selected = "", multiple = FALSE))
  ), # End of fluid row
  
  
  fluidRow(
    column(12,
           # App title ----
           titlePanel(h4("Results for Individual animal")),
    )
  ),
  
  
  
  fluidRow(
    box(
      column(2,
             titlePanel(
               h5("Tracking data for individual animal", align = 'left'))),  
      
      column(2,
             
             title = "Choose data aggregation products",  status = 'primary', solidHeader = TRUE,
             collapsible = FALSE, collapsed = FALSE,
             checkboxGroupInput('Data_for_map', h6('Select aggregation:'), # Data for map gets rendered:
                                c('All data','Geohash 5', 'Geohash 3', 'Space use'),
                                # c('All data','Geohash 6','Geohash 5', 'Geohash 3', 'Space use'),
                                selected = 'All data'),
             #      c('Geohash', 'Space use'))
      ),
      
      column(8,
             # Output: map ----
             leafletOutput("mymap")),
    ),
    # MOVED THIS BELLOW: 
    # column(4,
    #        # Output: plotly safegraph
    #        plotOutput(outputId = "indiv_scatterplot_ghm_sg"))
  ),
  #   column(2,
  #          # Output: plotly safegraph
  #          plotOutput(outputId = "indiv_scatterplot_sg"))
  # ),
  # 
  
  fluidRow(
    
    column(6,
           # Output: plotly safegraph
           plotOutput(outputId = "indiv_scatterplot_ghm_sg"))
  ,
  
    column(6,
           # Output: plotly safegraph
           plotOutput(outputId = "dbbmm_niche_week_size_id"))# ,
    # column(2,
    #        # Output: plotly safegraph
    #        plotOutput(outputId = "niche_week_id"))
  ),
  
  
  column(12,
         # App title ----
         titlePanel(h4("Species level results")),
         sidebarLayout(
           
           # strong("is fun") 
           
           sidebarPanel(      (HTML(
             '</b>-Space use models:</b> We estimated weekly space use for each individual using dynamic Brownian Bridge Movement. We then aggregated these individual space use models to the species level. This allowed making predictions of weekly space use under human modification and human mobility. We calculated these models using both an interactive effect, as well as additive effects of human mobility and modification separately. In this app, interactive models will be displayed if a significant interactive effect is shown. Otherwise, additive space use models are displayed. <br/>
             <br>
             </b>-Niche breadth models:</b> We estimated weekly niche breadth after controlling for individual animal space use.
           For all species, we considered the following variables as components of an individual’s realized niche: NDVI, LST and elevation. <br/>'
             
           ))),
           mainPanel(
           )
         )
  ),
  fluidRow(
    column(12,
           titlePanel(
             # h4(paste0("Model output for study species ", output$species_id ), align = 'center')
             h4(textOutput(outputId = 'titulo'), align = 'center')
           )
           # In the server part.
           
    )
  ),
  
  fluidRow(
    column(6,
           # Output: Additive space model ----
           plotOutput(outputId = "space_model")),
    column(6, 
           # Output: Niche model ----
           plotOutput(outputId = "niche_model"))
  ),
  
  fluidRow(
    box(
      column(4,
             titlePanel(
               # h5(paste0("Tracking data for species ", input$Species), align = 'left'))),  
               h4(textOutput(outputId = 'tracking_data_species_titulo'), align = 'center')
               
             )),
      column(8,
             # Output: map ----
             leafletOutput("mymap_species")),
    )
  ),
  
  
  fluidRow(
    
    box(
      title = "Data sharing options",  status = 'primary', solidHeader = TRUE, 
      collapsible = FALSE, collapsed = FALSE,
      uiOutput('tab')# )# ,
      # column(12,
      #        # Commented this out: ####
      #        # Input: add Identifier
      #        # textInput("Identifier", label = h4("Identifier:"), value = "", width = NULL, placeholder = "Type your first and last name"),
      #        # textInput("AHP", label = h4("Animal handling permit:"), value = "", width = NULL, placeholder = "Please add your animal handling permit"),
      #        # # status selector ----
      #        radioButtons("status", h4("Select a data sharing option:"), 
      #                     c("Metadata information only",
      #                       "Aggregated data visible only",
      #                       'Raw data visible only',
      #                       "Data downloadable"))
      # ),
      # # Commented this out ####           
      # # column(8,
      # #        # Output table ----
      # #        h4("Your data sharing permissions"),
      # #        dataTableOutput("table"),
      # #        fluidRow(
      # #          column(6, actionButton("submit", strong("Submit entries")), 
      # #                 style='padding-left:1px; padding-right:1px')
      # #        ),
      # #        hr(),
      # # Help box with instructions
      # fluidRow(
      #   column(12,
      #          actionButton("action", label = "Instructions", icon = icon("question-circle")),
      #          uiOutput("HelpBox")
      #          
      #   )
      # ) # Uncommenrted # End of column 8 
      
    ) # End of box
    
  ) # End of fluid row
  
) # Close Fluid page


server <- function(session, input, output){
  
  # add help box with instructions
  output$HelpBox = renderUI({
    if (input$action %% 2){
      helpText(HTML(
        
        # "Start by typing <I>your name</I> in the <b>Identifier</b> box, then select your <b>Animal handling permits</b> and </b>data sharing permission settings</b>. <br/>
        "Please fill out the Google Survey indicating the data sharing permission you wish to choose for your animal tracking data. <br/>
        In addition, we kindly ask you to provide information on the animal handling permits related to your animal tracking data. This information is necessary to participate in our study and to provide to journals during submission of our article.
      You can choose between the following options.
               -Public summary: A public summary of the study and the contact person identified as data owner               is provided. <br/>
              -Public tracks: Some or all tracks of your individual animals are visible to the public, but cannot be downloaded. In addition, only certain time periods of your animal tracks may be publicly visible, such as the last year of data. Such data sharing may be desired when data owners are still performing analysis on their data. <br/>
               -Public data: The public can view, and download the animal tracking data. Such data sharing is suggested when no endangered species are being worked with and when there are no institutional limitations about data sharing.<br/>
               -Private data: You may decide to not share any of your data, with only your summary information publicly available.<br/>
               "))
    } else {
      return()
    }
  })
  
  url <- a(HTML("</b>Submit your data sharing permissions here</b>"), href="https://docs.google.com/forms/d/e/1FAIpQLSdsDfN_UAHKUV-KhO2y5VNruls3WIfGrqvS6c3v6bEyjrk_MQ/viewform")
  output$tab <- renderUI({
    tagList(HTML("</b>URL link:</b>"), url)
  })
  
  
  # update species list option depending on checkbox
  observe({
    #browser()
    updateSelectInput(session, "Species",
                      label = "Species:",
                      choices =sort(unique(locs$taxon_canonical_name))
    )
  })
  
  # update species list option depending on checkbox
  observe({
    #browser()
    updateSelectInput(session, "StudyID",
                      label = "StudyID:",
                      choices = sort(unique((locs %>% filter(taxon_canonical_name == input$Species))$study_id))
    )                    
  })
  
  # To add:
  # Add study ID as a filter?!
  # update species list option depending on checkbox
  observe({
    #browser()
    updateSelectInput(session, "Individual",
                      label = "Individual:",
                      choices = sort(unique((locs %>% filter(taxon_canonical_name == input$Species & study_id == input$StudyID ))$individual_id))
    )                    
  })
  
  species_metadata = reactive({
    map_id <- locs %>% dplyr::filter(taxon_canonical_name == input$Species) 
    ddply(map_id,'taxon_canonical_name', function(x){
      data.frame(
        studies = length(unique(x$study_id)),
        inds = length(unique(x$individual_id)),
        n_locs = nrow(x),
        time_range_begin = range(x$timestamp)[1],
        time_range_end = range(x$timestamp)[2]
      )
    })
    
    
  })
  
  # subset data based on Individual input
  dat_id <- reactive({
    map_id <- locs %>% dplyr::filter(taxon_canonical_name == input$Species & study_id == input$StudyID & individual_id == input$Individual) 
    map_id <- map_id %>% arrange(as.Date(timestamp))
    map_id
  })
  
  # subset data using all individuals of a species and just display geohash lvl 3 and make the map of the USA:
  species_lvl_dat_id <- reactive({
    map_id <- locs %>% dplyr::filter(taxon_canonical_name == input$Species) 
    map_id <- map_id %>% arrange(as.Date(timestamp))
    map_id
  })
  
  
  # dbbmm 
  dbbmm_size_id <- reactive({
    map_id <- dbbmm_size %>% dplyr::filter(species == input$Species & study_id == input$StudyID & ind_id == input$Individual)
    map_id <- map_id %>% arrange(year, wk)
    map_id
  })
  
  niche_breadth_id <-  reactive({
    map_id <- niche_det %>% dplyr::filter(scientificname == input$Species & studyid == input$StudyID & individual == input$Individual)  %>% drop_na(total)
    map_id <- map_id %>% arrange(year, week)
    map_id
  })
  
  
  # Create geohash sf object for level 3
  geo_3 = reactive({
    dat <-  locs %>% dplyr::filter(taxon_canonical_name == input$Species & study_id == input$StudyID & individual_id == input$Individual) 
    geo_3_sf =   gh_to_sf( dat$geohash_3 )
    geo_3_sf$geohash_3 = rownames(geo_3_sf)
    n_sum_g3 = dat %>% group_by(geohash_3) %>% summarize(n_points =  n())
    geo_3_sf_2 = left_join(geo_3_sf, n_sum_g3, by = 'geohash_3')
    # paste0(str(geo_3_sf_2))  
    return(geo_3_sf_2)
  })
  
  
  # Create geohash object for level 5
  geo_5 = reactive({
    dat <- locs %>% dplyr::filter(taxon_canonical_name == input$Species & study_id == input$StudyID & individual_id == input$Individual) 
    geo_5_sf =   gh_to_sf( dat$geohash_5 )
    geo_5_sf$geohash_5 = rownames(geo_5_sf)
    n_sum_g5 = dat %>% group_by(geohash_5) %>% summarize(n_points =  n())
    geo_5_sf_2 = left_join(geo_5_sf, n_sum_g5, by = 'geohash_5')
    return(geo_5_sf_2)
    
  })
  
  
  geo_3_species = reactive({
    dat <- locs %>% dplyr::filter(taxon_canonical_name == input$Species) 
    geo_3_sf =   gh_to_sf( dat$geohash_3 )
    geo_3_sf$geohash_3 = rownames(geo_3_sf)
    n_sum_g3 = dat %>% group_by(geohash_3) %>% summarize(n_points =  n())
    geo_3_sf_2_species = left_join(geo_3_sf, n_sum_g3, by = 'geohash_3')
    return(geo_3_sf_2_species)
    
  })
  
  
  
  # Create geohash object for level 5
  # geo_6 = reactive({
  #   dat <- locs %>% dplyr::filter(taxon_canonical_name == input$Species & study_id == input$StudyID & individual_id == input$Individual) 
  #   geo_6_sf =   gh_to_sf( dat$geohash_6 )
  #   geo_6_sf$geohash_6 = rownames(geo_6_sf)
  #   n_sum_g6 = dat %>% group_by(geohash_6) %>% summarize(n_points =  n())
  #   geo_6_sf_2 = left_join(geo_6_sf, n_sum_g6, by = 'geohash_6')
  #   return(geo_6_sf_2)
  #   
  # })
  
  # 
  # int_model_space = list.files(indir_space_models, full.names = T, pattern = input$Species)
  # load( int_model_space) # get the newest:
  # area_int <- out$model
  # area_ghmq <- quantile(out$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)
  
  
  space_use_interactive = reactive({
    req(input$Species)
    if(file.exists(list.files(indir_space_interaction_models, full.names = T, pattern = input$Species))){  
      int_model_space = list.files(indir_space_interaction_models, full.names = T, pattern = input$Species)
      load( int_model_space) # get the newest:
      area_int <- out$model
      rm(out)
      return(area_int)
    }
  })
  
  # browser()
  space_use_interactive_data = reactive({
    req(input$Species)
    # tmp = area_mod_summary %>% filter(species == input$Species)
    # if(tmp$inter_sig == 'Y'){
    if(file.exists(list.files(indir_space_interaction_models, full.names = T, pattern = input$Species))){  
      int_model_space = list.files(indir_space_interaction_models, full.names = T, pattern = input$Species)
      load( int_model_space) # get the newest:
      # area_int_data <- data.frame(out$data$ghm_scale)
      area_ghmq <- quantile(out$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)
      rm(out)
      return(area_ghmq)
      
    }
  })
  
  # browser()
  space_use = reactive({
    req(input$Species) # requirement function: the input the user provides we require the user 
    # sp = species %in% input$taxa # Filter by the user defined taxa
    # Load the newest model:  
    single_sp_add_area_mod = list.files(indir_space_models, full.names = T, pattern = input$Species)
    
    # Add if else interactive or additive:
    # area_mod_summary %>% filter(inter_sig == 'Y')
    
    load(single_sp_add_area_mod[length(single_sp_add_area_mod)])  
    area_add <- out
    space_use <- area_add$model
    rm(out)
    return(space_use)
  })
  
  niche_breadth =  reactive({
    req(input$Species) # requirement function: the input the user provides we require the user 
    # sp = species %in% input$taxa # Filter by the user defined taxa
    # Load the newest model:  
    
    # Load the newest model:  
    single_sp_cont_niche_mod = list.files(indir_niche_models, full.names = T, pattern = input$Species)
    load(single_sp_cont_niche_mod[length(single_sp_cont_niche_mod)])  
    
    niche_breadth <- out$model
  })
  
  
  # Add significane statements:
  signif_add_ghm =  reactive({
    req(input$Species) # requirement function: the input the user provides we require the user 
    signif_tmp = area_mod_summary %>% dplyr::filter(species == input$Species) %>% dplyr::select(ghm_sig)
    
    if(signif_tmp$ghm_sig == 'Y'){
      signif <- 'Significant effect '
    }else{
      signif <- 'Non significant effect '
    }
    signif_add_ghm <- signif
  })
  
  signif_add_sg =  reactive({
    req(input$Species) # requirement function: the input the user provides we require the user 
    signif_tmp = area_mod_summary %>% dplyr::filter(species == input$Species) %>% dplyr::select(sg_sig)
    
    if(signif_tmp$sg_sig == 'Y'){
      signif <- 'Significant effect '
    }else{
      signif <- 'Non significant effect '
    }
    signif_add_sg <- signif
  })
  
  niche_add_sg = reactive({
    req(input$Species) # requirement function: the input the user provides we require the user 
    signif_tmp = niche_mod_summary %>% dplyr::filter(species == input$Species) %>% dplyr::select(cont_sg_sig)
    
    if(signif_tmp$cont_sg_sig == 'Y'){
      signif_niche_sg <- 'Significant effect '
    }else{
      signif_niche_sg <- 'Non significant effect '
    }
    niche_add_sg <- signif_niche_sg
  })
  
  niche_add_ghm  = reactive({
    req(input$Species) # requirement function: the input the user provides we require the user 
    signif_tmp = niche_mod_summary %>% dplyr::filter(species == input$Species) %>% dplyr::select(cont_ghm_sig)
    
    if(signif_tmp$cont_ghm_sig == 'Y'){
      signif_niche_ghm <- 'Significant effect '
    }else{
      signif_niche_ghm <- 'Non significant effect '
    }
    niche_add_ghm <- signif_niche_ghm
  })
  
  # if(area_mod_summary %>% filter(species == input$Species) %>% dplyr::select(ghm_sig) == 'Y'){
  #   if(area_mod_summary %>% filter(species == input$Species) %>% dplyr::select(sg_sig) == 'Y'){
  # 
  
  
  
  # Plot Space Use model:
  output$space_model <- renderPlot({
    # make_additive_area_model(wildlife(), indir_space_models, area_mod_summary)
    
    if(area_mod_summary %>% filter(species == input$Species) %>% dplyr::select(inter_sig) == 'Y'){
      
      # int_model_space = list.files(indir_space_models, full.names = T, pattern = input$Species)
      # load( int_model_space) # get the newest:
      # area_int <- out$model
      # area_ghmq <- quantile(out$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)
      
      ## get observed quantiles of ghm to set "low" and "high" human mod
      #  area_ghmq <- quantile(space_use_interactive_data(), probs = c(0.10, 0.90), na.rm = T)
      
      # Conditional Effects Plot for interaction
      area_ce_int <- conditional_effects(x=space_use_interactive(), 
                                         effects = "sg_norm:ghm_scale",
                                         int_conditions = list(ghm_scale = space_use_interactive_data()),
                                         re_formula = NA)
      pal <- c("#7552A3", "#CEBEDA")
      (area_int_ce_plot <-  plot(area_ce_int, plot = FALSE,
                                 line_args = list("se"=T,
                                                  "alpha" = 0.2))[[1]] +
          scale_color_manual(values = pal, name = "Modification",
                             labels = c("High", "Low")) +
          scale_fill_manual(values = pal, name = "Modification",
                            labels = c("High", "Low")) +
          xlab("Mobility") +
          ylab("Space Use")+
          theme_cowplot()  +
          theme(# legend.position = "none",
            legend.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            aspect.ratio = 1) +
          scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
          labs(title = 'Significant interactive space use model') + labs(subtitle = 'Conditional effects of modification and mobility from \n interaction models with significant effects on species \n space use. Low modification reflects the \n 10th pecentile and high modification \n the 90th percentile values of human modification') + theme(legend.position = 'bottom') + 
          NULL
      )
      # 
      # 
      # p + labs(title="Two separate labs calls for each title") +
      #   labs(subtitle=st) +
      #   theme(legend.position="bottom"),
      # 
      area_int_ce_plot
      
    }else{
      
      # If else: 
      
      # Make a variable that is significant or not significant and paste it on the header. ####
      # Add Home ranges ####
      
      
      #-- Mobility --#
      # Get conditional effects
      area_sg <- conditional_effects(x=space_use(),
                                     effects = "sg_norm",
                                     re_formula = NA)
      
      (area_sg_ce_plot <-  plot(area_sg, plot = F,
                                line_args = list("se" = T,
                                                 "color" = col,
                                                 "fill" = col))[[1]] +
          # scale_color_manual(values = palnew[3])+
          # theme_minimal() +
          xlab("Mobility") +
          ylab("Space Use")+
          # ggtitle("Puma concolor")+
          theme_cowplot()  +
          theme(legend.title = element_blank(),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 12),
                axis.text = element_text(size = 12),
                aspect.ratio = 1) +
          scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
          ggtitle(paste0('Additive space use model \n', signif_add_sg())) +  
          theme(legend.position = 'bottom') + 
          # labs(x  = "Day of year", y = "Mobility") +
          NULL
      )
      
      #-- Modification --#
      
      # Get conditional effects
      area_ghm <- conditional_effects(x=space_use(),
                                      effects = "ghm_scale",
                                      re_formula = NA)
      
      (area_ghm_ce_plot <-  plot(area_ghm, plot = F,
                                 line_args = list("se" = T,
                                                  "color" = col,
                                                  "fill" = col))[[1]] +
          # scale_color_manual(values = palnew[3])+
          # theme_minimal() +
          xlab("Modification") +
          ylab("Space Use")+
          # ggtitle("Puma concolor")+
          theme_cowplot()  +
          theme(legend.title = element_blank(),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 12),
                axis.text = element_text(size = 12),
                aspect.ratio = 1) +
          scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
          ggtitle(paste0('Additive space use model \n', signif_add_ghm())) +
          #  if(area_mod_summary %>% filter(species == input$Species) %>% dplyr::select(ghm_sig) == 'Y'){
          #  ggtitle(paste0(signif_add_ghm, ' additive space use model')) +
          # labs(x  = "Day of year", y = "Mobility") +
          NULL
        # }else{
        #   ggtitle(paste0('Not significant additive space use model')) +
        #     # labs(x  = "Day of year", y = "Mobility") +
        #     NULL
        #   
        # }
      )
      
      # Patchwork both together
      # area_sg_ce_plot + area_ghm_ce_plot + labs(subtitle = 'Additive effects of human modification and mobility \n on species space use ')
      
      # Add significance statements:       
      # if(area_mod_summary %>% filter(species == input$Species) %>% dplyr::select(ghm_sig) == 'Y'){
      #   if(area_mod_summary %>% filter(species == input$Species) %>% dplyr::select(sg_sig) == 'Y'){
      # 
      
      grid.arrange(area_sg_ce_plot,area_ghm_ce_plot,nrow=1,top=textGrob("Additive effects of human modification and mobility \n on species space use"))# , face = 'bold', size = 14))
      
    } # End of else
    
    
    #
    
    
    
    
  })
  
  # Plot Niche model:
  output$niche_model <- renderPlot({
    #---- Niche Plots ----#
    
    
    #-- Mobility --#
    
    # Get conditional effects
    niche_sg_cont <- conditional_effects(x=niche_breadth(),
                                         effects = "sg_norm",
                                         re_formula = NA)
    
    (niche_sg_ce_plot_cont <-  plot(niche_sg_cont, plot = F,
                                    line_args = list("se" = T,
                                                     "color" = col,
                                                     "fill" = col))[[1]] +
        # scale_color_manual(values = palnew[3])+
        # theme_minimal() +
        xlab("Mobility") +
        ylab("Niche Breadth")+
        # ggtitle("Puma concolor")+
        theme_cowplot()  +
        theme(legend.position = "none",
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 12),
              aspect.ratio = 1) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
        # ggtitle(paste0('Area controlled niche model \n')) +
        # ggtitle(paste0('Additive space use model \n', signif_add_sg())) +  
        ggtitle(paste0('Area controlled niche model \n', niche_add_sg() )) +
        # labs(x  = "Day of year", y = "Mobility") +
        NULL
    )
    
    #-- Modification --#
    
    # Get conditional effects
    niche_ghm_cont <- conditional_effects(x=niche_breadth(),
                                          effects = "ghm_scale",
                                          re_formula = NA)
    
    (niche_ghm_ce_plot_cont <-  plot(niche_ghm_cont, plot = F,
                                     line_args = list("se" = T,
                                                      "color" = col,
                                                      "fill" = col))[[1]] +
        # scale_color_manual(values = palnew[3])+
        # theme_minimal() +
        xlab("Modification") +
        ylab("Niche Breadth")+
        # ggtitle("Puma concolor")+
        theme_cowplot()  +
        theme(legend.title = element_blank(),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 12),
              aspect.ratio = 1) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
        # ggtitle(paste0('Area controlled niche model \n' )) +
        ggtitle(paste0('Area controlled niche model \n', niche_add_ghm() )) +
        # ggtitle(paste0('Area controlled niche model \n', niche_add_ghm() )) +
        #        geom_rug(data = dat_id(), aes(x = dat_id()$ghm),inherit.aes=F, col="steelblue",alpha=0.1, size=1.5) + 
        # labs(x  = "Day of year", y = "Mobility") +
        NULL
    )
    
    
    #  niche_sg_ce_plot_cont + niche_ghm_ce_plot_cont
    grid.arrange(niche_sg_ce_plot_cont,niche_ghm_ce_plot_cont,nrow=1,top=textGrob("Additive effects of human modification and mobility \n on species niche breadth")) ##,  face = 'bold', size = 14))
    
  })
  
  
  # Add map of tracking data:
  
  # This render fluid could be 3 interactive leaflets:
  # This works but not in ctreactive 
  # locs_sub = locs[20000:21000,]
  # geo_3_sf =   gh_to_sf( locs_sub$geohash_3 )
  # geo_3_sf$geohash_3 = rownames(geo_3_sf)
  # n_sum_g3 = locs[20000:21000,] %>% group_by(geohash_3) %>% summarize(n_points =  n())
  # geo_3_sf_2 = left_join(geo_3_sf, n_sum_g3, by = 'geohash_3')
  # 
  # pal = colorNumeric(palette = 'viridis', 
  #                    domain = geo_3_sf_2$n_points)
  # leaflet() %>%
  #   addProviderTiles(providers$Stamen.Terrain,
  #                    options = providerTileOptions(noWrap = TRUE)) %>%
  #   addPolygons(data = geo_3_sf_2, color = ~pal(geo_3_sf_2$n_points),popup = paste('Number of locs', geo_3_sf_2$n_points, "<br>")
  #   ) %>% 
  #   addScaleBar(position='topright',
  #               options=scaleBarOptions(maxWidth=200,imperial=FALSE))
  # 
  # 
  
  # Add an observer so everytime you add a layer 
  # browser()
  output$mymap <- renderLeaflet({
    
    #  if(input$Data_for_map %in% 'All data'){
    #        pal = colorFactor(palette = 'viridis', dat_id()$year)
    #       leaflet() %>%
    #         addProviderTiles(providers$Stamen.Terrain,
    #                          options = providerTileOptions(noWrap = TRUE)) %>%
    #        
    #   #      pal = colorNumeric(palette = 'viridis',domain = geo_5()$n_points) %>%
    #   # addPolygons(data = geo_5(), popup = paste('Number of locs', geo_5()$n_points, "<br>"), color = ~pal(geo_5()$n_points) ) %>%
    # #        addPolygons(data = geo_5(), color = ~pal(geo_5()$n_points) ) %>%
    #   # #
    #   # pal = colorNumeric(palette = 'viridis',domain = geo_3()$n_points) %>%
    # #  addPolygons(data = geo_3(), color = ~pal(geo_3()$n_points) ) %>%
    #         addCircleMarkers(data = dat_id(), lng = ~lon, lat = ~lat,
    #                          popup = paste('Date',dat_id()$date,"<br>",
    #                                        'Geohash 5',dat_id()$geohash_5,"<br>",
    #                                        'Geohash 3',dat_id()$geohash_3,"<br>",
    #                                        'Human modification',round(dat_id()$ghm,2),"<br>",
    #                                        'Human device count',dat_id()$safegraph_daily_count,"<br>"), 
    #                          #clusterOptions = markerClusterOptions(),
    #                          color = ~pal(year), # fill = ~col_level,
    #                          radius = 6, stroke = FALSE, fillOpacity = 0.8,
    #                          group = 'All data') %>% 
    #         addScaleBar(position='topright',
    #                     options=scaleBarOptions(maxWidth=200,imperial=FALSE)) %>%
    # pal <- colorNumeric(palette =c('purple','yellow'), domain = dat_id()$year)
    
    
    if(input$Data_for_map %in% c('Geohash 5', 'Geohash 3')){
      
      # pal = colorFactor(palette = 'viridis', dat_id()$year, discrete = TRUE)
      # pal <- c("purple", "yellow")
      # pal <- colorFactor(c("purple", "yellow"), domain = c(2019, 2020))
      # pal <- colorNumeric(palette =c('purple','yellow'), domain = c(2019, 2020))
      pal <- colorNumeric(palette =c('#440154FF','#FDE725FF'), domain = c(2019, 2020))
      cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}
      pal_ind_geo_3 <- colorNumeric(cm.cols1(100), domain=geo_3()$n_points)
      pal_ind_geo_5 <- colorNumeric(cm.cols1(100), domain=geo_5()$n_points)
      
      
      leaflet() %>%
        addProviderTiles(providers$Stamen.Terrain,
                         options = providerTileOptions(noWrap = TRUE)) %>%
        
        addCircleMarkers(data = dat_id(), lng = ~lon, lat = ~lat,
                         popup = paste('Date',dat_id()$date,"<br>",
                                       'Year', dat_id()$year,"<br>",
                                       #  'Geohash 6',dat_id()$geohash_6,"<br>",                                       
                                       'Geohash 5',dat_id()$geohash_5,"<br>",
                                       'Geohash 3',dat_id()$geohash_3,"<br>",
                                       'Human modification',round(dat_id()$ghm,2),"<br>",
                                       'Human device count',dat_id()$safegraph_daily_count,"<br>"), 
                         #clusterOptions = markerClusterOptions(),
                         color = ~pal(year), # fill = ~col_level,
                         radius = 6, stroke = FALSE, fillOpacity = 0.8,
                         group = 'All data') %>% 
        addScaleBar(position='topright',
                    options=scaleBarOptions(maxWidth=200,imperial=FALSE)) %>%
        # addPolygons(data = geo_6(), popup = paste('Number of locs', geo_6()$n_points, "<br>"),
        #             group = 'Geohash 6') %>%
        addPolygons(data = geo_5(), popup = paste('Number of locs', geo_5()$n_points, "<br>"),
                    group = 'Geohash 5',   fillColor = ~pal_ind_geo_5(n_points), fillOpacity = 0.5) %>%
        addPolygons(data = geo_3(), popup = paste('Number of locs', geo_3()$n_points, "<br>"),
                    group = 'Geohash 3',   fillColor = ~pal_ind_geo_3(n_points), fillOpacity = 0.5)   # %>%
      # addLegend("bottomright", pal = ~pal(year), values = ~dat_id()$year,
      #           title = "Year")
      
      
      
      
      # addLegend("bottomright", pal = ~pal(geo_3()$n_points), values = geo_3()$n_points,
      #           title = "N. locations geohash 3",
      #           labFormat = labelFormat(prefix = "$"),
      #           opacity = 1, group = 'group_1') %>% 
      # addLegend("bottomright", pal = ~pal(year), title="Year",
      #            values = dat_id()$year) # %>%
      # addLayersControl(overlayGroups = c("group_1","group_2"),
      #                  options = layersControlOptions(collapsed = FALSE))
    }else{
      # pal = colorFactor(palette = 'viridis', dat_id()$year)
      pal <- colorNumeric(palette =c('#440154FF','#FDE725FF'), domain = c(2019, 2020))
      leaflet() %>%
        addProviderTiles(providers$Stamen.Terrain,
                         options = providerTileOptions(noWrap = TRUE)) %>%
        
        addCircleMarkers(data = dat_id(), lng = ~lon, lat = ~lat,
                         popup = paste('Date',dat_id()$date,"<br>",
                                       #   'Geohash 6',dat_id()$geohash_6,"<br>",
                                       'Geohash 5',dat_id()$geohash_5,"<br>",
                                       'Geohash 3',dat_id()$geohash_3,"<br>",
                                       'Human modification',round(dat_id()$ghm,2),"<br>",
                                       'Human device count',dat_id()$safegraph_daily_count,"<br>"), 
                         #clusterOptions = markerClusterOptions(),
                         color = ~pal(year), # fill = ~col_level,
                         radius = 6, stroke = FALSE, fillOpacity = 0.8,
                         group = 'All data') %>% 
        addScaleBar(position='topright',
                    options=scaleBarOptions(maxWidth=200,imperial=FALSE))   # %>%
      #   addLegend("bottomright", pal = ~pal(year), values = ~dat_id()$year,
      #             title = "Year")# %>%
      # # addLegend("bottomright", pal = ~pal(year), title="Year",
      #           values = dat_id()$year)
    }
    #   }
  })
  
  
  
  observeEvent(input$Data_for_map, {
    
    # if('All data ' %in% input$Data_for_map) {
    #   leafletProxy('mymap', data = dat_id()) %>%
    #     showGroup('All data') %>%
    #     hideGroup('Geohash 5') %>%
    #     hideGroup('Geohash 3')
    # }
    
    
    if('Geohash 5' %in% input$Data_for_map) {
      leafletProxy('mymap', data = geo_5()) %>%
        showGroup('Geohash 5') %>%
        hideGroup('Geohash 3')
      # hideGroup('Geohash 6')
    }else{
      leafletProxy('mymap', data = geo_5()) %>%
        hideGroup('Geohash 5')
      
    }
    if('Geohash 3' %in% input$Data_for_map) {
      leafletProxy('mymap', data = geo_3()) %>%
        showGroup('Geohash 3') %>%
        hideGroup('Geohash 5')
      # hideGroup('Geohash 6')
    }else{
      leafletProxy('mymap', data = geo_3()) %>%
        hideGroup('Geohash 3')
    }
    # if('Geohash 6' %in% input$Data_for_map) {
    #   leafletProxy('mymap', data = geo_6()) %>%
    #     showGroup('Geohash 6') %>%
    #     hideGroup('Geohash 5')
    #   hideGroup('Geohash 3')
    # }else{
    #   leafletProxy('mymap', data = geo_6()) %>%
    #     hideGroup('Geohash 6')
    # }
  }, ignoreNULL = FALSE)
  
  #  qpal <- colorQuantile("YlOrRd", domain = geo_3_species()$n_points, n=7)
  
  # pal <- reactive({
  #   colorQuantile("YlGn", geo_3_species()$n_points, n = 5)
  # })
  
  # pal <- colorNumeric(palette =c('purple','yellow'), domain = species_lvl_dat_id()$year)
  output$mymap_species <- renderLeaflet({
    
    # pal <- reactive({
    #   colorNumeric(
    #     palette = cm.cols1(100),
    #     domain = geo_3_species()$n_points
    #   )
    # })
    # pal <- colorNumeric("Reds", domain=species_lvl_dat_id()$n_points)
    cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}
    pal_species_geo_3 <- colorNumeric(cm.cols1(100), domain=species_lvl_dat_id()$n_points)
    # pal = colorNumeric(palette = 'viridis',domain = geo_3()$n_points) %>%
    leaflet() %>%
      addProviderTiles(providers$Stamen.Terrain,
                       options = providerTileOptions(noWrap = TRUE)) %>%
      
      addScaleBar(position='topright',
                  options=scaleBarOptions(maxWidth=200,imperial=FALSE)) %>%
      addPolygons(data = geo_3_species(), popup = paste('Number of locs', geo_3_species()$n_points, "<br>"),
                  group = 'Geohash 3',   fillColor = ~pal_species_geo_3(n_points), fillOpacity = 0.5 )    #  %>%
    # addLegend("bottomright", pal = ~pal(year), values = ~species_lvl_dat_id()$year,
    #           title = "Year")# %>%
  })
  
  
  
  # observe(
  #   
  #   leafletProxy('mymap', data = geo_5()) %>% clearMarkers() %>% 
  #   # pal = colorNumeric(palette = 'viridis',domain = geo_5()$n_points) %>% 
  #         addPolygons(data = geo_5(),
  #                     popup = paste('Number of locs', geo_5()$n_points, "<br>"))#,
  #                     # color = ~pal(geo_5()$n_points))
  # )
  
  # if(input$Data_for_map %in% 'Geohash'){
  #   leafletProxy('mymap') %>% 
  #     pal = colorNumeric(palette = 'viridis',domain = geo_5()$n_points) %>% 
  #       addPolygons(data = geo_5(),
  #                   popup = paste('Number of locs', geo_5()$n_points, "<br>"),
  #                   color = ~pal(geo_5()$n_points)            )
  # }
  # 
  # Why is this not working ####
  # if(input$Data_for_map %in% 'Geohash 5'){
  #    pal = colorNumeric(palette = 'viridis',
  #                       domain = geo_5()$n_points)
  # 
  #    leaflet() %>%
  #      addProviderTiles(providers$Stamen.Terrain,
  #                       options = providerTileOptions(noWrap = TRUE)) %>%
  # addPolygons(data = geo_5(),
  # popup = paste('Number of locs', geo_5()$n_points, "<br>"),
  # color = ~pal(geo_5()$n_points)
  #     ) %>%
  #     addScaleBar(position='topright',
  #                 options=scaleBarOptions(maxWidth=200,imperial=FALSE))
  # }
  # # # Why is this not working ####
  # if(input$Data_for_map == 'Geohash 3'){
  #   pal = colorNumeric(palette = 'viridis',
  #                      domain = geo_3()$n_points)
  # 
  #   leaflet() %>%
  #     addProviderTiles(providers$Stamen.Terrain,
  #                      options = providerTileOptions(noWrap = TRUE)) %>%
  #     addPolygons(data = geo_3(),
  #                 popup = paste('Number of locs', geo_3()$n_points, "<br>"),
  #                 color = ~pal(geo_3()$n_points)
  #                 ) %>%
  #     addScaleBar(position='topright',
  #                 options=scaleBarOptions(maxWidth=200,imperial=FALSE))
  # }
  # Add polygons shapefile dbmm
  # })
  
  # observe({
  #   pal = colorNumeric(palette = 'viridis',  domain = geo_5()$n_points)    
  #   leafletProxy("map", data = geo_5()) %>%
  #     clearShapes() %>%
  #     addPolygons(data = geo_3(),
  #                 popup = paste('Number of locs', geo_3()$n_points, "<br>"),
  #                 color = ~pal(geo_3()$n_points)
  #     )
  # })
  # 
  cols <- c('2019' = "#440154FF", '2020' = '#FDE725FF')
  
  # DBBMM size individual
  output$dbbmm_niche_week_size_id <- renderPlot({
    dbmm_week = ggplot( data = dbbmm_size_id(), aes(x = wk, y = (area / 1000000), color = col_level, fill = col_level)) + geom_point()  +
      labs(y = 'Space use in km2', x = 'Week of the year') +
      # scale_color_viridis(discrete = TRUE) +
      scale_colour_manual(values = cols) +
      theme_classic() + theme(axis.title = element_text(face = "bold", size = 14),
                              legend.position = c(0.90, 0.15), legend.title = element_blank())#  +
      # ggtitle('Weekly space use \n of an individual animal \n through time')
  # })
  
  # output$niche_week_id <- renderPlot({
    niche_week = ggplot( data = niche_breadth_id(), aes(x = week, y = total, color = col_level, fill = col_level)) + geom_point()  +
      labs(y = 'Total niche breadth', x = 'Week of the year') +
      # scale_color_viridis(discrete = TRUE) +
      scale_colour_manual(values = cols) +
      theme_classic() + theme(axis.title = element_text(face = "bold", size = 14),
                              legend.position="none", legend.title = element_blank()) #+
    
    
    grid.arrange(dbmm_week,niche_week,nrow=1,top=textGrob("Modeling responses: \n Weekly space use and niche breadth \n of an individual animal"))
    
     # ggtitle('Weekly niche breadth \n of an individual animal \n through time')
  })
  
  
  
  
  
  # Scatterplot GHM
  output$indiv_scatterplot_ghm_sg <- renderPlot({
    
    
    indiv_ghm = ggplot( data = dat_id(), aes(x = jday, y = ghm, color = col_level, fill = col_level)) + geom_point()  +
      labs(y = 'Global Human Modification', x = 'Julian day') +
      # scale_color_viridis(discrete = TRUE) +
      scale_colour_manual(values = cols) +
      theme_classic() + theme(axis.title = element_text(face = "bold", size = 14),
                              legend.position="none", legend.title = element_blank())
    
    indiv_sg = ggplot( data = dat_id(), aes(x = jday, y = safegraph_daily_count, color = col_level, fill = col_level)) + geom_point() +
      labs(y = 'Daily cellphone count', x = 'Julian day') +
      # scale_color_viridis(discrete = TRUE) +
      scale_colour_manual(values = cols) +
      theme_classic() + theme(axis.title = element_text(face = "bold", size = 14),
                              legend.position = c(0.90, 0.15),
                              legend.title = element_blank())
    
    grid.arrange(indiv_sg,indiv_ghm,nrow=1,top=textGrob("Static and dynamic human activities \n experienced by an individua animal"))
    
    
    
  })
  # 
  # # Scatterplot GHM
  # output$indiv_scatterplot_sg <- renderPlot({
  #   ggplot( data = dat_id(), aes(x = jday, y = safegraph_daily_count, color = col_level, fill = col_level)) + geom_point() +
  #     labs(y = 'Daily cellphone count', x = 'Julian day') +
  #     # scale_color_viridis(discrete = TRUE) +
  #     scale_colour_manual(values = cols) +
  #     theme_classic() + theme(axis.title = element_text(face = "bold", size = 14),
  #                             legend.position = c(0.90, 0.15),
  #                             legend.title = element_blank()) # + 
  #   # ggtitle('Human mobility experienced \n by an individual animal \n through time')# + guides(color = FALSE)
  # })
  
  # 
  
  
  
  output$titulo = renderText({
    # req(input$Species)
    # species = unique(input$Species)
    paste0("Model output for study species ", input$Species,'\n (', species_metadata()$studies, ' studies ', species_metadata()$inds, ' individuals)' )
    
    # output$txtOutput = renderText({
    #   paste0("The area of the circle is: ", pi*input$numInput^2)
    # 
    # return(title)
  })
  
  output$tracking_data_species_titulo  = renderText({
    # req(input$Species)
    # species = unique(input$Species)
    paste0("Tracking data for species  ", input$Species,'\n (', species_metadata()$studies, ' studies ', species_metadata()$inds, ' individuals)' )
    
    # output$txtOutput = renderText({
    #   paste0("The area of the circle is: ", pi*input$numInput^2)
    # 
    # return(title)
  })
  
  
  # FIX THIS 
  # output$metadata = renderText({
  #   paste0(species_metadata() )
  # })
  
  
  
} # End of server 

# Add map and 2 plots to explore individuals in detail: GHM and Safegraph through time colour coded by year:

# output$species_ID <- renderText({
#   temp_nseq <- (nseq %>% 
#                   filter(subproject_name == input$subproj))$nseq
#   paste("Total number of sequences:", temp_nseq)
# })


# output$species_id = renderText({
#   req(input$Species)
#   species = unique(input$Species)
#   return(species)
#   paste0('')
# })

values <- reactiveValues()
values$dat <- data.frame(Identifier = factor(),
                         AHP = factor(),
                         # Individual = numeric(),
                         Status = factor()
)

observeEvent(input$submit, {
  new_entry <- data.frame(Indentifier = input$Identifier,
                          AHP = input$AHP,
                          Status = factor(input$Status), 
                          levels = c("Metadata information only",
                                     "Aggregated data visible only",
                                     "Raw data visible only",
                                     "Data downloadable")
  )
  #  browser()
  # observeEvent(input$submit, {
  # new_entry <- data.frame(Indentifier = unique(values$dat$Identifier),
  #                              AHP = unique(values$dat$AHP),
  #                              Status = factor(values$dat$Status), 
  #                              levels = c("Metadata information only",
  #                                         "Aggregated data visible only",
  #                                         "Raw data visible only",
  #                                         "Data downloadable")
  # )
  
  # add new entry to dataframe
  values$dat <- rbind(values$dat, new_entry) 
  
  
  update_counter <- data.frame(Indentifier = unique(values$dat$Identifier),
                               AHP = unique(values$dat$AHP),
                               Status = unique(values$dat$Status))
  
  
  sheet_append(ss = SHEET_PATH, update_counter)
  
  showModal(modalDialog(
    title = "Data sharing agreement uploaded succesfully",
    easyClose = TRUE,
    footer = NULL
  ))
  
})

#Date = Date())

# observeEvent(input$submit, {

# Create a unique file name
# filename <- paste0(input$Species, "_", input$StudyID,"_", Sys.time(), "_.csv")
# filename <- gsub(" ", "_", filename, fixed=TRUE)
# filename <- gsub(":", "_", filename, fixed=TRUE)
# # Write the data to a temporary file locally
# filePath <- file.path(tempdir(), filename, fsep = "\\")
# write.csv(values$dat, filePath)
# # Upload the file to Google Drive
# drive_upload(media = filePath, path = DRIVE_PATH, name = filename)

# Keep track of which individuals have been processed
# extract desired information
# update_counter <- data.frame(Indentifier = unique(values$dat$Identifier),
#                              Species = unique(values$dat$Species),
#                              Individual = unique(values$dat$Individual))
# # add the data as new rows
# sheet_append(ss = SHEET_PATH, update_counter)
# 
# pop-up message

# observeEvent(input$submit, {
#                                        
#   # add the data as new rows
#   sheet_append(ss = SHEET_PATH, update_counter)
# 
# showModal(modalDialog(
#   title = "Data sharing agreement uploaded succesfully",
#   easyClose = TRUE,
#   footer = NULL
# ))
# 
# remove all rows after submitting
# empty_df <- values$dat[0, ]
# values$dat <- empty_df
#   
# })

# }

# 
# ui <- basicPage(
#   h1('Model output for study species'),
#   selectInput(inputId = 'taxa',
#               label = 'Choose wildlife',
#               'Names'),
#   plotOutput('plot')
# )

shinyApp(ui = ui, server = server)