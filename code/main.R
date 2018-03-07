# integrated model with data on individual growth and population
#   size distributions

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
source("./code/install_packages.R")

# source helper functions
source("./code/construct_fitted_matrix.R")
source("./code/helpers.R")
source("./code/mpm_est.R")
source("./code/param_modules.R")
source("./code/sim_dynamics.R")
source("./code/otolith_calc.R")

# load data
alldat <- read.csv("./data/VEFMAP_FISH_20171024.csv")

# load snags data set
snags_data <- read.csv("./data/SNAGS_FISH_20171205.csv")
snags_data$date_new <- format(dmy_hms(snags_data$surveydate), format = "%d/%m/%Y")
snags_data$YEAR <- sapply(strsplit(snags_data$date_new, "/"),
                          function(x) x[3])
snags_data$taxonname <- as.character(snags_data$taxonname)
snags_data$taxonname <- ifelse(snags_data$taxonname == "Yellowbelly",
                               "Golden perch",
                               snags_data$taxonname)
snags_data$taxonname <- factor(snags_data$taxonname)
snags_data2 <- data.frame(SYSTEM = rep("LOWERMURRAY", nrow(snags_data)),
                          SITE_CODE = paste0("Lm", snags_data$idsite),
                          Reach = rep(1, nrow(snags_data)),
                          geartype = factor(rep("EF/Boat"), nrow(snags_data)),
                          Event_Date = snags_data$date_new,
                          Pass.No = rep(1, nrow(snags_data)),
                          total_no_passes = rep(1, nrow(snags_data)),
                          seconds = snags_data$seconds,
                          Common.Name = snags_data$taxonname,
                          Scientific.Name = snags_data$Scientific.Name,
                          totallength = snags_data$totallength,
                          WEIGHT = snags_data$weight,
                          Total.Sampled = rep(1, nrow(snags_data)),
                          VEFMAP.Stage = rep(NA, nrow(snags_data)),
                          YEAR = as.integer(snags_data$YEAR))

# load ovens data and combine with alldat
ovens_data <- read.table("./data/vba_ovens_2008_2017.csv", sep = "\t", header = TRUE)
ovens_data$date_new <- format(dmy(ovens_data$date), format = "%d/%m/%Y")
ovens_data$YEAR <- sapply(strsplit(ovens_data$date_new, "/"),
                          function(x) x[3])
ovens_data$species <- as.character(ovens_data$species)
ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii ",
                             "Maccullochella peelii peelii",
                             ovens_data$species)
ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii",
                             "Maccullochella peelii peelii",
                             ovens_data$species)
ovens_data$species <- factor(ovens_data$species)
ovens_data$common_name <- alldat$Common.Name[match(ovens_data$species, alldat$Scientific.Name)]
ovens_data2 <- data.frame(SYSTEM = rep("OVENS", nrow(ovens_data)),
                          SITE_CODE = paste0("Ov", ovens_data$site),
                          Reach = rep(1, nrow(ovens_data)),
                          geartype = ovens_data$gear_type,
                          Event_Date = ovens_data$date_new,
                          Pass.No = rep(1, nrow(ovens_data)),
                          total_no_passes = rep(1, nrow(ovens_data)),
                          seconds = ovens_data$electro_seconds,
                          Common.Name = ovens_data$common_name,
                          Scientific.Name = ovens_data$species,
                          totallength = ovens_data$total_length_mm,
                          WEIGHT = ovens_data$weight_g,
                          Total.Sampled = ovens_data$no_collected,
                          VEFMAP.Stage = rep(NA, nrow(ovens_data)),
                          YEAR = as.integer(ovens_data$YEAR))
alldat <- rbind(alldat, ovens_data2, snags_data2)

# load flow data
# source("./code/load-flow-data.R")

# set systems of interest
system_sub <- c("BROKEN",
                "LOWERMURRAY",
                "CAMPASPE",
                "GLENELG",
                "GOULBURN",
                "LODDON",
                "THOMSON",
                "OVENS")
alldat <- alldat[-which(is.na(match(alldat$SYSTEM, system_sub))), ]

# clean up common names
alldat$Common.Name <- tolower(alldat$Common.Name)
alldat$Common.Name <- gsub(" ", "", alldat$Common.Name)
alldat$Common.Name <- gsub("-[[:digit:]]*", "", alldat$Common.Name)
alldat$Common.Name <- gsub("/", "", alldat$Common.Name)
alldat$Common.Name <- gsub("sp\\.", "", alldat$Common.Name)
alldat$Common.Name <- gsub("flatheadedgudgeon", "flatheadgudgeon", alldat$Common.Name)
alldat$Common.Name <- gsub("europeancarp", "carp", alldat$Common.Name)
alldat$Common.Name <- gsub("redfinperch", "redfin", alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "weatherloach",
                             "orientalweatherloach",
                             alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "rainbowfish",
                             "murrayriverrainbowfish",
                             alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "hardyhead",
                             "unspeckedhardyhead",
                             alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "blackfish",
                             "riverblackfish",
                             alldat$Common.Name)

# remove spp not of interest
sp_to_rm <- c("carpgoldfishhybrid",
              "carpgudgeoncomplex",
              "commonyabby",
              "easternsnakeneckedturtle",
              "freshwatershrimp",
              "glenelgspinyfreshwatercrayfish",
              "goldfish",
              "hypseleostris",
              "murrayspinycrayfish",
              "nofish",
              "unidentifiedcod",
              "unknown",
              "",
              "codspp.",
              "yabbie",
              "longneckturtle",
              "gudgeons",
              "murraycray",
              "murraycodtroutcodhybrid")
alldat <- alldat[-which(!is.na(match(alldat$Common.Name, sp_to_rm))), ]
if (any(is.na(alldat$Common.Name)))
  alldat <- alldat[-which(is.na(alldat$Common.Name)), ]

# fill missing weights
sp_tmp <- unique(alldat$Common.Name)
for (i in seq_along(sp_tmp)) {
  
  # subset to single species
  dat_tmp <- alldat[which(alldat$Common.Name == sp_tmp[i]), ]
  
  # check to make sure some weights are missing      
  if (any(is.na(dat_tmp$WEIGHT))) {
    
    # use species-specific equation if it exists, generic otherwise
    if (length_weight_conversion[[sp_tmp[i]]]$n) {
      coefs <- length_weight_conversion[[sp_tmp[i]]]$coef
    } else {
      coefs <- length_weight_conversion$generic
    }
    
    # subset NA observations
    na_sub <- which(is.na(dat_tmp$WEIGHT))
    
    # estimate weight from length
    dat_tmp$WEIGHT[na_sub] <- exp(coefs[1] + coefs[2] * log(dat_tmp$totallength[na_sub]))
    
    # return estimated data to full data set
    alldat[which(alldat$Common.Name == sp_tmp[i]), ] <- dat_tmp
    
  }
}

# reformat dates
alldat$Date <- dmy(alldat$Event_Date)

# arrange alldat to have more consistent set of variables
alldat <- data.frame(date = alldat$Date,
                     year = alldat$YEAR,
                     site = alldat$SITE_CODE,
                     system = alldat$SYSTEM,
                     reach = alldat$Reach,
                     species = alldat$Common.Name,
                     length = alldat$totallength,
                     weight = alldat$WEIGHT,
                     abundance = alldat$Total.Sampled,
                     intensity = (alldat$total_no_passes * alldat$seconds))

# set up species names for plots
sp_names <- data.frame(full = c("Australian bass",
                                "Australian grayling",
                                "Australian smelt",
                                "Black bream",
                                "Bony bream",
                                "Brown trout",
                                "Carp",
                                "Carp gudgeon",
                                "Common galaxias",
                                "Dwarf flathead gudgeon",
                                "Eastern gambusia",
                                "Eel",
                                "Estuary perch",
                                "Flathead galaxias",
                                "Flathead gudgeon",
                                "Golden perch",
                                "Long-finned eel",
                                "Luderick",
                                "Mountain galaxias",
                                "Murray cod",
                                "Murray river rainbowfish",
                                "Obscure galaxias",
                                "Oriental weather loach",
                                "Pouched lamprey",
                                "Pygmy perch",
                                "Rainbow trout",
                                "Redfin",
                                "River blackfish",
                                "River garfish",
                                "Roach",
                                "Sea mullet",
                                "Short finned eel",
                                "Short headed lamprey",
                                "Silver perch",
                                "Southern pygmy perch",
                                "Tench",
                                "Trout cod",
                                "Tupong",
                                "Two spined blackfish",
                                "Unspecked hardyhead",
                                "Variegated pygmy perch",
                                "Western carp gudgeon",
                                "Yarra pygmy perch",
                                "Yellow eyed mullet"),
                       code = as.character(levels(alldat$species)))

# clean up NAs in reach column
if (any(is.na(alldat$reach))) {
  alldat <- alldat[-which(is.na(alldat$reach)), ]
}

# subset species
sp_sub <- c("goldenperch", "murraycod",
            "troutcod", "silverperch")
alldat <- alldat[which(!is.na(match(alldat$species, sp_sub))), ]

# add length/weight to age calculations
alldat$age <- rep(NA, nrow(alldat))
for (i in seq_along(unique(alldat$species))) {
  
  row_sub <- which(alldat$species == unique(alldat$species)[i])
  if (any(rownames(params_all) == unique(alldat$species)[i])) {
    alldat$age[row_sub] <- length_to_age(alldat$length[row_sub],
                                         params = params_all[which(rownames(params_all) == unique(alldat$species)[i]), ])
  } else {
    alldat$age[row_sub] <- rep(NA, length(row_sub))
  }
  
}

# do this for all ages, then just focus on 1-3 year olds for some analyses
catch_curve_fun <- function(x, sp) {
  
  x_sub <- x[which(x$species == sp), ]
  
  bins <- c(0, seq_len(max(x_sub$age, na.rm = TRUE) + 1)) + 0.5
  # bins <- c(bins[1:7], bins[length(bins)])
  obs_unique <- paste(x_sub$system, paste0("reach", x_sub$reach), x_sub$year, sep = "_")
  
  obs_vec <- unique(obs_unique)
  
  out <- matrix(NA, nrow = length(obs_vec), ncol = (length(bins) - 1))
  for (i in seq_along(obs_vec)) {
    dat_sub <- x_sub[which(obs_unique == obs_vec[i]), ]
    out[i, ] <- hist(dat_sub$age, breaks = bins, plot = FALSE)$counts
  }
  colnames(out) <- seq_len(ncol(out))
  rownames(out) <- obs_vec
  site_split <- strsplit(obs_vec, "_")
  site_info <- data.frame(system = sapply(site_split, function(x) x[1]),
                          reach = sapply(site_split, function(x) as.numeric(substr(x[2],
                                                                                   start = nchar(x[2]),
                                                                                   stop = nchar(x[2])))),
                          year = sapply(site_split, function(x) as.numeric(x[3])))
  
  out <- list(age_dist = out,
              info = site_info,
              bins = bins)
  
}

# do this for all ages, then just focus on 1-3 year olds for some analyses
size_dist_fun <- function(x, sp, nbins = 5) {
  
  x_sub <- x[which(x$species == sp), ]
  
  bins <- c(0, 100, 250, 400, 600, 1000, 1500, 15000, max(x_sub$weight, na.rm = TRUE))
  # bins <- exp(seq(0, log(max(x_sub$weight, na.rm = TRUE)), length = nbins))
  # bins <- quantile(x_sub$weight, p = seq(0, 1, length = nbins), na.rm = TRUE)
  # bins[1] <- 0
  # bins <- seq(0, max(x_sub$weight, na.rm = TRUE), length = nbins)
  obs_unique <- paste(x_sub$system, paste0("reach", x_sub$reach), x_sub$year, sep = "_")
  
  obs_vec <- unique(obs_unique)
  
  out <- matrix(NA, nrow = length(obs_vec), ncol = (length(bins) - 1))
  for (i in seq_along(obs_vec)) {
    dat_sub <- x_sub[which(obs_unique == obs_vec[i]), ]
    out[i, ] <- hist(dat_sub$weight, breaks = bins, plot = FALSE)$counts
  }
  colnames(out) <- seq_len(ncol(out))
  rownames(out) <- obs_vec
  site_split <- strsplit(obs_vec, "_")
  site_info <- data.frame(system = sapply(site_split, function(x) x[1]),
                          reach = sapply(site_split, function(x) as.numeric(substr(x[2],
                                                                                   start = nchar(x[2]),
                                                                                   stop = nchar(x[2])))),
                          year = sapply(site_split, function(x) as.numeric(x[3])))
  
  out <- list(size_dist = out,
              info = site_info,
              bins = bins)
  
}

mc_catch_curve <- catch_curve_fun(alldat, sp = "murraycod")
mc_size_tmp1 <- size_dist_fun(alldat, sp = "murraycod", nbins = 5)
mc_size_tmp2 <- size_dist_fun(alldat, sp = "murraycod", nbins = 6)
mc_size_tmp3 <- size_dist_fun(alldat, sp = "murraycod", nbins = 12)
tc_catch_curve <- catch_curve_fun(alldat, sp = "troutcod")
gp_catch_curve <- catch_curve_fun(alldat, sp = "goldenperch")
sp_catch_curve <- catch_curve_fun(alldat, sp = "silverperch")

mc_size1 <- list(t(mc_size_tmp1$size_dist[grep("MURRAY", mc_size_tmp1$info$system), ]),
                 t(mc_size_tmp1$size_dist[grep("OVENS", mc_size_tmp1$info$system), ]),
                 t(mc_size_tmp1$size_dist[grep("BROKEN_reach3", rownames(mc_size_tmp1$size_dist)), ]),
                 t(mc_size_tmp1$size_dist[grep("BROKEN_reach5", rownames(mc_size_tmp1$size_dist)), ]))
mc_size2 <- list(t(mc_size_tmp2$size_dist[grep("MURRAY", mc_size_tmp2$info$system), ]),
                 t(mc_size_tmp2$size_dist[grep("OVENS", mc_size_tmp2$info$system), ]),
                 t(mc_size_tmp2$size_dist[grep("BROKEN_reach3", rownames(mc_size_tmp2$size_dist)), ]),
                 t(mc_size_tmp2$size_dist[grep("BROKEN_reach5", rownames(mc_size_tmp2$size_dist)), ]))
mc_size3 <- list(t(mc_size_tmp3$size_dist[grep("MURRAY", mc_size_tmp3$info$system), ]),
                 t(mc_size_tmp3$size_dist[grep("OVENS", mc_size_tmp3$info$system), ]),
                 t(mc_size_tmp3$size_dist[grep("BROKEN_reach3", rownames(mc_size_tmp3$size_dist)), ]),
                 t(mc_size_tmp3$size_dist[grep("BROKEN_reach5", rownames(mc_size_tmp3$size_dist)), ]))

mc_oti <- read.csv("./data/MC.csv")
mc_oti_new <- data.frame(growth = mc_oti$Otolith_growth,
                         id = mc_oti$ID,
                         year = mc_oti$Year)  

# fit model with correct density dependence form
size_mod_mc_nodens1 <- estimate_mpm(pop_samp = mc_size1,
                                    site = seq_len(length(mc_size1)),
                                    growdat = mc_oti_new,
                                    mat_type = "stage",
                                    dens_depend = "none",
                                    greta_settings = list(nsamples = 10,
                                                          nwarmup = 10,
                                                          inits = "random"))
