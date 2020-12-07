
library(raster)
library(sp)
library(stringr)
library(gdistance)
library(RSAGA)
library(rgdal)

rm(list=ls())

############## START OF PARAMETERS ##############
#Directories and files settings
base_dir              <- "G:/gregalis/"
target_dir            <- paste(base_dir, "res1/", sep="")  # where the lineage model grids and working data go
SDM_model_base        <- paste(base_dir, "sdm/", sep="")
SDM_model             <- paste(SDM_model_base, "inverse_gregalis-05.asc", sep="")
cost_in               <- paste(SDM_model_base, "inverse_gregalis-05.asc", sep="")
lineage_site_filename <- paste(base_dir, "tabs/", "clades.csv", sep="")
lineage_field_name    <- "Clade"   # the column for lineage name in the site data

aa <- rsaga.set.env(workspace="C://Program Files/saga-7.2.0_x64/", cmd="saga_cmd.exe", path="C://Program Files/saga-7.2.0_x64/", modules="C://Program Files/saga-7.2.0_x64//tools", version="7.2.0", parallel=FALSE)
aa$version = rsaga.get.version(env=aa)

#Various filter settings
lin_exclude_list      <- c()   # this list allows for skipping at the lineage level
use_use               <- 0     # use column "use" for records filter ("1") or not ("0")

#Grid settings
grid_resolution       <- 2000
spRef                 <- "+proj=moll +lon_0=30 +x_0=3335846.22854 +y_0=-336410.83237 +datum=WGS84 +units=m +no_defs"
output_model_extent   <- extent(3797130.398597037, 11593130.39859704, 3685333.673568079, 7773333.673568079)    # the maximum extent for all lineage models

#Points settings
ppRef                 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#Model settings
buffer_dist         <- 5000    # the buffer distance in map units (presumably decimal degrees)
additional_buffer   <- 0      # how much (as a proportion) the output grids should extend beyond the buffered points
distance_method     <- "model-cost"           # determines whether distance is calculated as euclidean or model-weighted cost distance
## so far, can be "euclidian" or "model-cost"
weight_function     <- "inverse"    ## determines whether lineage weight is calculated as 1/distance or 1/(distance^2), or simply closest distance
## so far, can be "inverse" or "inverse_square", "inverse_cube", "inverse_quad"

min_SDM_value       <- 0.00005  # values below this for the SDM are set to zero, to simplify calculation in areas of essentially unsuitable habitat
scale_to            <- "one"      ## determines whether lineage weights within a model group sum to the SDM value or to 1 or no weights at all
## can be "model" or "one" or "none"


############### END OF PARAMETERS ###############

cat ("\n\n*************************************** \n Lineage Distribution Estimation Tool ")
cat ("\n    Dan Rosauer, with changes made by Andrey Lissovsky \n    September 2012 - February 2020")
cat ("\n***************************************\n")


time_start <- Sys.time()

  # Load the sequence site data
  cat ("Loading the lineage locations...\t")
  all_lineages_sites <- read.csv(lineage_site_filename)
  orig_rowcount <- nrow(all_lineages_sites)

  # filter the records
  excluded_lineages <- grep("excluded", all_lineages_sites[, lineage_field_name])  # remove records with 'excluded' in the lineage name
  if (length(excluded_lineages >0)) {all_lineages_sites <- all_lineages_sites[-excluded_lineages, ]}
  if (use_use == 1) {  all_lineages_sites <- all_lineages_sites[which(all_lineages_sites$Use %in% c(-1,1)), ]}    # remove records where 'use' <> -1 or 1
  all_lineages_sites <- all_lineages_sites[which(is.finite(all_lineages_sites[, "latitude"]) & is.finite(all_lineages_sites$longitude)), ]    # remove without a numeric latitude and longitude

  #all_lineages_sites$model_group <- str_trim(all_lineages_sites$model_group)
  all_lineages_sites[, lineage_field_name] <- str_trim(all_lineages_sites[, lineage_field_name])

  groupLineages     <- unique(all_lineages_sites[, lineage_field_name])

  cat (orig_rowcount, "rows read,", nrow(all_lineages_sites), "valid records loaded\n")

  # turn lineage points into a SpatialPoints object
  allLineagePoints <- SpatialPointsDataFrame(all_lineages_sites[, c("longitude", "latitude")], data=all_lineages_sites, proj4string = CRS(ppRef))
  allLineagePoints <- spTransform(allLineagePoints, spRef)
  allLineagePoints <- crop(allLineagePoints, output_model_extent)

    # load the SDM model raster
  #********************************?*******************************
    SDM.ras <- raster(SDM_model)
    projection(SDM.ras) <- CRS(spRef)
    group_extent <- extent(SDM.ras)

    # get a list of the lineages in this group
    lineages <- groupLineages
    if (length(lin_exclude_list) > 0) {
      lineages <- setdiff(lineages,lin_exclude_list)
    }

    cat ("\nLineages in", ":")
    for (lineage in lineages) {
      cat ("\n   ", lineage)
    }

    if (length(lineages) > 1) { # proceed with lineage models if there are multiple
                                # lineages - otherwise just copy the SDM for the model group

       thisGroupPoints <- allLineagePoints
       thisGroupPoints <- crop(thisGroupPoints, group_extent)

      # calculate a new extent
      group_points_extent <- extent(thisGroupPoints)
      buffer_ratio  <- 1 + additional_buffer
      extent_buffer <- buffer_dist * buffer_ratio

      # new extent is the same as points layer + a buffer, but where the extended buffer
      # goes beyond the extent of the maxent model, limit to the output model extent.
      xmin <- max(group_points_extent@xmin - extent_buffer, output_model_extent@xmin)
      ymin <- max(group_points_extent@ymin - extent_buffer, output_model_extent@ymin)
      xmax <- min(group_points_extent@xmax + extent_buffer, output_model_extent@xmax)
      ymax <- min(group_points_extent@ymax + extent_buffer, output_model_extent@ymax)
      thisGroupExtent <- extent(xmin, xmax, ymin, ymax)

      ### generate a weight grid for each lineage

      if (distance_method == "model-cost") {
        cat("\n\nCreating a transition matrix based on the SDM", "\n")

#        model_cost.ras  <- -1 * log(SDM.ras)     # this is the original version of model cost  
        
      }

      cat ("\nLooping through the lineages ", "to generate weight grids\n")
      count <- 0

      for (lineage in lineages) {

        count <- count +1

        thisLineagePoints <- allLineagePoints[allLineagePoints@data[, lineage_field_name] == lineage, ]
        # create a cost distance layer for the current lineage
        if (lineage == "0") {
          cat ("\nCreating distance layer for sequenced locations of unknown lineage")
        } else {
          cat ("\nCreating distance layer for lineage", lineage)
        }

        if (distance_method == "model-cost") {                                   ## STEP 5b
          # calculates the least cost distance to the nearest lineage point
          # the result is written directly to lineage_dist_gridname
        
        mif_path <- paste(target_dir, lineage, ".mif", sep="")
        writeOGR(thisLineagePoints, mif_path, lineage, driver="MapInfo File", dataset_options="FORMAT=MIF", check_exists=FALSE, overwrite_layer=TRUE)
        cost_out <- paste(target_dir, lineage, ".sdat", sep="")

        rsaga.geoprocessor("grid_analysis",0,list(DEST_TYPE ="0", DEST_POINTS=mif_path, COST=cost_in, ACCUMULATED=cost_out, THRESHOLD=0.000000), env=aa)
        lin_dist.ras <- raster(cost_out)
        projection(lin_dist.ras) <- CRS(spRef)


        } else {
          thisLineagePoints.ras <- rasterize(thisLineagePoints, SDM.ras, 1)
          lin_dist.ras <- distance(thisLineagePoints.ras) / 1000     ## STEP 5a
        }

        # change zero values to a very small non-zero value, to avoid nodata in division
        lin_dist.ras[lin_dist.ras < min_SDM_value] <- min_SDM_value

        if (weight_function == "inverse_square") {                 ## STEP 6b
          lin_weight.ras <- 1 / (lin_dist.ras ^ 2)
        } else if (weight_function == "inverse_cube") {
          lin_weight.ras <- 1 / (lin_dist.ras ^ 3)
        } else if (weight_function == "inverse_quad") {
          lin_weight.ras <- 1/(lin_dist.ras ^ 4)
        } else if (weight_function == "inverse") {
          lin_weight.ras <- 1/lin_dist.ras
        }

        # add the results to stack of distance layers for this group
        if (count == 1) {
          weight.stack <- stack(lin_weight.ras)
        } else {
          weight.stack <- stack(weight.stack, lin_weight.ras)
        }
        names(weight.stack)[count] <- lineage
      }

      if (scale_to == "model") {
        weight_sum.ras <- sum(weight.stack)
        model.stack <- weight.stack / weight_sum.ras
        model.stack <- model.stack * SDM.ras
      } else if (scale_to == "one") {
        weight_sum.ras <- sum(weight.stack)
        model.stack <- weight.stack / weight_sum.ras
      } else if (scale_to == "none") {
        model.stack <- weight.stack
      }

      names(model.stack) <- names(weight.stack)

    } else {
      SDM.ras[SDM.ras < min_SDM_value] <- min_SDM_value
      model.stack <- stack(SDM.ras)
      names(model.stack)[1] <- lineages[1]
    }

    # write results to file for this group
    model_names <- names(model.stack)

    cat("\n")

    for (i in 1:nlayers(model.stack)) {
      model.ras <- model.stack[[i]]
      model_filename <- paste("lin_model_", model_names[i], ".asc", sep="")
      model_path <- paste(target_dir, model_filename, sep="")
      cat("writing model ascii for:", model_filename, "\n")
      writeRaster(model.ras, model_path, overwrite=T, datatype='FLT4S', NAflag=-9999, prj=T)
    }
	model.max <- max(model.stack)
	model_path_max <- paste(target_dir, "max_prob.asc", sep="")
  writeRaster(model.max, model_path_max, overwrite=T, datatype='FLT4S', NAflag=-9999, prj=T)

  model.sum <- sum(model.stack)
  model_path_sum <- paste(target_dir, "sum.asc", sep="")
  writeRaster(model.sum, model_path_sum, overwrite=T, datatype='FLT4S', NAflag=-9999, prj=T)
  
    cat ("\nLineage models done", "\n")
    time_diff <- difftime(Sys.time(), time_start, units='mins')
    cat ("\nTime elapsed:", round(time_diff,2), "minutes\n")

