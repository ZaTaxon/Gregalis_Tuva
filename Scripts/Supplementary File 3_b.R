library(raster)
library(sp)
library(stringr)
library(gdistance)
library(RSAGA)
library(rgdal)

############### PARAMETERS ###############
#Directories and files settings
base_dir              <- "G:/gregalis/"
target_dir            <- paste(base_dir, "distances/", sep="")  # where the lineage model grids and working data go
SDM_model_base        <- paste(base_dir, "sdm/", sep="")
cost_in               <- paste(SDM_model_base, "inverse_gregalis.asc", sep="")
lineage_site_filename <- paste(base_dir, "tabs/", "clades_single_new.csv", sep="")
lineage_field_name    <- "Clade"   # the column for individual name for each point

aa <- rsaga.set.env(workspace="C://Program Files/saga-7.2.0_x64/", cmd="saga_cmd.exe", path="C://Program Files/saga-7.2.0_x64/", modules="C://Program Files/saga-7.2.0_x64//tools", version="7.2.0", parallel=FALSE)
aa$version = rsaga.get.version(env=aa)

#Various filter settings
lin_exclude_list      <- c("A82", "A81", "A43", "A42", "A41", "A11", "A1")     # this list allows for skipping at the individual level. Delimiter comma, IDs in double quotes, as c("A82", "A81", "A43")
use_use               <- 0       # use column "use" for records filter ("1") or not ("0")
additional_buffer     <- 0       # how much (as a proportion) the output grids should extend beyond the buffered points
buffer_dist           <- 5000    # the buffer distance in map units (presumably decimal degrees)

#Grid settings
grid_resolution       <- 2000
spRef                 <- "+proj=moll +lon_0=30 +x_0=3335846.22854 +y_0=-336410.83237 +datum=WGS84 +units=m +no_defs"
output_model_extent   <- extent(3797130.398597037, 11593130.39859704, 3685333.673568079, 7773333.673568079)    # the maximum extent for all lineage models

#Points settings
ppRef                 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
############### END OF PARAMETERS ###############


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

  all_lineages_sites[, lineage_field_name] <- str_trim(all_lineages_sites[, lineage_field_name])

  groupLineages     <- unique(all_lineages_sites[, lineage_field_name])

  cat (orig_rowcount, "rows read,", nrow(all_lineages_sites), "valid records loaded\n")

  # turn lineage points into a SpatialPoints object
  allLineagePoints <- SpatialPointsDataFrame(all_lineages_sites[, c("longitude", "latitude")], data=all_lineages_sites, proj4string = CRS(ppRef))
  allLineagePoints <- spTransform(allLineagePoints, spRef)
  allLineagePoints <- crop(allLineagePoints, output_model_extent)
   
    group_extent <- output_model_extent
    
    res.cost <- array(dim = c(orig_rowcount, orig_rowcount), dimnames = list(groupLineages,groupLineages))
    mif_all_path <- paste(target_dir, "all", ".mif", sep="")
    writeOGR(allLineagePoints, mif_all_path, "lineage", driver="MapInfo File", dataset_options="FORMAT=MIF", check_exists=FALSE, overwrite_layer=TRUE)
    
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

      ### generate a weight grid for each lineage  START OF STEP 4

      cat ("\nLooping through the lineages ", "to generate accumalated cost grids\n")
      count <- 0
      xs <- c(1:orig_rowcount)
      
      for (x in xs) {

        count <- count +1

        thisLineagePoints <- allLineagePoints[x, ]
        lineage <- thisLineagePoints@data[1,1]
#       check intersection of SDM and the points
#        d <- intersect(SDM.ras,thisLineagePoints)

        # create a cost distance layer for the current lineage
        if (lineage == "0") {
          cat ("\nCreating accumulated cost layer for sequenced locations of unknown lineage")
        } else {
          cat ("\nCreating accumulated cost layer for lineage", lineage)
        }

          # calculates the least cost distance to the nearest lineage point
        
        mif_path <- paste(target_dir, x, ".mif", sep="")
        writeOGR(thisLineagePoints, mif_path, lineage, driver="MapInfo File", dataset_options="FORMAT=MIF", check_exists=FALSE, overwrite_layer=TRUE)
        cost_out <- paste(target_dir, x, ".sdat", sep="")

        rsaga.geoprocessor("grid_analysis",0,list(DEST_TYPE ="0", DEST_POINTS=mif_path, COST=cost_in, ACCUMULATED=cost_out, THRESHOLD=0.000000), env=aa)
        lin_dist.ras <- raster(cost_out)
        projection(lin_dist.ras) <- CRS(spRef)

        cat ("\nCreating least cost path lines for lineage", lineage)
        
        dist_mif_path <- paste(target_dir, x, "d.mif", sep="")
        rsaga.geoprocessor("grid_analysis",5,list(SOURCE=mif_all_path, DEM=cost_out, LINE=dist_mif_path), env=aa)  
      
        cost.dist <- readGDAL(cost_out)
        r.cost.dist <- as(cost.dist,"RasterLayer")
        v.dst <- readOGR(mif_all_path)
      
        pos <- cellFromXY(r.cost.dist, v.dst@coords)
        print(v.dst@data$dist <- r.cost.dist[pos])
        res.cost[x,] <- r.cost.dist[pos]
	 }

    } else {
        stop("Nothing to analyse: no multiple points set.")
    }

    write.table(res.cost,file=paste(target_dir, "lcp", ".csv", sep=""), append = FALSE, col.names = NA, sep = ";")
    cat("\n")
    res.euclid <-pointDistance(allLineagePoints@coords, allLineagePoints@coords, type='Euclidean', lonlat = F, allpairs = TRUE)
    res.euclid <- as.array(res.euclid)
    dimnames(res.euclid) <- list(groupLineages,groupLineages)
    write.table(res.euclid,file=paste(target_dir, "euclid", ".csv", sep=""), append = FALSE, col.names = NA, sep = ";")
    cat("Pearson correlation LCP*Euclid: ", cor(c(res.cost), c(res.euclid), method="pearson"))

    cat ("\nLeast cost path calculations done", "\n")
    time_diff <- difftime(Sys.time(), time_start, units='mins')
    cat ("\nTime elapsed:", round(time_diff,2), "minutes\n")
