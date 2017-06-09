#' Analyse geomagnetic data from archaeological context
#'
#' @param anomalies_sdf    Spatialdataframe: Point data containing the Anomalies
#' @param magnetic_raster    RasterLayer: Geomagnetic data (nT)
#' @param searchradius    The Radius around every anomaly that is under observation
#' @param dipolfactor    1. version of defining a dipole: The factor used to descripe a dipole; A dipol is defined by |minimal value * dipolfactor| > maximal value
#' @param dipol_minima    Additional defining of a dipole: All anomalies with min nT < dipol_minima in the profile are defined as dipoles
#' @param angle_steps    Profiles are generated around the anomalie points in steps of angle_steps degrees
#' @param cut_value     nT-value at which the width of the amplitude will be determined
#' @param method    method for calculating the profile values  avg = Average width of all profiles; median = Median width of all profiles
#' @param get_dipol    Identifing the dipoles (TRUE)
#' @param get_profile_values    Calculating the values of all profiles (TRUE)
#' @return A spatial data frame with the results of the dipole identifing (0/1) or/and the width and heigth of the amplitudes
#' @examples
#' sum(1:10)
analyseMagnetic <- function(anomalies_sdf,magnetic_raster,searchradius=2,dipolfactor=2,angle_steps=10,get_profile_values=FALSE,get_dipol=FALSE,method="avg",cut_value=1.5,dipol_minima=-4){
  ########Start

  #searcradius was diameter first
  searchradius <- searchradius * 2
  #empty column for result of dipol test
  if(is.null(anomalies_sdf@data$dipol_kB) && get_dipol == TRUE){
    anomalies_sdf@data$di_kB <- 0
  }
  #empty column for profile values
  if (is.null(anomalies_sdf@data$an_breite) && get_profile_values == TRUE)
  {
    anomalies_sdf@data$an_breite <- 0
    anomalies_sdf@data$an_hoehe <- 0
  }
  #For each point of an anomaly profiles are generated to analyse the values.
  #The amount of profiles is angle_steps / 180, because for every angle between 0 and 180 in angle_steps steps
  #a profile will be generated
  #For the dipol the test will compare the highest and lowest values of all profile of an anomalie
  #For getting the profile values there will be two methods; the profile values are the width of each anomalie
  #the width is defined by the distance between at the cutting value
  #avg will calculate the average distance of all profiles
  #median will calculate the median value of all profiles

  #total amount anomalies
  n_anomalies <- length(anomalies_sdf)
  #resolution of inputraster = the distance between the points in the profile
  x_resolution <- xres(magnetic_raster)


  #The tests has to made for every point in anomalies
  k <- 1
  while (k <= n_anomalies){
    #Because of the relativly long time of calculating it counts backwarts and shows every 10th value
    if (k %% 10 == 0) {
      print(n_anomalies - k)
    }

    #If median was choosen, declare the a vector to buffer the values
    if (method == "median"){
      median_vector = numeric()
    }

    #coordinaten of the first anomaly
    coord <- data.frame(anomalies_sdf@coords[k,1], anomalies_sdf@coords[k,2])
    colnames(coord) <- c("x", "y")
    rownames(coord) <- NULL

    #The amount of points that are needed in a profile
    n_newpoints <- round(searchradius / x_resolution)
    #if the value is even, add one. therefore the profile is centric
    if(n_newpoints %% 2 == 0){
      n_newpoints <- n_newpoints + 1
    }


    #For each point a coordinate should be added
    #Starting at the most left point of the profile
    startx <- coord$x - (searchradius / 2)
    #Generating the points
    i <- 2
    increase <- 0
    while (i <= n_newpoints + 2){
      coord[i,] <- c(startx + increase, coord$y[1])
      increase <- increase + x_resolution
      i <- i + 1
    }
    #The result is a profile from east to west throught the anomaly
    #The test should be done for every profile in angle_steps
    #Starting with the degree 0
    dregee <- 0
    #For every angle between 0 and 180
    while (dregee < 180)
    {
      #if it is the first run, no recalculating of the profile is needed
      if (dregee == 0){
        #coord_trans is needed for the calculations
        coord_trans <- coord

      } else {
        #If the degree is > 0, the coodinates have to be rotated
        #Center is the middle point
        #This is the first value of the data frame
        #For the calculating the new coordinates use translation and rotation
        x_mittel <- coord$x[1]
        y_mittel <- coord$y[1]
        coord_trans <- coord
        for (z in 1:nrow(coord)){
          coord_trans[z,] <- c(
            x_mittel + (coord$x[z] - x_mittel) * cos(dregee / 180 * pi) - sin(dregee / 180 * pi) * (coord$y[z] - y_mittel),
            y_mittel + (coord$x[z] - x_mittel) * sin(dregee / 180 * pi) + (coord$y[z] - y_mittel) * cos(dregee / 180 * pi))
          #http://www.matheboard.de/archive/460078/thread.html

        }

      }

      #Getting the value of the magnetic raster at everypoint of the profile
      mag_profile <- data.frame(coordinates(coord_trans), extract(magnetic_raster, coord_trans))
      names(mag_profile) <- c("x", "y", "value")
      #First value can be deleted
      mag_profile <- mag_profile[-1, ]
      #A dipol is defined by |minimal value * dipolfactor| > maximal value
      #na.omit, they can be generated be boundary effects
      if (get_dipol == TRUE && abs(min(na.omit(mag_profile$value))) * dipolfactor > max(na.omit(mag_profile$value)) &&
          min(na.omit(mag_profile$value)) < 0)
      {
        #writing the result in spatial data frame
        anomalies_sdf@data$di_kB[k] <- 1
      }
      #The dipol can also be defined as a minimum < dipol_minima
      if (get_dipol == TRUE && min(na.omit(mag_profile$value)) < dipol_minima){
        #writing the result in spatial data frame
        anomalies_sdf@data$di_kB[k] <- 2
      }

      #If choosen, get width and height values of each profile
      if (get_profile_values == TRUE){
        mag_profile <- na.omit(mag_profile)
        #Getting the highest value (mor eor less the middle of the profile) and the starting point
        x_row <- which(mag_profile$value == max(na.omit(mag_profile$value)))[1]
        #result buffer
        distance_right <- 0
        distance_left <- 0
        #Going step by step from the middle point to the right until the nT value is lower the cutting value
        w <- x_row
        #Stopp if value has reached
        stopp <- FALSE

        while(w <= nrow(mag_profile) && stopp == FALSE){
          #HIER PRÜFEN OB DER ERSTE WERT UNTER DEM CUT VALUE LIEGT, WENN JA DANN DISTANCE = 0
          #If the cutting value has reached, claculate the nT value at the cutting value
          if(mag_profile$value[w] < cut_value){
            #getting the distance between the highest point and the value before the cutting value had reached
            x1 <- mag_profile$x[w-1]
            y1 <- mag_profile$y[w-1]
            x2 <- mag_profile$x[w]
            y2 <- mag_profile$y[w]
            #getting the slope
            m <- (y2-y1) / (x2-x1)
            #calculating the distance between the two points
            distance <- sqrt(  ((x2 - x1) ^2) + ((y2 - y1) ^2))
            #generate dataframe with values
            #xvalue is the distance between the two points
            xw <- c(0,distance)
            #yvalue is the nT magnetic value at this points
            yw <- c(mag_profile$value[w-1], mag_profile$value[w])
            #to get the magnetic value at the cutting point, a linear regression is used
            fm <- lm(xw ~ yw)
            #predict the xvalue(distance) for the yvalue(cutting value)
            a <- predict(fm, data.frame(yw = c(cut_value)), se.fit = TRUE)$fit
            #The distance from the middle point to the cutting value is the distance to the point before
            #reaching the cutting value and the predicted distance
            distance_right <- (w - x_row - 1) * distance + a
            if(distance_right < 0 || w == x_row)
            {distance_right <- 0}
            stopp <- TRUE
          }
          w <- w + 1
        }
        #Same will be done for the left site of the graph
        w <- x_row
        stopp <- FALSE
        while(w > 0 && stopp == FALSE){
          if(mag_profile$value[w] < cut_value){
            #HIER PRÜFEN OB DER ERSTE WERT UNTER DEM CUT VALUE LIEGT, WENN JA DANN DISTANCE = 0
            x1 <- mag_profile$x[w+1]
            y1 <- mag_profile$y[w+1]
            x2 <- mag_profile$x[w]
            y2 <- mag_profile$y[w]
            m <- (y2 - y1) / (x2 - x1)
            stopp <- TRUE
            distance <- sqrt(((x2 - x1) ^2) + ((y2 - y1) ^2))
            xw <- c(0, distance)
            yw <- c(mag_profile$value[w], mag_profile$value[w+1])
            fm <- lm(xw ~ yw)
            a <- predict(fm, data.frame(yw = c(cut_value)), se.fit = TRUE)$fit
            distance_left <- (x_row - w) * distance - a
            if(distance_left<0 || w == x_row)
            {distance_left <- 0}
          }
          w <- w - 1
        }

        #if the cutting value had been reached in all cases  distance_left>0 &&distance_right>0
        #if the method "avg" was choosen, use the mean value of all profiles
        if (distance_left > 0 && distance_right > 0 && method == "avg"){
          if (anomalies_sdf@data$an_breite[k] == 0){
            anomalies_sdf@data$an_breite[k] <- distance_left + distance_right
          } else if (anomalies_sdf@data$an_breite[k] > 0){
            anomalies_sdf@data$an_breite[k] <- (anomalies_sdf@data$an_breite[k] + distance_left + distance_right) / 2
          }
        }

        #If the method is median, write the result to the buffer
        if (distance_left > 0 && distance_right > 0 && method == "median" && dregee == 0){
          median_vector <- c(distance_right+distance_left)
        } else if (distance_left > 0 && distance_right > 0 && method == "median" && dregee > 0){
          median_vector <- c(median_vector, distance_right + distance_left)
        }

        #The heigth (magnetic value) is alle the time the highest value
        if (max(na.omit(mag_profile$value))-cut_value > 0){
         anomalies_sdf@data$an_hoehe[k] <- max(na.omit(mag_profile$value))-cut_value
        }
      }

      #increasing the angle
      dregee <- dregee + angle_steps
    }
    #for median, after calculating all values, get the median
    if (method == "median"){
      anomalies_sdf@data$an_breite[k] <- median(median_vector)
    }

    k <- k + 1
  }


  #return value
  return(anomalies_sdf)

}




