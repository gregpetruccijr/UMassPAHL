#' Expands the Noldus data created by 'read.noldus'
#'
#' \code{trim.noldus} Will return the Noldus data expanded to 0.01 second data for aligning with accelerometer data
#'
#' @param noldus_data The Noldus data imported using read.noldus
#'
#' @param duration The duration of observations of interest in seconds
#'
#' @return Returns the expanded Noldus data.
#'
#' @examples
#' session_noldus = trim.noldus(noldus_data)
#'

trim.noldus = function(noldus_data, duration){
  noldus_data <- noldus_data[noldus_data$Event_Type=="State start",] # remove all stop times since they're 0 duration
  secs <- round(noldus_data$Duration_sf*100) # round duration of obs to nearest hundredth of a second (hence * 100)

  big.DO.1 <- data.frame(Time_Relative_sf = rep(noldus_data$Time_Relative_sf,times = secs),
                         Duration_sf = rep(noldus_data$Duration_sf,times = secs),
                         Behavior=rep(noldus_data$Behavior,times=secs),
                         METs=rep(noldus_data$Modifier_1,times=secs),
                         Modifier_2=rep(noldus_data$Modifier_2,times=secs),
                         Modifier_3 = rep(noldus_data$Modifier_3,times = secs),
                         Modifier_4 = rep(noldus_data$Modifier_4,times = secs)) #creates DF that only has Beh, METs, Mod
  big.DO.1$Behavior <- as.character(big.DO.1$Behavior)  # transforms Beh and Mod to character strings
  big.DO.1$Modifier_2 <- as.character(big.DO.1$Modifier_2)
  big.DO.1$Modifier_3 <- as.character(big.DO.1$Modifier_3)
  big.DO.1$Modifier_4 <- as.character(big.DO.1$Modifier_4)
  #big.DO.1$Modifier_5 <- sapply(big.DO.1$Behavior, function(x) ifelse(any(str_detect(x,c('Lying','Sitting'))),'Sedentary','Active'))

  if(nrow(big.DO.1) > (duration *100)){
    excess = nrow(big.DO.1)-(duration*100)
    big.DO.1 = big.DO.1[1:(nrow(big.DO.1)-excess),]
  }

  # Handles files that are not a full 1 hour of data
  if(nrow(big.DO.1) < duration){

    # Rounds down to nearest second, removes excess milliseconds of observations
    # Round down because can't create observations that didn't exist
    excess = floor(nrow(big.DO.1)/100)
    big.DO.1 = big.DO.1[1:(excess*100),]
  }
  return(big.DO.1)
}

