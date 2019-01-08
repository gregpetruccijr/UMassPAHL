#' Combines direct observation data with trimmed one-second count data from the accelerometer
#'
#' \code{DO.info} Extracts the behavior and activity (Modifier 2) info from trim.noldus and appends observations to accelerometer output
#'
#' @param ag_data_1sec The one-second count accelerometer data exported from trim.ag
#'
#' @param noldus_expanded The Noldus data created from trim.noldus
#'
#' @return Returns the original 1-second accelerometer data with appended behaviors and activities.
#'
#' @examples
#' session_data = DO.info(Session_AG_1sec_hipdata, noldus_expanded)
#'

DO.info <- function(ag_data_1sec, noldus_expanded) {
  # Pre-dispose columns to be filled
  ag_data_1sec$Behavior <- NA
  ag_data_1sec$Behavior.2 = NA
  ag_data_1sec$Modifier2 <- NA
  ag_data_1sec$Modifier2.2 <- NA


  if ((nrow(noldus_expanded) / 100) > nrow(ag_data_1sec)) {
    excess.amount <-
      (((nrow(
        noldus_expanded
      ) / 100) - nrow(ag_data_1sec)) * 100)
    noldus_expanded <-
      noldus_expanded[1:(nrow(noldus_expanded) - excess.amount), ]
  } else {
    NULL
  }

  milli <- dim(noldus_expanded)[1]
  secs <- floor(milli / 100)

  for (i in 1:secs) {
    if (i == 1) {
      temp.inds <- noldus_expanded[i:100,]

      # Find behaviors in first 1-second window, if a transition behavior.2 is filled
      temp.behav <- as.data.frame(table(temp.inds$Behavior))

      if (nrow(temp.behav) >= 2) {
        temp.behav = arrange(temp.behav, desc(Freq))
        temp.behav1 <- paste(temp.behav[1, 1])
        temp.behav2 = paste(temp.behav[2, 1])
        ag_data_1sec$Behavior.2[i] = paste(as.character(temp.behav.2))

      } else{
        temp.behav.mode <- max(temp.behav$Freq)
        temp.behav.inds <- which(temp.behav$Freq == temp.behav.mode)
        temp.behav1 <- temp.behav$Var1[temp.behav.inds]
      }

      # Find activities (modifier2) in first 1-second window, if a transition modifier2.2 is filled
      temp.mod2 <- as.data.frame(table(temp.inds$Modifier_2))

      if (nrow(temp.mod2) >= 2) {
        temp.mod2 = arrange(temp.mod2, desc(Freq))
        temp.mod2.1 <- paste(temp.mod2[1, 1])
        temp.mod2.2 = paste(temp.mod2[2, 1])
        ag_data_1sec$Modifier2.2[i] <-
          paste(as.character(temp.mod2.2))

      } else{
        temp.mode <- max(temp.mod2$Freq)
        temp.mod2.inds <- which(temp.mod2$Freq == temp.mode)
        temp.mod2.1 <- temp.mod2$Var1[temp.mod2.inds]

      }

      ag_data_1sec$Behavior[i] <- paste(as.character(temp.behav1))
      ag_data_1sec$Modifier2[i] <- paste(as.character(temp.mod2.1))

    } else {

      temp.inds <- noldus_expanded[(((i - 1) * 100) + 1):(i * 100),]

      # Find behaviors in 1-second window, if a transition behavior.2 is filled
      temp.behav <- as.data.frame(table(temp.inds$Behavior))

      if (nrow(temp.behav) >= 2) {
        temp.behav = arrange(temp.behav, desc(Freq))
        temp.behav1 <- paste(temp.behav[1, 1])
        temp.behav2 = paste(temp.behav[2, 1])
        ag_data_1sec$Behavior.2[i] = paste(as.character(temp.behav2))

      } else{
        temp.behav.mode <- max(temp.behav$Freq)
        temp.behav.inds <- which(temp.behav$Freq == temp.behav.mode)
        temp.behav1 <- temp.behav$Var1[temp.behav.inds]
      }

      # Find activities (modifier2) in first 1-second window, if a transition modifier2.2 is filled
      temp.mod2 <- as.data.frame(table(temp.inds$Modifier_2))

      if (nrow(temp.mod2) >= 2) {
        temp.mod2 = arrange(temp.mod2, desc(Freq))
        temp.mod2.1 <- paste(temp.mod2[1, 1])
        temp.mod2.2 = paste(temp.mod2[2, 1])
        ag_data_1sec$Modifier2.2[i] <-
          paste(as.character(temp.mod2.2))

      } else{
        temp.mode <- max(temp.mod2$Freq)
        temp.mod2.inds <- which(temp.mod2$Freq == temp.mode)
        temp.mod2.1 <- temp.mod2$Var1[temp.mod2.inds]

      }

      ag_data_1sec$Behavior[i] <- paste(as.character(temp.behav1))
      ag_data_1sec$Modifier2[i] <- paste(as.character(temp.mod2.1))

    }
  }

  if(any(is.na(ag_data_1sec$Behavior))){
    na_indices = which(is.na(ag_data_1sec$Behavior))
    ag_data_1sec = ag_data_1sec[-na_indices,]
  }

  return(ag_data_1sec)
}
