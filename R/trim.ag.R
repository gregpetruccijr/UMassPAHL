#' Trims ActiGraph data to give you data during specific intervals that you want to analyize
#'
#' \code{trim.ag} Will return ActiGrpah data integreated into the epoch that you desire.
#'
#' @param ag_data The data that results from running read.ag
#' @param vid_times A spreadhseet that has the video start and stop times for each participant
#' @param time_index The index start and stop times of the particpant you want to analyize, from the origrinal spreadsheet
#' @param raw Indicate whether the ag_data file is raw or count data
#'
#' @return Returns the AcrtiGraph data trimmed to the specific time that you want to analyize.
#'
#' @examples
#' actigraph_data_epoch = ag_epochr(ag_data_1sec = ag_data_1sec_hip, epoch = 60)

## Trimmed ActiGraph data to include only data between sesion start and session end indices----
trim.ag = function(ag_data,vid_times,time_index, raw = F){
  real_date = mdy(vid_times$Date[time_index])
  real_date = as.character(real_date)
  real_start = vid_times$Session.Start.Time[time_index]
  real_end = vid_times$Session.Stop.Time[time_index]

  session_ag_data = filter(ag_data, Date == real_date)
  real_start_index = which(str_detect(session_ag_data$Time,real_start))
  real_end_index = which(str_detect(session_ag_data$Time,real_end))

  if(raw == F){
    session_ag_data = session_ag_data[real_start_index:real_end_index,]

    if(nrow(session_ag_data) > 3600){
      excess = nrow(session_ag_data)-3600
      session_ag_data = session_ag_data[1:(nrow(session_ag_data)-excess),]
    }

    # Handles sessions that were not collected for an hour
    if(nrow(session_ag_data) < 3600){
      excess = floor(nrow(session_ag_data)/10)
      session_ag_data = session_ag_data[1:(excess*10),]
    }
  } else {
    session_ag_data = session_ag_data[real_start_index[1]:(real_end_index[1]-1),]

    # Handles 80 Hz raw data, need to change multiplier if other sampling frequency
    if(nrow(session_ag_data) > 3600*80){
      excess = nrow(session_ag_data)-3600*80
      session_ag_data = session_ag_data[1:(nrow(session_ag_data)-excess),]
    }

  }

  return(session_ag_data)
}
