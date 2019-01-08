#' Reads in the Noldus Observer XT output
#'
#' \code{read.ag} Will return ActiGrpah data; can handle count or raw data.
#'
#' @param filepath The filepath where the ActiGraph file exists. 
#' @param calibration_file  The filepath for calibration data.

#' @return Returns the AcrtiGraph data.
#'
#' @examples
#' actigraph_data = read.ag(filepath)

## Read in Raw or Count ActiGraph data, export to environment with Timestamps and VM calculated----
read.ag = function(filepath, calibration_file){
  
  # Accounts for exported data that includes the header or not
  check_data = read.csv(filepath,header = F,skip = 10)[1,]
  if(is.numeric(check_data[1,])){
    file_data = fread(filepath,header = F, skip = 10, stringsAsFactors = F)
  } else {
    file_data = fread(filepath,header = T, skip = 10, stringsAsFactors = F)
  }
  
  ag_header = read.csv(filepath,header = F,stringsAsFactors = F, nrows = 10)
  device_serial = str_split(ag_header[2,],'Number: ')[[1]][2]
  start = str_split(ag_header[3,],'Time ')[[1]][2]
  date = str_split(ag_header[4,],'Date ')[[1]][2]
  frequency = as.numeric(str_split(str_split(ag_header[1,],'at ')[[1]][3],' Hz')[[1]][1])
  
  epoch_full = str_split(ag_header[5,],' 00:')[[1]][2]
  epoch_temp = str_split(epoch_full,':')
  epoch = (as.numeric(epoch_temp[[1]][1])*60) + (as.numeric(epoch_temp[[1]][2]))
  
  file_length = nrow(file_data)
  date_time_start = mdy_hms(paste(date,start, sep = ' '))
  
  if(epoch < 1){
    # For Raw data
    date_time_end = date_time_start + (file_length/frequency)
    Timestamp = seq(from = date_time_start,to = date_time_end, by = 1/frequency)
    Timestamp = Timestamp[1:length(Timestamp)-1]
    C = g.calibrate(filepath,use.temp = F, printsummary=F)
    
    if(C$offset[1] == 0 & C$scale[1] == 1){
      
      device_serial_index = str_which(calibration_file$Serial,device_serial)
      
      if(length(device_serial_index) == 0){
        
      } else {
        C$offset[1] = calibration_file$Offset_X[device_serial_index]
        C$offset[2] = calibration_file$Offset_Y[device_serial_index]
        C$offset[3] = calibration_file$Offset_Z[device_serial_index]
        
        C$scale[1] = calibration_file$Scale_X[device_serial_index]
        C$scale[2] = calibration_file$Scale_Y[device_serial_index]
        C$scale[3] = calibration_file$Scale_Z[device_serial_index]
        
      }
    }
    file_data = mutate(file_data,
                       calX = `Accelerometer X`*C$scale[1] + C$offset[1],
                       calY = `Accelerometer Y`*C$scale[2] + C$offset[2],
                       calZ = `Accelerometer Z`*C$scale[3] + C$offset[3],
                       VM = sqrt(`Accelerometer X`^2 + `Accelerometer Y`^2 + `Accelerometer Z`^2),
                       VMcorrG = sqrt(abs(`Accelerometer X`^2 + `Accelerometer Y`^2 + `Accelerometer Z`^2-1)),
                       ENMO = sqrt(calX^2 + calY^2 + calZ^2)-1)
    
    file_data = mutate(file_data, ENMO = ifelse(ENMO < 0,0,ENMO))
    
  } else {
    # For Count data
    date_time_end = date_time_start + (file_length)
    Timestamp = seq(from = date_time_start, to = date_time_end, by = epoch)
    Timestamp = Timestamp[1:length(Timestamp)-1]
    file_data = mutate(file_data, VM = sqrt(Axis1^2 + Axis2^2 + Axis3^2))
  }
  
  Dates = as.character(date(Timestamp))
  Hour = as.character(hour(Timestamp)) %>% str_pad(width = 2, side = 'left',pad = '0')
  Minute = as.character(minute(Timestamp)) %>% str_pad(width = 2, side = 'left',pad = '0')
  Second = as.character(second(Timestamp)) %>% str_pad(width = 2, side = 'left',pad = '0')
  
  ag_data = cbind(Dates,paste(Hour,Minute,Second,sep = ':'),file_data, stringsAsFactors = F)
  
  if(epoch<1){
    colnames(ag_data) = c('Date','Time','AxisX','AxisY','AxisZ', 'CalibratedX','CalibratedY','CalibratedZ','VM', 'VMcorrG', 'ENMO')
  }else{
    colnames(ag_data) = c('Date','Time','Axis1','Axis2','Axis3','VM')
  }
  
  return(ag_data)
}

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
