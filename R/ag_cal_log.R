#' Calibration Finder for specific ActiGraph devices using serial number
#'
#' \code{ag_cal_log} Will return the serial number of the selected file with the recommended calibration settings specific to that device
#'
#' @param filepath The filepath where the noldus file exists

#' @return Returns the recommended calibration offsets and scaling factor that should be applied to raw accelerometer data
#'
#' @examples
#' acc_calibration_settings = ag_cal_log(filepath)
#'

ag_cal_log = function(filepath){

  ag_header = read.csv(filepath,header = F,stringsAsFactors = F, nrows = 10)

  device_serial = str_split(ag_header[2,],'Number: ')[[1]][2]

  C = g.calibrate(filepath, use.temp = F, printsummary=F)

  ag_cal_output = data.frame(device_serial, t(C$offset), t(C$scale))

  colnames(ag_cal_output) = c('Serial','Offset_X','Offset_Y','Offset_Z','Scale_X','Scale_Y','Scale_Z')

  return(ag_cal_output)
}
