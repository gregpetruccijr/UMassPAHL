#' Reads in the Noldus Observer XT output
#'
#' \code{read.noldus} Will return the Noldus excel file exported from Observer XT
#'
#' @param filepath The filepath where the noldus file exists

#' @return Returns the Noldus excel file.
#'
#' @examples
#' noldus_data = read.noldus(filepath)
#'

read.noldus = function(filepath){
  if(str_detect(filepath,'.csv')){
    file_data = read.csv(filepath, header = T)
  } else {
    file_data = read_xlsx(filepath, col_names = T)
    file_data$Modifier_1 = as.numeric(file_data$Modifier_1)
  }
}
