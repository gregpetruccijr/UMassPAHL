#' Re-integreates ActiGraph data into whatever epoch you need, producing the same results that ActiLife does. 
#'
#' \code{read.ag} Will return ActiGrpah data integreated into the epoch that you desire.
#'
#' @param ag_data_1_sec The ActiGraph data in the original epoch you have.
#' @param epoch The new epoch length that you want your data to be in.

#' @return Returns the AcrtiGraph data in the new epoch.
#'
#' @examples
#' actigraph_data_epoch = ag_epochr(ag_data_1sec = ag_data_1sec_hip, epoch = 60)

## Re-integrates 1-second count data from ActiGraph to larger Epoch----
ag_epochr = function(ag_data_1sec,epoch = 60){
  rows = nrow(ag_data_1sec)
  
  new_rows = floor(rows/epoch)
  
  new_rows_index = (seq(1:new_rows)*epoch)+1
  new_rows_index = c(1,new_rows_index)
  new_rows_index = new_rows_index[-length(new_rows_index)]
  
  Date = ag_data_1sec$Date[new_rows_index]
  Time = ag_data_1sec$Time[new_rows_index]
  count_data = data.frame(Axis1 = rep(0,new_rows), Axis2 = rep(0,new_rows), Axis3 = rep(0,new_rows))
  
  for (i in 1:new_rows){
    
    if(i+1 <= new_rows){
      count_data$Axis1[i] = sum(ag_data_1sec$Axis1[new_rows_index[i]:(new_rows_index[i+1]-1)])
      count_data$Axis2[i] = sum(ag_data_1sec$Axis2[new_rows_index[i]:(new_rows_index[i+1]-1)])
      count_data$Axis3[i] = sum(ag_data_1sec$Axis3[new_rows_index[i]:(new_rows_index[i+1]-1)])
      
    } else {
      count_data$Axis1[i] = sum(ag_data_1sec$Axis1[new_rows_index[i]:(new_rows_index[i]+epoch-1)])
      count_data$Axis2[i] = sum(ag_data_1sec$Axis2[new_rows_index[i]:(new_rows_index[i]+epoch-1)])
      count_data$Axis3[i] = sum(ag_data_1sec$Axis3[new_rows_index[i]:(new_rows_index[i]+epoch-1)])
      
      # Handles excess if not exact 1 hour length
      # if((new_rows_index[i]+epoch-1) < rows){
      #   excess$Axis1 = sum(ag_data_1sec$Axis1[(new_rows_index[i]+epoch-1):rows])
      #   excess$Axis2 = sum(ag_data_1sec$Axis2[(new_rows_index[i]+epoch-1):rows])
      #   excess$Axis3 = sum(ag_data_1sec$Axis3[(new_rows_index[i]+epoch-1):rows])
      # }
    }
    
  }
  
  epoch_data = data.frame(Date = Date, Time = Time, count_data)
  epoch_data$VM = sqrt(epoch_data$Axis1^2 + epoch_data$Axis2^2 + epoch_data$Axis3^2)
  
  return(epoch_data)
}
