#' logscale
#' Create logarithmically spaced scales
#' @param scale_min an integer indicating the minimum scale to be resovled
#' @param scale_max an integer indicating the maximum scale to be resolved
#' @param scale_ratio a double indiciting the ratio by which scale successive 
#' scales. For example, scale_ratio = 2 would create a scales increasing by 
#' a power of 2.
#' @export
logscale = function(scale_min, scale_max, scale_ratio){
  number_of_scales =  ceiling(log(scale_max/scale_min)/log(scale_ratio));
  scales = rep(NA, number_of_scales)
  scales[1] = scale_min;
  for (i in 2:number_of_scales){
    scales[i] = scales[i-1]*scale_ratio
  }
  
  return(unique(floor(scales)))
  
}