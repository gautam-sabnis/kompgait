#' Filter data for per line level analyses 
#' 
#' Filters the strains by setting the minimum number of animals per
#' strain
#' @param data_stride a csv file containing data per stride
#' @param data_animal a csv file containing data per animal 
#' @return a list containing filtered data_per_stride and data_per_animal 
#' @examples 
#' filter_data(data_stride = data_per_stride, data_animal = 
#' data_per_animal, threshold = 5)
#' @export 

filter_data <- function(data_stride, data_animal_mean, data_animal_disp, threshold){

	data_per_stride <- data_stride
	data_per_animal_mean <- data_animal_mean
	data_per_animal_disp <- data_animal_disp

	Strains_thresh <- names(table(data_per_animal_mean$Strain))[table(data_per_animal_mean$Strain) >= threshold]
	data_per_animal_mean <- data_per_animal_mean[data_per_animal_mean$Strain %in% Strains_thresh,] 
	data_per_animal_disp <- data_per_animal_disp[data_per_animal_disp$Strain %in% Strains_thresh,] 
	data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains_thresh,]
	data_per_animal_mean$Strain <- droplevels(data_per_animal_mean$Strain)
	data_per_animal_disp$Strain <- droplevels(data_per_animal_disp$Strain)
	data_per_stride$Strain <- droplevels(data_per_stride$Strain)

	return(list(data_per_stride, data_per_animal_mean, data_per_animal_disp))

}