#' Select control animals from control strain 
#' 
#' Choose control animals from the control strain to be examined during a specific time frame (in
#' days)
#' @param CtrlStrain specify the control strain
#' @param tw specify the time period in days 
#' @return a data.frame containing containing Strain, MouseID, Type 
#' @examples 
#' komp_select_controls(CtrlStrain="C57BL/6NJ", tw = 21)
#' @export 

komp_select_controls <- function(CtrlStrain, tw){
	control.list <- subset(data_per_stride, Strain == CtrlStrain)[,c('MouseID','Strain','TestDate')]
	mutant.list <- subset(data_per_stride, Strain != CtrlStrain)[,c('MouseID','Strain','TestDate')]
	mutant.list$Type <- 'Mutant'
	all.control.dates.df <- NULL
	control.dates.df.tmp <- NULL
	all.control.mouseid.df <- NULL

	for (s in unique(mutant.list$Strain)){
    	mutant.dates.list <- unique(subset(mutant.list,Strain == s)$TestDate)
    	control.dates.list <- NULL
    	mouse.id.list <- NULL
    	for (d in mutant.dates.list){
        	start.d <- as.Date(d) - tw
        	end.d <- as.Date(d) + tw
        	control.dates <- unique(subset(control.list, as.Date(TestDate) >= start.d & as.Date(TestDate) <= end.d)$TestDate)
        	control.dates.list <- c(control.dates.list, format(control.dates, format = '%Y-%m-%d'))
        	mouse.id <- unique(subset(control.list, as.Date(TestDate) >= start.d & as.Date(TestDate) <= end.d)$MouseID)
    		mouse.id.list <- c(mouse.id.list, as.character(mouse.id))
    	}
    	if (identical(mouse.id.list,character(0))) next
    	control.dates.df <- data.frame(Strain = s, TestDate = format(control.dates.list, format = '%Y-%m-%d'), Type = 'Control')
    	all.control.dates.df <- rbind(all.control.dates.df, control.dates.df[,c('Strain','Type')])
    	control.dates.df.tmp <- c(control.dates.df.tmp, format(control.dates.list, format = '%Y-%m-%d'))

    	control.mouseid.df <- data.frame(Strain = s, MouseID = mouse.id.list, Type = 'Control')
    	all.control.mouseid.df <- rbind(all.control.mouseid.df,control.mouseid.df[,c('Strain','MouseID','Type')])
	}

	all.control.dates.df <- cbind(all.control.dates.df, TestDate = format(control.dates.df.tmp, format = '%Y-%m-%d'))
	dates.df <- rbind(all.control.dates.df,mutant.list[,c('Strain','TestDate','Type')])
	dates.df$TestDate <- as.Date(dates.df$TestDate, format = "%Y-%m-%d")
	mouseid.df <- rbind(all.control.mouseid.df, mutant.list[,c('Strain','MouseID','Type')])
	return(mouseid.df)
}
