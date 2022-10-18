#' Multivariate outlier detection  
#' 
#' Use a PCA-based projection outlier algorithm to identify outlier #' mutant lines
#' @param data data.frame containing data at the strain level 
#' @return a plot showing outlier mutant lines in red 
#' @examples 
#' komp_outlier_strain(data = data_per_strain)
#' @export 

komp_outlier_strain <- function(data){

	df <- data
	tmp <- mvoutlier::pcout(df[,names(df) %in% Phenos.lin])
	df.out <- data.frame(Distance1 = tmp$x.dist1, Distance2 = tmp$x.dist2, Strain = df[,'Strain'],
		Outlier = as.numeric(!tmp$wfinal01),WeightL = tmp$wloc, WeightS = tmp$wscat)
	df.out[df.out$Strain == 'C57BL/6NJ', 'Outlier'] <- -1
	df.out$Outlier <- as.factor(df.out$Outlier)

	p <- ggplot(df.out, aes(x = Distance1, y=Distance2)) + geom_point(alpha=0.8, aes(color=Outlier), size=4) + ggrepel::geom_text_repel(aes(label=ifelse(Outlier %in% c(-1,1),as.character(Strain),'')),size=5,box.padding=1) + labs(x = TeX("Distance$_\\1$"), y = TeX("Distance$_\\2")) + ggtitle('KOMP mutant outliers') + scale_color_manual(values=c("black","grey50","red")) + theme_bw(base_size = 15) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title = element_text(size = 22), legend.position='none')

	return(list(p, df.out))

}

#' Multivariate outlier detection  
#' 
#' Use a PCA-based projection outlier algorithm to identify outlier #' animal
#' @param data_per_strain a data.frame containing data at the strain #' level 
#' @param data_per_animal a data.frame containing data at the animal #' level 
#' @return a plot showing outlier animals in red 
#' @examples 
#' komp_outlier_strain(data_per_strain = data_per_strain, 
#' data_per_animal = data_per_animal)
#' @export 

komp_outlier_animal <- function(data_per_strain, data_per_animal){
	df <- data.frame(Strain = data_per_animal$Strain, MouseID = data_per_animal$MouseID, Sex = data_per_animal$Sex)
	df <- cbind(df, apply(data_per_animal[,names(data_per_animal) %in% Phenos.lin],2,function(x) (x - mean(x))/sd(x)))
	tmp <- mvoutlier::pcout(df[,names(df) %in% Phenos.lin],makeplot=FALSE)
	df.out <- data.frame(Distance1 = tmp$x.dist1, Distance2 = tmp$x.dist2, 
	Label = paste0(df[,'MouseID']," (", df[,'Strain'], ")"), MouseID = df[,'MouseID'], Strain = df[,'Strain'], Outlier = as.numeric(!tmp$wfinal01))
	df.out[(df.out$Distance2 >= 5.9 | df.out$Distance1 > 10), 'Outlier'] <- -1
	df.out[df.out$MouseID == 'J80962', 'Outlier'] <- -1
	df.out[df.out$MouseID == 'J76119', 'Outlier'] <- -1
	df.out$Outlier <- as.factor(df.out$Outlier)
	options(ggrepel.max.overlaps = Inf)
	p <- ggplot(df.out, aes(x = Distance1, y=Distance2)) + geom_point(alpha=0.8, aes(color=Outlier), size=4) + ggrepel::geom_text_repel(aes(label=ifelse(Outlier==-1,as.character(Label),'')),size=6,box.padding=1) + labs(x = TeX("Distance$_\\1$"), y = TeX("Distance$_\\2")) + ggtitle('KOMP Outliers') + scale_color_manual(values=c("red","black","grey50")) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title = element_text(size = 24), plot.title = element_text(size = 26), legend.position='none') + ggtitle('KOMP animal outliers')
	return(list(p, df.out))

}
