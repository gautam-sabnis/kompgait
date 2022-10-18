#' Post hoc comparisons  
#' 
#' Comparing univariate and multivariate analyses to each other and 
#' with existing findings   
#' @param M1 a list obtained using komp_lmer_mean(data, CtrlStrain,
#' model = "M1")
#' @param M2 a list obtained using komp_lmer_mean(data, CtrlStrain,
#' model = "M2")
#' @param M3 a list obtained using komp_lmer_mean(data, CtrlStrain,
#' model = "M3")
#' @param MvOut a list obtained using komp_outlier_strain(data)
#' @return  Plot comparing the number of pink mutants 
#' across different models and multivariate outliers
#' @examples 
#' post_hoc_comparison(M1, M2, M3, MvOut, type = "bar_plot")
#' @export 

post_hoc_comparison <- function(M1, M2, M3, MvOut, type){

	if (type == "bar_plot"){
		Strains1 <- M1[[3]]
		Strains2 <- M2[[3]]
		Strains3 <- M3[[3]]
		StrainsMv <- MvOut[[2]][MvOut[[2]][,"Outlier"] == 1, "Strain"]

		df <- data.frame(Analysis = c('M1','M2','M3','Mvoutlier'), Proportion = c(length(intersect(Strains1,StrainsMv))/length(StrainsMv),length(intersect(Strains2,StrainsMv))/length(StrainsMv),length(intersect(Strains3,StrainsMv))/length(StrainsMv),length(intersect(StrainsMv,StrainsMv))/length(StrainsMv)), Total = c(length(Strains1), length(Strains2), length(Strains3), length(StrainsMv))) 
		p <- ggplot(df, aes(x=Analysis,y=Proportion)) + geom_bar(stat='identity') + theme_bw(base_size=22) +  theme( axis.text.x.top = element_text(vjust=0.5)) + labs(y = 'Proportion') + geom_text(label = paste0(df$Total),vjust = -.01, hjust = 0.5,size=7)
		return(p)
	} else if (type == "upset"){
		Strains1 <- M1[[3]]
		Strains2 <- M2[[3]]
		Strains3 <- M3[[3]]
		StrainsMv <- MvOut[[2]][MvOut[[2]][,"Outlier"] == 1, "Strain"]

		AllStrains <- union(union(union(Strains1,Strains2),Strains3),StrainsMv)
		df <- data.frame(Strains = AllStrains, M1 = rep(0,length(AllStrains)), M2 = rep(0,length(AllStrains)), M3 = rep(0,length(AllStrains)), Mvout = rep(0,length(AllStrains)))
		df$M1 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains1, 1, 0))
		df$M2 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains2, 1, 0))
		df$M3 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains3, 1, 0))
		df$Mvout <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% StrainsMv, 1, 0))
		p <- UpSetR::upset(df, text.scale = 2, point.size = 2, line.size = 1)
		return(p)
	} else if (type == "table"){
		Strains1 <- M1[[3]]
		Strains2 <- M2[[3]]
		Strains3 <- M3[[3]]
		StrainsMv <- MvOut[[2]][MvOut[[2]][,"Outlier"] == 1, "Strain"]

		AllStrains <- union(union(union(Strains1,Strains2),Strains3),StrainsMv)
		df <- data.frame(Strains = AllStrains, M1 = rep(0,length(AllStrains)), M2 = rep(0,length(AllStrains)), M3 = rep(0,length(AllStrains)), Mvout = rep(0,length(AllStrains)))
		df$M1 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains1, 1, 0))
		df$M2 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains2, 1, 0))
		df$M3 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains3, 1, 0))
		df$Mvout <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% StrainsMv, 1, 0))
		return(df)
	}
}