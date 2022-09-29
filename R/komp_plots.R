#' Make a vignette plot of the results from per line analyses 
#' 
#' Compare mutant versus control for phenotypes identified as 
#' significantly different by models M1-3
#' @param data_animal a csv file containing data per animal
#' @param phenotype full trait name
#' @param phenoname trait name in the dataframe
#' @param esizeGen2 effect size matrix - output of modeling function
#' @param lodGen2 lod scores matrix - output of modeling function
#' @param model character "M1-3" 
#' @return a plot 
#' @examples 
#' vignette_plot(phenotype = "Stride Length", phenoname = "stride 
#' _length", esizeGen2 = dispM3[[9]], lodGen2 = dispM3[[10]], model = 
#' "M3")
#' @export 

vignette_plot <- function(data_animal, phenotype, phenoname, esizeGen2, lodGen2, model){

	data_per_animal <- data_animal
	if (model == 'M2' || model == 'M3'){
		Phenos.lin <- setdiff(Phenos.lin,"speed");
		Phenos.lin.Nomen <- setdiff(Phenos.lin.Nomen,"Speed");
	} else {
		Phenos.lin <- Phenos.lin;
		Phenos.lin.Nomen <- Phenos.lin.Nomen;
	}
	esize <- names(sort(abs(esizeGen2[lodGen2[,paste0(phenotype)] > 2.98,paste0(phenotype)]), decreasing=TRUE)[1:5])
	esize <- esize[!is.na(esize)]
	if (length(esize)!=0){
			dflist <- list();
		invisible(lapply(seq(length(esize)), function(x) {CtrlIDs <- unique(subset(controlids.df,Strain == esize[x])$MouseID);
		dftmp <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('Strain','Sex','TestDate',paste0(phenoname))];
		dftmp['Genotype'] <- ifelse(dftmp$Strain == "C57BL/6NJ", 'Control','Mutant');
		dftmp$Genotype <- relevel(factor(dftmp$Genotype), ref = "Control");
		dftmp$Strain <- droplevels(dftmp$Strain);
		dflist[[x]] <<- dftmp}))
		dftmp <- do.call(rbind, dflist)
		dftmp$Strain <- rep(esize, sapply(seq(length(esize)), function(x) dim(dflist[[x]])[1]))
		p <- ggplot(dftmp, aes_string(y=phenoname, x='Genotype')) + geom_boxplot(outlier.shape=NA,width=0.7) + geom_jitter(aes(color=Genotype), width=0.1, size = 5, shape = 1, stroke=1.5) + scale_color_manual(values=c("#6a51a3", "#d94801")) + theme_bw(base_size = 40) + theme(legend.position='none') + labs(y = paste0(phenotype)) + facet_grid(~Strain)
		return(p)
	}
}

#' Make a bar plot showing the number of pink strains across different
#' analyses
#' A plot showing the number of pink strains found across different  
#' analyses (M1-3, MvOut)
#' @param M1 full trait name
#' @param M2 trait name in the dataframe
#' @param M3 effect size matrix - output of modeling function
#' @param MvOut lod scores matrix - output of modeling function
#' @return a plot 
#' @examples 
#' plot_model_comparisons(M1 = aveM1[[3]], M2 = aveM2[[3]], M3 = 
#' aveM3[[3]], MvOut = mv_outlier_strains[[2]])
#' @export 

plot_model_comparisons <- function(M1, M2, M3, MvOut){

	df.out <- MvOut
	Strains1 <- M1[[3]]
	Strains2 <- M2[[3]]
	Strains3 <- M3[[3]]
	StrainsMv <- unique(df.out[df.out$Outlier == 1, 'Strain'])

	df <- data.frame(Analysis = c('M1','M2','M3','Mvoutlier'), Proportion = c(length(intersect(Strains1,StrainsMv))/length(StrainsMv),length(intersect(Strains2,StrainsMv))/length(StrainsMv),length(intersect(Strains3,StrainsMv))/length(StrainsMv),length(intersect(StrainsMv,StrainsMv))/length(StrainsMv)), Total = c(length(Strains1), length(Strains2), length(Strains3),length(StrainsMv))) 
    p <- ggplot(df, aes(x=Analysis,y=Proportion)) + geom_bar(stat='identity') + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title = element_text(size = 24), plot.title = element_text(size = 26), legend.position='none') + labs(y = 'Proportion') + geom_text(label = paste0(df$Total),vjust = -.01, hjust = 0.5,size=7)
    return(p)
}

#' Make an UpSet plot
#' 
#' A plot showing the common pink strains found across different  
#' analyses (M1-3, MvOut)
#' @param M1 full trait name
#' @param M2 trait name in the dataframe
#' @param M3 effect size matrix - output of modeling function
#' @param MvOut lod scores matrix - output of modeling function
#' @return a plot 
#' @examples 
#' UpSet_plot(M1 = aveM1[[3]], M2 = aveM2[[3]], M3 = 
#' aveM3[[3]], MvOut = mv_outlier_strains[[2]])
#' @export 

UpSet_plot <- function(M1, M2, M3, MvOut){

	df.out <- MvOut
	Strains1 <- M1[[3]]
	Strains2 <- M2[[3]]
	Strains3 <- M3[[3]]
	StrainsMv <- unique(df.out[df.out$Outlier == 1, 'Strain'])
	AllStrains <- union(union(union(Strains1,Strains2),Strains3),StrainsMv)
	df <- data.frame(Strains = AllStrains, M1 = rep(0,length(AllStrains)), M2 = rep(0,length(AllStrains)), M3 = rep(0,length(AllStrains)), Mvout = rep(0,length(AllStrains)))
	Strains1 <- c("Arpc5l-/+","Bcl11b-/+","Rab3gap2-/-","Rad51ap1-/-","Ric8a-/+","Zwint-/+")
	df$M1 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains1, 1, 0))	
	df$M2 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains2, 1, 0))
	df$M3 <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% Strains3, 1, 0))
	df$Mvout <- sapply(seq(nrow(df)), function(x) ifelse(df$Strains[x] %in% StrainsMv, 1, 0))

	p <- UpSetR::upset(df, text.scale = 3, point.size = 2, line.size = 1, order.by = "freq")
	return(p)

}

#' Make a vignette plot for the outlier animal
#' 
#' A plot that compares the outlier animal to other animals in 
#' belonging to the same strain tested within 3 weeks 
#' analyses (M1-3, MvOut)
#' @param data_animal data_per_animal
#' @param df_out a dataframe containing outlier animals
#' @param Mutant Strain that the outlier animal belongs to
#' @param phenotypes body length adjusted phenotypes to use for 
#''	plotting
#' @return a plot 
#' @examples 
#' plot_outlier_animal_comparison(data_animal = data_per_animal,
#' df_out = outlier_animals[[2]], Mutant = "Sfxn5-/+", phenotypes = 
#' c("Speed", "Limb Duty Factor", "Step Length")) 
#' @export 

plot_outlier_animal_comparison <- function(data_animal, df_out, Mutant, phenotypes){

	data_per_animal <- data_animal
	CtrlIDs <- unique(subset(controlids.df,Strain == Mutant)$MouseID)
	dfa <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs, ]
	dfa$Genotype <- ifelse(dfa$Strain == 'C57BL/6NJ', 'Control', 'Mutant')
	dfa$MouseID <- droplevels(dfa$MouseID)
	
	if (Mutant == "Alg11"){
		dfa <- dfa[dfa$TestDate %in% '2016-09-06', ]
	} else {
		dfa <- dfa[dfa$TestDate %in% names(which(table(dfa$TestDate, dfa$Genotype)[,2] >= 1)), ]	
	}

	dfa$TestDate <- droplevels(dfa$TestDate)
	CtrlIDs <- setdiff(unique(dfa[dfa$Genotype == 'Control', 'MouseID']), unique(df_out$MouseID))
	df <- data_per_stride[data_per_stride$Strain == Mutant, c('MouseID','BodyLength','Sex',Phenos.lin)]
	df <- rbind(df,data_per_stride[data_per_stride$MouseID %in% CtrlIDs, c('MouseID', 'BodyLength','Sex',Phenos.lin)])
	df$MouseID <- droplevels(df$MouseID)
	df$Outlier <- as.factor(ifelse(df$MouseID %in% df_out$MouseID, 1, ifelse(df$MouseID %in% CtrlIDs, -1, 0))) 
	if (!("0" %in% levels(df$Outlier))){df$Outlier <- factor(df$Outlier, levels = c("0",levels(df$Outlier)))}
	df2 <- df[,names(df) %in% c(Phenos.lin)]
	df2 <- data.frame(do.call(cbind,lapply(seq(length(Phenos.lin)), function(p) as.numeric(resid(lm(df[[Phenos.lin[p]]] ~ BodyLength, data = df))))))
	names(df2) <- Phenos.lin.Nomen
	df2 <- cbind(id = 1:dim(df)[1], df2)
	df.melt <- reshape::melt(df2, id.vars = 'id')
	df.melt <- cbind(df.melt, MouseID = rep(rep(names(table(df$MouseID)), as.numeric(table(df$MouseID))),length(Phenos.lin)))
	df.melt$Outlier <- as.factor(ifelse(df.melt$MouseID %in% df_out$MouseID, 1, ifelse(df.melt$MouseID %in% CtrlIDs, -1, 0))) 
	if (!("0" %in% levels(df.melt$Outlier))){df.melt$Outlier <- factor(df.melt$Outlier, levels = c("0",levels(df.melt$Outlier)))}
	p <- ggplot(df.melt[df.melt$variable %in% phenotypes,], aes(x=MouseID,y=value)) + geom_boxplot(outlier.shape=NA,aes(fill = Outlier),alpha = 0.4) + facet_wrap(~ variable, scales = 'free') + ggtitle(paste0(Mutant)) + scale_fill_manual(name = 'Genotype', values =c("0" = "black","1" = "#d94801", "-1" = "#6a51a3"), labels = c("0"='Mutant',"-1"='Control',"1"="Mutant (Outlier)"), drop = FALSE) + theme_bw(base_size=22) + theme(axis.text.x=element_text(angle=90,size=16), legend.position='none') + labs(y = 'Residuals') 

	return(p)

}


