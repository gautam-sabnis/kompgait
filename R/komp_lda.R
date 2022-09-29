#' Linear discriminant analysis
#' 
#' Project the data from the 9D gait space to a 2D space
#' 
#' @param data_animal a csv file containing data per animal
#' @param Mutants a vector of mutants  
#' @param controlids.df data frame containing appropriate controls 
#' @return a plot containing LDA space and gait feature 
#' contributions to the LD axes and vignette plots   
#' @examples 
#' komp_lda(data_animal = data_animal, Mutants = 
#' c('Arpc5l-/+','Fam120b-/+'))
#' @export 

komp_lda <- function(data_animal, Mutants, controlids.df){
	data_per_animal <- data_animal
	df <- data.frame()
	CtrlStrain <- "C57BL/6NJ"
	for (m in seq(Mutants)){
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
		df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
		df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
    	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
		df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 1)), ]
		df1$Strain <- droplevels(df1$Strain)
		df <- rbind(df,df1)
	}
	df <- unique(df)
	df <- df[complete.cases(df),]
	df$Outlier <- as.factor(sapply(seq(nrow(df)), function(x) df.out[df.out$MouseID == df$MouseID[x], 'Outlier']))
	df[df$Strain == CtrlStrain, 'Outlier'] <- 0 
	MutStrain <- setdiff(levels(df$Strain),CtrlStrain)
	df$Strain <- factor(df$Strain, levels = c(CtrlStrain,MutStrain[1],MutStrain[2]), ordered=TRUE)
	FReffects <- 'BodyLength'
	formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
   	fits <- lapply(formulas, lm, data = df)
   	df_resid <- data.frame(sapply(seq(Phenos.lin), function(x) resid(fits[[x]])))
   	colnames(df_resid) <- Phenos.lin
   	df_resid <- cbind(Strain = df$Strain, df_resid)
   	df_lda <- df_resid
   	df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])
	fit_lda <- lda(Strain ~ ., data = df_lda)
	lda_values <- predict(fit_lda, df_lda[,-1])
	C <- ggplot(data = data.frame(Strain = df$Strain, lda_values$x), aes(x=LD1,y=LD2,shape=Strain,color=Strain,fill=Strain)) + geom_point(size = 5,aes(color=Strain)) + stat_ellipse(geom = "polygon", alpha = 1/3, aes(fill=Strain)) + theme_bw(base_size = 25) + theme(legend.position = "top") +scale_color_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) + scale_fill_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) 
	
	PCLD_df <- as.data.frame(fit_lda$scaling)
	rownames(PCLD_df) <- Phenos.lin
	PCLD1 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,1]))
	PCLD1$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 
		'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
	CX <- ggplot(PCLD1, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = 'Loadings') + ggtitle('LD1')
	PCLD2 <- data.frame(Phenos = Phenos.lin.Nomen, value = abs(PCLD_df[,2]))
	PCLD2$Phenos <- factor(PCLD1$Phenos, levels = c('Speed','Limb Duty Factor', 'Step Length', 'Step Width', 'Stride Length', 'TS', 'Base Tail LD', 'Tip Tail LD', 'Nose LD'))
	CY <- ggplot(PCLD2, aes(x = Phenos, y = value)) + geom_bar(stat = 'identity', color = 'black') + theme_bw(base_size = 16) +
	theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + labs(x = NULL, y = NULL) + ggtitle('LD2')

	phenotype <- PCLD1$Pheno[which.max(PCLD1$value)]
	phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
	p1 <- ggplot(df_lda, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a")))


	phenotype <- PCLD2$Pheno[which.max(PCLD2$value)]
	phenoname <- Phenos.lin[which(Phenos.lin.Nomen == phenotype)]
	p2 <- ggplot(df_lda, aes_string(x = 'Strain', y = phenoname, fill = 'Strain')) + geom_boxplot(alpha = 1/3, outlier.shape = NA) + theme_bw(base_size = 16) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(y = paste0(phenotype," ", "(Residuals)")) + scale_fill_manual(values = c("C57BL/6NJ" = "#e41a1c",assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a")))
	

	layout_matrix <- rbind(c(1,1,1,2,3),c(1,1,1,4,5))
	p <- gridExtra::grid.arrange(C,CX,CY,p1,p2,layout_matrix = layout_matrix)

	return(p)
}


komp_lda_supp <- function(Mutants){
	df <- data.frame()
	CtrlStrain <- "C57BL/6NJ"
	for (m in seq(Mutants)){
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[m])$MouseID)
		df1 <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,]
		df1['Genotype'] <- ifelse(df1$Strain == CtrlStrain, 'Control','Mutant')
    	df1$Genotype <- relevel(factor(df1$Genotype), ref = "Control")
		df1 <- df1[df1$TestDate %in% names(which(table(df1$TestDate, df1$Genotype)[,2] >= 1)), ]
		df1$Strain <- droplevels(df1$Strain)
		df <- rbind(df,df1)
	}
	df <- unique(df)
	df <- df[complete.cases(df),]
	df$Outlier <- as.factor(sapply(seq(nrow(df)), function(x) df.out[df.out$MouseID == df$MouseID[x], 'Outlier']))
	df[df$Strain == CtrlStrain, 'Outlier'] <- 0 
	MutStrain <- setdiff(levels(df$Strain),CtrlStrain)
	df$Strain <- factor(df$Strain, levels = c(CtrlStrain,MutStrain[1],MutStrain[2]), ordered=TRUE)
	FReffects <- 'BodyLength'
	formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
   	fits <- lapply(formulas, lm, data = df)
   	df_resid <- data.frame(sapply(seq(Phenos.lin), function(x) resid(fits[[x]])))
   	colnames(df_resid) <- Phenos.lin
   	df_resid <- cbind(Strain = df$Strain, df_resid)
   	df_lda <- df_resid
   	df_lda <- data.frame(Strain = df$Strain, df[,names(df) %in% Phenos.lin])
	fit_lda <- lda(Strain ~ ., data = df_lda)
	lda_values <- predict(fit_lda, df_lda[,-1])
	C <- ggplot(data = data.frame(Strain = df$Strain, lda_values$x), aes(x=LD1,y=LD2,shape=Strain,color=Strain,fill=Strain)) + geom_point(size = 5,aes(color=Strain)) + stat_ellipse(geom = "polygon", alpha = 1/3, aes(fill=Strain)) + theme_bw(base_size = 25) + theme(legend.position = 'top') + scale_color_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) + scale_fill_manual(values = c(assign(CtrlStrain,"#e41a1c"),assign(MutStrain[1],"#377eb8"),assign(MutStrain[2],"#4daf4a"))) 
	return(C)
}