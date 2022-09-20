#' Fit a GLMM to linear gait traits 
#' 
#' Fit a generalized (gamma) linear mixed model to identify significant differences in trait's dispersion value (sd) between
#' mutant and control strain
#' @param CtrlStrain specify the control strain
#' @param model specify the model  
#' @return a list containing the results of the model 
#' @examples 
#' komp_lmer_mean(CtrlStrain = "C57BL/6NJ", model = "M3")
#' @export 

komp_glmer_dispersion <- function(CtrlStrain="C57BL/6NJ", model){	
	if (model == 'M2' || model == 'M3'){
		Phenos.lin <- setdiff(Phenos.lin,"speed");
		Phenos.lin.Nomen <- setdiff(Phenos.lin.Nomen,"Speed");
	} else {
		Phenos.lin <- Phenos.lin;
		Phenos.lin.Nomen <- Phenos.lin.Nomen;
	}

	Mutants <- setdiff(unique(data_per_animal$Strain),"C57BL/6NJ")
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	pvalSex <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	esizeSex <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.lin)))
	
	for (s in 1:(length(Mutants))){
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		if (model == 'M2' || model == 'M3'){
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength','speed',Phenos.lin)]	
		} else {
			df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength',Phenos.lin)]
		
		}
		df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    	df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
		if (sum(table(df$TestDate,df$Genotype)[,1] & table(df$TestDate,df$Genotype)[,2] >= 1) >= 1){
    		df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 1)), ]
    	} else {
    		df <- df
    	}
		df$Strain <- droplevels(df$Strain)
		df$TestDate <- droplevels(df$TestDate)
		df$Sex <- relevel(factor(df$Sex), ref = "Male")
    	df$TestAge <- as.factor(df$TestAge)
    	df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)
    	if (model=='M3')
    		{FReffects.lmer <- 'BodyLength + speed + Sex + Genotype + (1|TestDate)';
    		 FReffects <- 'BodyLength + speed + Sex + Genotype';
    		 df$speed <- (df$speed - mean(df$speed,na.rm=TRUE))/sd(df$speed,na.rm=TRUE)}
    	else if (model == 'M2')
    		{FReffects.lmer <- 'speed + Sex + Genotype + (1|TestDate)';
    		 FReffects <- 'speed + Sex + Genotype'}
    	else {FReffects.lmer <- 'BodyLength + Sex + Genotype + (1|TestDate)';
    		 FReffects <- 'BodyLength + Sex + Genotype'}
    	if (length(levels(df$TestDate)) >= 2){
    		formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects.lmer), simplify=TRUE)) 
    		fits <- sapply(formulas, function(x) glmer(formula=x, data = df, family = Gamma(link = "log"), control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))),simplify=FALSE)
    		fits_removeGenotype <- sapply(seq(length(fits)), function(x) update(fits[[x]], .~ .- Genotype))
    		fits_removeSex <- sapply(seq(length(fits)), function(x) update(fits[[x]], .~ .- Sex))
  			
    		pvalGen[s,] <- sapply(seq(length(fits)), function(x) anova(fits[[x]],fits_removeGenotype[[x]])[["Pr(>Chisq)"]][2]) 
			esizeGen[s,] <- sapply(seq(length(fits)), function(x) anova(fits[[x]],fits_removeSex[[x]])[["Pr(>Chisq)"]][2]) 
			if (min(table(df$Genotype,df$Sex)['Mutant',]) > 1){
				pvalGen[s,] <- sapply(seq(length(fits)), function(x) anova(fits[[x]],fits_removeGenotype[[x]])[["Pr(>Chisq)"]][2])  
				esizeSex[s,] <- sapply(seq_along(fits), function(x) S(fits[[x]])$fixed.effects['SexFemale','Estimate']) 
			} else {
				pvalSex[s,] <- rep(1,length(Phenos.lin)) 
				esizeSex[s,] <- rep(0,length(Phenos.lin)) 
			}
			
    	} else {
    		formulas <- unname(sapply(Phenos.lin ,function(x) paste(x, "~", FReffects), simplify=TRUE)) 
    		fits <- lapply(formulas, lm, data = df)
    		pvalGen[s,] <- sapply(seq_along(fits), function(x) Anova(fits[[x]])['Genotype','Pr(>F)']) 
			esizeGen[s,] <- sapply(seq_along(fits), function(x) unname(fits[[x]]$coefficients['GenotypeMutant'])) 
			pvalSex[s,] <- sapply(seq_along(fits), function(x) Anova(fits[[x]])['Sex','Pr(>F)']) 
			esizeSex[s,] <- sapply(seq_along(fits), function(x) unname(fits[[x]]$coefficients['SexFemale'])) 

    	}
	}

	pvalGen <- pvalGen[complete.cases(pvalGen),]
	esizeGen <- esizeGen[complete.cases(esizeGen),]
	pvalSex <- pvalSex[complete.cases(pvalSex),]
	esizeSex <- esizeSex[complete.cases(esizeSex),]
	
	#Genotype
	Mutants <- setdiff(unique(data_per_animal$Strain),"C57BL/6NJ")
	esizeGen <- as.matrix(esizeGen)
	rownames(pvalGen) <- Mutants
	colnames(pvalGen) <- Phenos.lin.Nomen
	rownames(esizeGen) <- Mutants
	colnames(esizeGen) <- Phenos.lin.Nomen
	pvalGenadj <- sapply(Phenos.lin.Nomen, function(r) p.adjust(pvalGen[,r],method = 'fdr'))
	pvalGenSignif <- apply(pvalGenadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodGen <- apply(pvalGenadj,2, function(x)-log(x))
    rownames(pvalGenSignif) <- Mutants
    
    rownames(lodGen) <- Mutants
	colnames(lodGen) <- Phenos.lin.Nomen
	lodGen <- lodGen[order(row.names(lodGen)),]
	esizeGen <- esizeGen[order(row.names(esizeGen)),]
	pvalGenSignif <- pvalGenSignif[order(row.names(pvalGenSignif)),]
	lodGen2 <- lodGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	esizeGen2 <- esizeGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	pvalGenSignif2 <- pvalGenSignif[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	Mutants <- rownames(lodGen2)
	row <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][1]
	col <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][2]
	esizeGen2[row,col] <- tail(sort(abs(esizeGen2)),1)[1]

	summ.df.Gen <- data.frame(Phenos.lin.Nomen, 
	pink = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99 & abs(esizeGen2[,x]) >= 0.2))), 
	black = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99))))

	pinkStrainsGen <- names(which(apply((lodGen2 > 2.99 | lodGen2 > 2.99) & (abs(esizeGen2) >= 0.2 | abs(esizeGen2) >= 0.2),1,any)))
	blackStrainsGen <- names(apply(lodGen2 > 2.99,1,any))

	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Gen.p <- Heatmap((lodGen2), row_names_gp = gpar(fontsize = 16, fontface='italic',
		col=ifelse(apply(((lodGen2 > 2.99 | lodGen2 > 2.99) & (abs(esizeGen2) >= 0.2 | abs(esizeGen2 >= 0.2))),1, function(x) any(x == TRUE)), '#dd3497', 'black')),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalGenSignif2[i,j]),x, y, gp = gpar(fontsize = 8))}, column_title=paste0(model,': Mean Phenotypes'))

	col_fun <- circlize::colorRamp2(c(min(esizeGen2),0,max(esizeGen2)),c("#542788","#f7f7f7","#d73027"))
	ht.Gen.e <- Heatmap((esizeGen2), row_names_gp = gpar(fontsize = 16, fontface="italic"),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(round(min(esizeGen2),2),0,round(max(esizeGen2),2)),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeGen2[i, j]) * 0.8 * (max(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen2[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Gen <- ht.Gen.p + ht.Gen.e

	#Sex
	Mutants <- setdiff(unique(data_per_animal$Strain),c("C57BL/6NJ"))
	esizeSex <- as.matrix(esizeSex)
	rownames(pvalSex) <- Mutants
	colnames(pvalSex) <- Phenos.lin.Nomen
	rownames(esizeSex) <- Mutants
	colnames(esizeSex) <- Phenos.lin.Nomen
	pvalSexadj <- sapply(Phenos.lin.Nomen, function(r) p.adjust(pvalSex[,r],method = 'fdr'))
	pvalSexSignif <- apply(pvalSexadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodSex <- apply(pvalSexadj,2, function(x)-log(x))
    rownames(pvalSexSignif) <- Mutants
    
    rownames(lodSex) <- Mutants
	colnames(lodSex) <- Phenos.lin.Nomen
	lodSex <- lodSex[order(row.names(lodSex)),]
	esizeSex <- esizeSex[order(row.names(esizeSex)),]
	pvalSexSignif <- pvalSexSignif[order(row.names(pvalSexSignif)),]
	lodSex2 <- lodSex[sapply(seq(nrow(lodSex)), function(x) any(lodSex[x,] >= 2.98)),]
	esizeSex2 <- esizeSex[sapply(seq(nrow(lodSex)), function(x) any(lodSex[x,] >= 2.98)),]
	pvalSexSignif2 <- pvalSexSignif[sapply(seq(nrow(lodSex)), function(x) any(lodSex[x,] >= 2.98)),]
	Mutants <- rownames(lodSex2)
	row <- which(abs(esizeSex2)==max(abs(esizeSex2)), arr.ind=TRUE)[1:2][1]
	col <- which(abs(esizeSex2)==max(abs(esizeSex2)), arr.ind=TRUE)[1:2][2]
	esizeSex2[row,col] <- tail(sort(abs(esizeSex2)),1)[1]

    summ.df.Sex <- data.frame(Phenos.lin.Nomen, 
	pink = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodSex2[,x] > 2.99 & abs(esizeSex2[,x]) >= 0.1))), 
	black = sapply(seq(Phenos.lin.Nomen), function(x) sum(as.numeric(lodSex2[,x] > 2.99))))

	pinkStrainsSex <- names(which(apply((lodSex2 > 2.99 | lodSex2 > 2.99) & (abs(esizeSex2) >= 0.1 | abs(esizeSex2) >= 0.1),1,any)))
	blackStrainsSex <- names(apply(lodGen2 > 2.99,1,any))
	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Sex.p <- Heatmap((lodSex2), row_names_gp = gpar(fontsize = 16, fontface='italic',
		col=ifelse(apply(((lodSex2 > 2.99 | lodSex2 > 2.99) & (abs(esizeSex2) >= 0.2 | abs(esizeSex2 >= 0.2))),1, function(x) any(x == TRUE)), '#dd3497', 'black')),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,layer_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalSexSignif2[i,j]),x, y, gp = gpar(fontsize = 8))}, column_title=paste0(model,': Variance Phenotypes'))

	col_fun <- circlize::colorRamp2(c(min(esizeSex2),0,max(esizeSex2)),c("#542788","#f7f7f7","#d73027"))
	ht.Sex.e <- Heatmap((esizeSex2), row_names_gp = gpar(fontsize = 16, fontface="italic"),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(round(min(esizeSex2),2),0,round(max(esizeSex2),2)),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, layer_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeSex2[i, j]) * (max(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeSex2[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Sex <- ht.Sex.p + ht.Sex.e


	return(list(ht.Gen,summ.df.Gen,pinkStrainsGen,blackStrainsGen,ht.Sex,summ.df.Sex,pinkStrainsSex,
		blackStrainsSex,esizeGen2,lodGen2))

}