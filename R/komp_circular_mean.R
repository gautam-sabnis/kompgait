#' Fit a circular-linear regression model
#'
#' Fit a circular-linear regression model to circular posture 
#' phenotypes 
#' @param data data frame containing per animal data
#' @param CtrlStrain specify the control strain
#' @param model specify the model  
#' @return a list containing the results of the model 
#' @examples 
#' komp_circular_mean(data, CtrlStrain = "C57BL/6NJ", model = "M3")
#' @export 

komp_circular_mean <- function(data, CtrlStrain="C57BL/6NJ", model){

	Phenos.circ.Nomen <- c('Base Tail Phase', 'Tip Tail Phase', 'Nose Phase')
	data_per_animal <- data
	Mutants <- setdiff(unique(data_per_animal$Strain), "C57BL/6NJ")
	pvalGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.circ)))
	esizeGen <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.circ)))
	pvalSex <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.circ)))
	esizeSex <- data.frame(matrix(0,nrow=length(Mutants),ncol=length(Phenos.circ)))
	names(pvalGen) <- Phenos.circ
	names(esizeGen) <- Phenos.circ
	names(pvalSex) <- Phenos.circ
	names(esizeSex) <- Phenos.circ
	for (s in 1:length(Mutants)){
		setTimeLimit(60)
		cat("Analyzing Strain", paste0(Mutants[s]), "\n")
		CtrlIDs <- unique(subset(controlids.df,Strain == Mutants[s])$MouseID)
		df <- data_per_animal[data_per_animal$MouseID %in% CtrlIDs,c('MouseID','Strain','Sex','TestAge','TestDate','BodyLength','speed',Phenos.circ)]	
		df[names(df) %in% Phenos.circ] <- lapply(df[names(df) %in% Phenos.circ], function(x) 
			circular(x, type = 'directions', units = 'radians'))
		df['Genotype'] <- ifelse(df$Strain == CtrlStrain, 'Control','Mutant')
    	df$Genotype <- relevel(factor(df$Genotype), ref = "Control")
		#df <- df[df$TestDate %in% names(which(table(df$TestDate, df$Genotype)[,2] >= 1)), ]
		df$Strain <- droplevels(df$Strain)
		df$TestDate <- droplevels(df$TestDate)
		df$Sex <- relevel(factor(df$Sex), ref = "Male")
    	df$TestAge <- as.factor(df$TestAge)
    	df$BodyLength <- (df$BodyLength - mean(df$BodyLength))/sd(df$BodyLength)
    	df[,sapply(df,is.numeric)] <- apply(df[,sapply(df,is.numeric)],2, function(x) (x-mean(x))/sd(x))
    	if (model=='M3')
    	{X <- cbind(df$BodyLength, df$Genotype, df$speed, df$Sex)} else if (model == 'M2')
    	{X <- cbind(df$speed, df$Genotype, df$Sex)} else 
    	{X <- cbind(df$BodyLength, df$Genotype, df$Sex)}

    	fits <- lapply(Phenos.circ, function(p) {
		tryCatch({
			tmp <- suppressWarnings(lm.circular(y=df[[paste0(p)]], x = X, init=rep(0,ncol(X)), type='c-l'));
			pvalGen[s,p] <<- tmp$p.values[2];
			esizeGen[s,p] <<- tmp$coefficients[2];
			if (model == "M3"){
				pvalSex[s,p] <<- tmp$p.value[4];
				esizeSex[s,p] <<- tmp$coefficients[4];
			} else {
				pvalSex[s,p] <<- tmp$p.value[3];
				esizeSex[s,p] <<- tmp$coefficients[3];
			}
			}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
		})
    }	
	
	pvalGen[pvalGen == 0] <- NA #For lines where the algorithm didn't converge
	
	#Genotype
	esizeGen <- as.matrix(esizeGen)
	rownames(pvalGen) <- Mutants
	colnames(pvalGen) <- Phenos.circ.Nomen
	rownames(esizeGen) <- Mutants
	colnames(esizeGen) <- Phenos.circ.Nomen
	pvalGenadj <- sapply(Phenos.circ.Nomen, function(r) p.adjust(pvalGen[,r],method = 'fdr'))
	pvalGenSignif <- apply(pvalGenadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodGen <- apply(pvalGenadj,2, function(x) -log(x))
    lodGen[is.na(lodGen)] <- 0
    rownames(pvalGenSignif) <- Mutants
    rownames(lodGen) <- Mutants
	colnames(lodGen) <- Phenos.circ.Nomen
	lodGen <- lodGen[order(row.names(lodGen)),]
	esizeGen <- esizeGen[order(row.names(esizeGen)),]
	pvalGenSignif <- pvalGenSignif[order(row.names(pvalGenSignif)),]
	lodGen2 <- lodGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	esizeGen2 <- esizeGen[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]
	pvalGenSignif2 <- pvalGenSignif[sapply(seq(nrow(lodGen)), function(x) any(lodGen[x,] >= 2.98)),]

	#Mutants <- rownames(lodGen2)
	row <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][1]
	col <- which(abs(esizeGen2)==max(abs(esizeGen2)), arr.ind=TRUE)[1:2][2]
	esizeGen2[row,col] <- tail(sort(abs(esizeGen2)),1)[1]
	summ.df.Gen <- data.frame(Phenos.circ.Nomen, 
	pink = sapply(seq(Phenos.circ.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99 & abs(esizeGen2[,x]) >= 0.2))), 
	black = sapply(seq(Phenos.circ.Nomen), function(x) sum(as.numeric(lodGen2[,x] > 2.99))))

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
    	heatmap_legend_param = list(at = c(min(-1,round(min(esizeGen2),2)),0,max(1,round(max(esizeGen2),2))),title = "Standardized Effect Size", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.circle(x = x, y = y, r = (esizeGen2[i, j]) * 1 * (min(unit.c(width, height))),
            gp = gpar(fill = col_fun((esizeGen2[i,j]))))}, rect_gp = gpar(type = "none"))

	ht.Gen <- ht.Gen.p + ht.Gen.e

	#Sex
	pvalSex[pvalSex == 0] <- NA #For lines where the algorithm didn't converge
	esizeSex <- as.matrix(esizeSex)
	rownames(pvalSex) <- Mutants
	colnames(pvalSex) <- Phenos.circ.Nomen
	rownames(esizeSex) <- Mutants
	colnames(esizeSex) <- Phenos.circ.Nomen
	pvalSexadj <- sapply(Phenos.circ.Nomen, function(r) p.adjust(pvalSex[,r],method = 'fdr'))
	pvalSexSignif <- apply(pvalSexadj, 2, function(j) symnum(j, corr = FALSE, na = FALSE, cutpoints = c(0, 
    	0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "+", " ")))
    lodSex <- apply(pvalSexadj,2, function(x) -log(x))
    lodSex[is.na(lodSex)] <- 0
    rownames(pvalSexSignif) <- Mutants
    
    rownames(lodSex) <- Mutants
	colnames(lodSex) <- Phenos.circ.Nomen
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

    summ.df.Sex <- data.frame(Phenos.circ.Nomen, 
	pink = sapply(seq(Phenos.circ.Nomen), function(x) sum(as.numeric(lodSex2[,x] > 2.99 & abs(esizeSex2[,x]) >= 0.1))), 
	black = sapply(seq(Phenos.circ.Nomen), function(x) sum(as.numeric(lodSex2[,x] > 2.99))))

	pinkStrainsSex <- names(which(apply((lodSex2 > 2.99 | lodSex2 > 2.99) & (abs(esizeSex2) >= 0.1 | abs(esizeSex2) >= 0.1),1,any)))
	blackStrainsSex <- names(apply(lodGen2 > 2.99,1,any))
	col_fun <- colorRamp2(c(2.98,2.99,4,5),c("#ffffbf","#abdda4","#66c2a5","#3288bd"))
	ht.Sex.p <- Heatmap((lodSex2), row_names_gp = gpar(fontsize = 16, fontface='italic',
		col=ifelse(apply(((lodSex2 > 2.99 | lodSex2 > 2.99) & (abs(esizeSex2) >= 0.2 | abs(esizeSex2 >= 0.2))),1, function(x) any(x == TRUE)), '#dd3497', 'black')),row_names_side = "left", column_names_gp = gpar(fontsize = 16),
    	heatmap_legend_param = list(at = c(2.98,2.99,4,5),title = "-log10(q-value)", title_position = "leftcenter-rot",  
        border = "black",legend_height = unit(4, "cm"), just = c("right", "top")), col = col_fun, 
    	cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,layer_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))
        grid.text((pvalSexSignif2[i,j]),x, y, gp = gpar(fontsize = 8))}, column_title=paste0(model,': Mean Phenotypes'))

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