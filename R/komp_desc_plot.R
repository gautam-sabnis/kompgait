#' Linear gait phenotype distributions across strains 
#' 
#' Creates a plot showing the phenotypic mean ± sd for all strains \code{\link{kompgait}}
#' @param df_animal data frame containing data at per animal level
#' @param df_strain data frame containing data at per strain level
#' @param phenotype choose phenotype for plotting
#' @param phenotype phenotype label
#' @param align horizontal or vertical plot
#' @return plot similar to Figure 2 in the manuscript 
#' @examples 
#' komp_desc_plot()
#' @export 

komp_desc_plot <- function(df_animal, df_strain, phenotype, phenoname, align){
	df_animal$sorted_strain = make_sort_column(df_animal, 'Strain', phenotype ,FUN = mean, decreasing=FALSE)
	me <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenotype)][1]
	std <- df_strain[df_strain$Strain == 'C57BL/6NJ', names(df_strain) %in% c(phenotype)][3]
	df <- as.data.frame(df_strain[,phenotype]) 
	df <- cbind(Strain = df_strain$Strain, df)

	if(align == 'h'){
		ggplot(df_animal, aes_string(x = "sorted_strain", y = phenotype)) + geom_point(alpha = 0.9, size = 3, aes(color=Sex)) + 
		geom_point(data = df, aes(x = Strain, y = df[,"mn"]), color = 'black', size = 4, alpha = 0.5) + 
		scale_color_manual(values=c("#E41A1C", "#377EB8")) + 
		geom_errorbar(data = df, aes(x = Strain, df[,"mn"], ymin = df[,"mn"]-df[,"std"],
		ymax = df[,"mn"]+df[,"std"]), color = 'grey30') + geom_hline(yintercept = me - std, 
		color = "black", size = 1) +  geom_hline(yintercept = me, color = "black", size = 1) +  
		geom_hline(yintercept = me + std, color = "black", size = 1) + labs(x = 'Strain', y = paste0(phenoname)) + 
		theme_bw(base_size = 16) + ggtitle('') + theme_bw(base_size = 18) + theme(legend.position = c(0.05,0.85), 
			legend.text = element_text(size=13),legend.background = element_blank(),
			axis.text.x = element_text(angle = 85, vjust = 0.55, size = 12), legend.key.size = unit(0.1, "cm"))        
    } else {
		ggplot(df_animal, aes_string(y = "sorted_strain", x = phenotype)) + geom_point(alpha = 0.9, size = 3, aes(color=Sex)) + 
		geom_point(data = df, aes(y = Strain, x = df[,"mn"]), color = 'black', size = 4, alpha = 0.5) + 
		scale_color_manual(values=c("#E41A1C", "#377EB8")) + 
		geom_errorbarh(data = df, aes(y = Strain, df[,"mn"], xmin = df[,"mn"]-df[,"std"],
		xmax = df[,"mn"]+df[,"std"]), color = 'grey30') + geom_vline(xintercept = me - std, 
		color = "black", size = 1,alpha = 0.6) +  geom_vline(xintercept = me, color = "black", size = 1,alpha = 0.5) +  
		geom_vline(xintercept = me + std, color = "black", size = 1,alpha = 0.6) + labs(y = 'Strain', x = paste0(phenoname)) +  
		theme_bw(base_size = 26) + ggtitle('') + theme(legend.position = c(0.60,0.05), 
			legend.text = element_text(size=14),legend.background = element_blank(),
			axis.text.y = element_text(size = 15,face="italic",hjust = 0.4), legend.key.size = unit(0.1, "cm"))
	}
}
