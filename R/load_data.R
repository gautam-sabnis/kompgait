load_data <- function(){

	Phenos.lin <- c("speed","limb_duty_factor","step_length1","step_width","stride_length",
	"temporal_symmetry", "base_tail_lateral_displacement","tip_tail_lateral_displacement",
	"nose_lateral_displacement")
	Phenos.lin.Nomen <- c("Speed","Limb Duty Factor","Step Length","Step Width","Stride Length",
	"TS","Base Tail LD","Tip Tail LD","Nose LD")
	Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase",
	"nose_lateral_displacement_phase")

	setwd("/Users/sabnig/my-projects/paper__komp")
	data_per_stride <- read.delim("Data/kompdf-corr", stringsAsFactors = TRUE)
	names(data_per_stride)[names(data_per_stride) == 'Mouse.ID'] <- 'MouseID'
	names(data_per_stride)[names(data_per_stride) == 'Date.of.Birth'] <- 'DOB'
	names(data_per_stride)[names(data_per_stride) == 'OFA_Date.of.test.New'] <- 'TestDate'
	names(data_per_stride)[names(data_per_stride) == 'OFA_Genotype'] <- 'Strain'
	names(data_per_stride)[names(data_per_stride) == 'speed_cm_per_sec'] <- 'speed'
	data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')] <- lapply(data_per_stride[,names(data_per_stride) %in% c('MouseID','Strain','Sex','TestAge')], function(x) as.factor(x))
	levels(data_per_stride$Strain)[1] <- "C57BL/6NJ"
	levels(data_per_stride$Strain)[119] <- "Mrps22<tm1.1(KOMP)Vlcg> -/+"


	#Remove Strains
	toMatch <- c("Esrrb", "<em2J>/J COIN", "IMPC")
	matches <- unique(grep(paste(toMatch, collapse = "|"), data_per_stride$Strain, value = TRUE))
	Strains <- setdiff(unique(data_per_stride$Strain), matches)
	data_per_stride <- data_per_stride[data_per_stride$Strain %in% Strains, ]

	#Focus on certain speed bins 
	data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% 
	vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]
	data_per_stride <- data_per_stride[data_per_stride$MouseID %in% 
	names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) >= 
		summary(as.numeric(table(data_per_stride$MouseID)))[2]]),]
	data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
	data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
	data_per_stride$Strain <- as.factor(data_per_stride$Strain)
	data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
	data_per_animal <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
	Strain <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Strain'][1])
	TestDate <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestDate'][1])
	TestAge <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'TestAge'][1])
	BodyLength <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'BodyLength'][1])
	Sex <- sapply(seq(dim(data_per_animal)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal$MouseID[x], 'Sex'][1])
	
	data_per_animal <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, data_per_animal)
	data_per_animal <- data_per_animal[,c("Strain","MouseID","Sex","TestAge","TestDate","BodyLength",Phenos.lin)]
}