#' Load and preprocess data 
#' 
#' Loads the data and preprocess the features. 
#' 1) Remove certain strains
#' 2) Focus on certain speed bins
#' @param data a csv file containing data per stride
#' @param type linear or circular 
#' @return a list containing data_per_stride and data_per_animal 
#' @examples 
#' load_data(data = data_per_stride)
#' @export 

load_data <- function(data, type){

	data_per_stride <- data
	Phenos.lin <- c("speed","limb_duty_factor","step_length1","step_width","stride_length","temporal_symmetry", "base_tail_lateral_displacement","tip_tail_lateral_displacement","nose_lateral_displacement")
	Phenos.lin.Nomen <- c("Speed","Limb Duty Factor","Step Length","Step Width","Stride Length","TS","Base Tail LD","Tip Tail LD","Nose LD")
	Phenos.circ <- c("base_tail_lateral_displacement_phase","tip_tail_lateral_displacement_phase","nose_lateral_displacement_phase")

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

	if (type == "circular"){
		require(circular)
		data_per_stride[names(data_per_stride) %in% Phenos.circ] <- lapply(data_per_stride[names(data_per_stride) %in% Phenos.circ], function(x) x*(2*pi))
		data_per_stride[names(data_per_stride) %in% Phenos.circ] <- lapply(data_per_stride[names(data_per_stride) %in% Phenos.circ],function(x) circular::circular(x, type = 'angles', units = 'radians'))	

		#Focus on certain speed bins 
		data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]
		data_per_stride <- data_per_stride[data_per_stride$MouseID %in% names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) >= summary(as.numeric(table(data_per_stride$MouseID)))[2]]),]
		data_per_stride$Strain <- vapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]), FUN.VALUE="1")
		data_per_stride$Strain <- vapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]), FUN.VALUE=character(1))
		data_per_stride$Strain <- as.factor(data_per_stride$Strain)
		data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)

		data_per_animal_mean <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.circ)], by = data_per_stride[c("MouseID")], FUN = function(x) mean(x)[[1]])
		BodyLength <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'BodyLength'][1])
		speed <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) mean(data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'speed']))
		Strain <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'Strain'][1])
		TestDate <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'TestDate'][1])
		TestAge <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'TestAge'][1])
		Sex <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'Sex'][1])
		data_per_animal_mean <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, speed, data_per_animal_mean)
		data_per_animal_mean[,names(data_per_animal_mean) %in% Phenos.circ] <- log(data_per_animal_mean[,names(data_per_animal_mean) %in% Phenos.circ])

		data_per_animal_disp <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.circ)], by = data_per_stride[c("MouseID")], FUN = angular.deviation)
		BodyLength <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'BodyLength'][1])
		speed <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) mean(data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'speed']))
		Strain <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'Strain'][1])
		TestDate <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'TestDate'][1])
		TestAge <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'TestAge'][1])
		Sex <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'Sex'][1])
		data_per_animal_disp <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, speed, data_per_animal_disp)
		data_per_animal_disp[,names(data_per_animal_disp) %in% Phenos.circ] <- log(data_per_animal_disp[,names(data_per_animal_disp) %in% Phenos.circ])

		return(list(data_per_animal_mean, data_per_animal_disp))


	} else {

		#Focus on certain speed bins 
		data_per_stride <- data_per_stride[data_per_stride$bingrpname %in% vapply(c(10,15,20,25,30), function(x) paste0("speed_",x,"_ang_vel_neg20"), FUN.VALUE= character(1)), ]
		data_per_stride <- data_per_stride[, c("Strain","MouseID","Sex","TestAge","TestDate","BodyLength",Phenos.lin)]
		data_per_stride <- data_per_stride[data_per_stride$MouseID %in% names(table(data_per_stride$MouseID)[table(data_per_stride$MouseID) >= summary(as.numeric(table(data_per_stride$MouseID)))[2]]),]
		data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub("<.*>", "", data_per_stride$Strain[x]))
		data_per_stride$Strain <- sapply(seq(nrow(data_per_stride)), function(x) gsub(" ", "", data_per_stride$Strain[x]))
		data_per_stride$Strain <- as.factor(data_per_stride$Strain)
		data_per_stride$MouseID <- droplevels(data_per_stride$MouseID)
		data_per_animal_mean <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = mean)
		Strain <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'Strain'][1])
		TestDate <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'TestDate'][1])
		TestAge <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'TestAge'][1])
		BodyLength <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'BodyLength'][1])
		Sex <- sapply(seq(dim(data_per_animal_mean)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_mean$MouseID[x], 'Sex'][1])
	
		data_per_animal_mean <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, data_per_animal_mean)
		data_per_animal_mean <- data_per_animal_mean[,c("Strain","MouseID","Sex","TestAge","TestDate","BodyLength",Phenos.lin)]

		data_per_animal_disp <- aggregate(x = data_per_stride[,names(data_per_stride) %in% c(Phenos.lin)], by = data_per_stride[c("MouseID")], FUN = var)
		Strain <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'Strain'][1])
		TestDate <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'TestDate'][1])
		TestAge <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'TestAge'][1])
		BodyLength <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'BodyLength'][1])
		Sex <- sapply(seq(dim(data_per_animal_disp)[1]), function(x) data_per_stride[data_per_stride$MouseID == data_per_animal_disp$MouseID[x], 'Sex'][1])
	
		data_per_animal_disp <- cbind(Strain, TestAge, TestDate, Sex, BodyLength, data_per_animal_disp)
		data_per_animal_disp <- data_per_animal_disp[,c("Strain","MouseID","Sex","TestAge","TestDate","BodyLength",Phenos.lin)]

		return(list(data_per_stride, data_per_animal_mean, data_per_animal_disp))


	}
	
	
}