# 
# 
# 1. Put your runfile's location in the quotation marks" 
# sets the working directory to where I have a copy of the test csv file
setwd("~/Documents/work files/R Related")
#  reads the selected columns of the csv file  160812 into a table called name, it ignores the headders.
name <- read.csv("160812.csv", header=F)[, c(1:3, 5, 6,9:11)]
# sets the runfile id to that of the csv file
runfile.id <- "160812"
#
# 2. Assign RM1 and RM2 the same character string as in your spreadsheet.
RM1.name <- "SEAL2" # <- the name of your RM1 exactly as it is in the spreadsheet
RM2.name <- "USGS40" # <- the name of your RM2 exactly as it is in the spreadsheet

# 3. For RM1 and RM2, give the true d15N values relative to AIR and stdev obtained from literature
RM1T.N <- 17.3	  #RM1 d15N mean
RM1Tsd.N <- 0.29		#RM1 d15N standard deviation
RM2T.N <- -4.52 	# RM2 d15N mean
RM2Tsd.N <- 0.06	#RM2 d15N standard deviation

# 4. For RM1 and RM2, give the true d13C values relative to VPDB and stdev obtained from literature
RM1T.C <- -13.3	# RM1 d13C mean
RM1Tsd.C <- 0.11		# RM1 d13C standard deviation
RM2T.C <- -26.39 	# RM2 d13C mean 
RM2Tsd.C <- 0.04	# RM2 d13C standard deviation

################################################
# 	Don't forget to change C RMs below		   #
################################################

#######################################################################################

################################################
#	Formatting columns from raw data file	   #
################################################
# sets up a vector of integers from 1 to the number of rows in the csv file
rownumber <- (1:nrow(name))
# puts the undriftcorrected data into the table name2. starting at row 9 
# and truncating the end statment and its row number at the end of the undrift corrected section
name2 <- name[9:((rownumber[name[,1]=="Drift Corrected"])-2),]

#puts the driftcorrected data from 4 rows beyond the drift corrected statment 
#to the end of the file into a table called name3 
name3 <- name[((rownumber[name[,1]=="Drift Corrected"])+4):nrow(name),]

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
runname <- runfile.id      
data <- data.frame(name2,  name3[,4:8], rep(runname, nrow(name2)), rep(0, nrow(name2)),rep(0, nrow(name2)),rep(0, nrow(name2)))

names(data) <- c("Ps", "ID", "Wt", "NugR", "d15NR", "CugR", "d13CR", "d18OR", "Nugdc", "d15Ndc", "Cugdc", "d13Cdc", "d18Odc", "Runfile", "pcC", "pcN", "CN")



#make numeric things numeric
for (i in c(1, 3:13, 15:17)){data[,i] <- as.numeric(as.character(data[,i]))}

#Add pcC, pcN and CN ratio						#NB these are based on drift corrected values
data$pcC <- data$Cugdc/data$Wt/10
data$pcN <- data$Nugdc/data$Wt/10 
data$CN <- data$Cugdc/data$Nugdc*14/12

data <- data[data$Ps!=7 & data$Ps!=38 & data$Ps!=55, ]