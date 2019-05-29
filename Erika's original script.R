

############################################################################
############################################################################
##	 																	  ##
##		Erika's normalization calculation and data-quality checking		  ##
##	 																	  ##
############################################################################
############################################################################


################################################
# 		MODIFY THIS SECTION					   #
################################################

# 1. Put your runfile's location in the quotation marks" Change nothing else in line 15.
#<- remove this comment if using lab basement computer ->  setwd("C:/Users/Sercon_1/Desktop") # Place csv file on desktop
setwd("~/Dropbox Oxford/Dropbox/dropbox AGRICURB/Runfiles/160812")
name <- read.csv("160812.csv", header=F)[, c(1:3, 5, 6,9:11)]

runfile.id <- "160812"

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

rownumber <- (1:nrow(name))
name2 <- name[9:((rownumber[name[,1]=="Drift Corrected"])-2),]
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

########################################################################################
# Create the archive copy
########################################################################################

archive <- data.frame(samplename=data$ID, sampleweight = data$Wt, RS="S", owner="AS", PreparedBy ="AS", Funding="Bogaard ERC", sampletype=0, Taxa="", site="", age.context=0, ugN=data$Nugdc, d15N=data$d15Ndc, ugC=data$Cugdc, d13C=data$d13Cdc, CN=data$CN, labno = paste(runfile.id, data$Ps, sep="-"))
list.of.salanines <- data$Ps[archive$samplename=="SALANINE"]
first.s <- list.of.salanines[1]
first.r <- first.s-8
list.of.ralanines <- c(2, first.r, list.of.salanines)
archive$RS <- "S"
archive$RS[list.of.ralanines] <- "R"


archive$sampletype[archive$samplename=="ALANINE" |archive$samplename=="SALANINE" | archive$samplename==paste(RM1.name) | archive$samplename==paste(RM2.name)] <- "Standard"

archive$sampletype[archive$samplename!="ALANINE" & archive$samplename!="SALANINE" & archive$samplename!=paste(RM1.name) & archive$samplename!=paste(RM2.name)] <- "Plant"

archive$age.context[archive$samplename=="ALANINE" |archive$samplename=="SALANINE" | archive$samplename==paste(RM1.name) | archive$samplename==paste(RM2.name)] <- "Modern"

archive$age.context[archive$samplename!="ALANINE" & archive$samplename!="SALANINE" & archive$samplename!=paste(RM1.name) & archive$samplename!=paste(RM2.name)] <- "Bronze Age"


names(archive) <- c("sample name", "sample weight", "R/S", "Owner", "Prepared by", "Funding", "Sample type", "Taxa", "Site", "age/context", "ug N", "d15N AIR", "ug C", "d13C VPDB", "C/N molar", "Lab number")


run.id <- matrix(nrow=2, ncol=16)
run.id[1, 1:2] <- c("Run number", paste(runfile.id))
run.id[2, 1:2] <- c("User", "Amy Styring")
run.id[1:2, 3:16] <- ""
run.id <- data.frame(run.id)
names(run.id) <- names(archive)
archive <- rbind(archive, run.id)
write.csv(archive, paste("archive", runfile.id, "csv", sep="."))

data2 <- read.csv(paste("archive", runfile.id, "csv", sep="."), header=F)

bottom <- nrow(data2)
run.id <- data2[(bottom-1):bottom, ]
just.data <- data2[1:(bottom-2), ]
final.archive <- rbind(run.id[2:17], just.data[2:17])
write.csv(final.archive, paste("final.archive", runfile.id, "csv", sep="."), row.names=FALSE)



################################################
# 		Calculate normalized N data  		   #
################################################



## Step 1 : Pre-calculate all of the means and stdevs of the standards
normd15N <- rep(0, length(data$d15Ndc))
d15Nsd <- rep(0, length(data$d15Ndc))
data <- data.frame(data, normd15N, d15Nsd)


RM1M <- mean(data$d15Ndc[data$ID==paste(RM1.name)])
RM2M <- mean(data$d15Ndc[data$ID==paste(RM2.name)]) 
RM1Msd <- sd(data$d15Ndc[data$ID==paste(RM1.name)])
RM2Msd <- sd(data$d15Ndc[data$ID==paste(RM2.name)])
alaninesd <- sd(data$d15Ndc[data$ID=="SALANINE" & data$Ps > 10])
mean.alanine <- mean(data$d15Ndc[data$ID=="SALANINE" & data$Ps > 10])

## Step 2: Normalize data and calculate uncertainties
## This is based on Kragten's spreadsheet
for (i in 1:nrow(data)){
  x1 <- data$d15Ndc[i]
  measuredcolumn <- c(RM1T.N, RM1M, RM2T.N, RM2M, x1)
  matrix1 <- cbind(c(RM1Tsd.N, 0, 0, 0, 0), 
                   c(0, RM1Msd, 0, 0, 0), 
                   c(0, 0, RM2Tsd.N, 0, 0), 
                   c(0, 0, 0,RM2Msd, 0),
                   c(0,0,0,0,alaninesd))
  matrix2 <- rbind(rep(RM1T.N, 5), rep(RM1M, 5), rep(RM2T.N, 5), rep(RM2M, 5), rep(x1, 5))
  matrix3 <- matrix1 + matrix2
  raw2true <- function(cv) {cv[1] + (cv[5]-cv[2])*((cv[1]-cv[3])/(cv[2]-cv[4]))}
  dsqr <- function(colno){(raw2true(matrix3[,colno])-raw2true(measuredcolumn))^2}
  normalizedd15N <- raw2true(measuredcolumn)
  finalerror <- sqrt(dsqr(1)+dsqr(2)+dsqr(3)+dsqr(4)+dsqr(5))
  data$d15Nsd[i]<- finalerror
  data$normd15N[i] <- normalizedd15N
}

################################################
# 		Calculate normalized C data  		   #
################################################



## Step 1 : Pre-calculate all of the means and stdevs of the standards
normd13C <- rep(0, length(data$d13Cdc))
d13Csd <- rep(0, length(data$d13Cdc))
data <- data.frame(data, normd13C, d13Csd)

RM1M <- mean(data$d13Cdc[data$ID==paste(RM1.name)]) # Make sure to select the correct RM!
RM2M <- mean(data$d13Cdc[data$ID==paste(RM2.name)]) # Make sure to select the correct RM!
RM1Msd <- sd(data$d13Cdc[data$ID==paste(RM1.name)])
RM2Msd <- sd(data$d13Cdc[data$ID==paste(RM2.name)])
alaninesd <- sd(data$d13Cdc[data$ID=="SALANINE" & data$Ps > 10])
mean.alanine <- mean(data$d13Cdc[data$ID=="SALANINE" & data$Ps > 10])

## Step 2: Normalize data and calculate uncertainties
## This is based on Kragten's spreadsheet
for (i in 1:nrow(data)){
  x1 <- data$d13Cdc[i]
  measuredcolumn <- c(RM1T.C, RM1M, RM2T.C, RM2M, x1)
  matrix1 <- cbind(c(RM1Tsd.C, 0, 0, 0, 0), 
                   c(0, RM1Msd, 0, 0, 0), 
                   c(0, 0, RM2Tsd.C, 0, 0), 
                   c(0, 0, 0,RM2Msd, 0),
                   c(0,0,0,0,alaninesd))
  matrix2 <- rbind(rep(RM1T.C, 5), rep(RM1M, 5), rep(RM2T.C, 5), rep(RM2M, 5), rep(x1, 5))
  matrix3 <- matrix1 + matrix2
  raw2true <- function(cv) {cv[1] + (cv[5]-cv[2])*((cv[1]-cv[3])/(cv[2]-cv[4]))}
  dsqr <- function(colno){(raw2true(matrix3[,colno])-raw2true(measuredcolumn))^2}
  normalizedd13C <- raw2true(measuredcolumn)
  finalerror <- sqrt(dsqr(1)+dsqr(2)+dsqr(3)+dsqr(4)+dsqr(5))
  data$d13Csd[i]<- finalerror
  data$normd13C[i] <- normalizedd13C
}
RM1M.C <- RM1M
RM2M.C <- RM2M
RM1Msd.C <- RM1Msd
RM2Msd.C <- RM2Msd		

####################################################################
# 		Now you have normalized data! BUT CAN YOU ACCEPT IT?	   #
####################################################################


write.csv(data, "Alldata.csv")


#######################################################################################

#########################################################################
#########################################################################
##		Diagnostic plots <- CHECK THESE BEFORE ACCEPTING DATA		   ##
#########################################################################
#########################################################################


####################################################################
# 		How good were the standards today?						   #
####################################################################

RM1 <- data[data$ID==paste(RM1.name),]
RM2 <- data[data$ID==paste(RM2.name),]
alanine <- data[data$ID=="SALANINE" & data$Ps > 10,]

samples <- data[data$ID != paste(RM1.name) & data$ID!= paste(RM2.name) & data$ID!= "SALANINE" & data$ID!= "ALANINE",]

d13C <- expression(paste(delta^{13},"C (\u2030) drift corrected"))
d15N <- expression(paste(delta^{15},"N (\u2030) drift corrected"))


###### For Carbon 
par(mfrow=c(1,3))
par(mar=c(4,5,5,1))
boxplot(RM1$d13Cdc)
par(new=T)
plot(rep(1, length(RM1$d13Cdc)), RM1$d13Cdc, axes=F, xlab=paste(RM1.name, "RM1"), ylab=d13C)

boxplot(RM2$d13Cdc)
par(new=T)
plot(rep(1, length(RM2$d13Cdc)), RM2$d13Cdc, axes=F, xlab=paste(RM2.name, "RM2"), ylab=d13C)
mtext("Drift corrected d13C for Two Reference Materials and S-Alanine", line=3)
mtext("Look at the spread of the measured values of References Materials and Alanines. What is the spread in permil?", cex=0.7, line=2)
mtext("Are the measured deltas of the alanines centred (randomly distributed) around their expected mean, in red?", cex=0.7, line=1)
mtext("Or is there a bias in the measurement? Are any of the alanines wrong?", cex=0.7, line=0)
boxplot(alanine$d13Cdc)
abline(h=-27.11000, col="red")
par(new=T)
plot(rep(1, length(alanine$d13Cdc)), alanine$d13Cdc, axes=F, xlab="S-Alanines", ylab=d13C)
dev.copy2pdf(file="plot1.pdf", encoding="WinAnsi")

###### For Nitrogen

###Assign NITROGEN RMs here, carbon RMs below


par(mfrow=c(1,3))
par(mar=c(4,5,5,1))
boxplot(RM1$d15Ndc)
par(new=T)
plot(rep(1, length(RM1$d15Ndc)), RM1$d15Ndc, axes=F, xlab=paste(RM1.name, "RM1"), ylab=d15N)

boxplot(RM2$d15Ndc)
mtext("Drift corrected d15N for Two Reference Materials and S-Alanine", line=3)
mtext("Look at the spread of the measured values of References Materials and Alanines. What is the spread in permil?", cex=0.7, line=2)
mtext("Are the measured deltas of the alanines centred (randomly distributed) around their expected mean, in red?", cex=0.7, line=1)
mtext("Or is there a bias in the measurement? Are any of the alanines wrong?", cex=0.7, line=0)
par(new=T)
plot(rep(1, length(RM2$d15Ndc)), RM2$d15Ndc, axes=F, xlab=paste(RM2.name, "RM2"), ylab=d15N)

boxplot(alanine$d15Ndc)
abline(h=-1.56, col="red")
par(new=T)
plot(rep(1, length(alanine$d15Ndc)), alanine$d15Ndc, axes=F, xlab="S-Alanines", ylab=d15N)
legend("topright", "true alanine = -1.63", lty=1, col="red")
dev.copy2pdf(file="plot2.pdf", encoding="WinAnsi")

#ID position of all dots

####################################################################
# 		What effect did drift correction have?					   #
####################################################################


par(mfrow=c(2,2))

xtext <- "Position in Run"
d15N <- expression(paste("Normalized ", delta^{15},"N (\u2030)"))
d13C <- expression(paste("Normalized ", delta^{13},"C (\u2030)"))
plot(data$Ps, data$normd13C, xlab=xtext, ylab=d13C, main="Change in d13C through run?")
plot(data$Ps, data$normd15N, xlab=xtext, ylab=d15N, main="Change in d15N through run?" )
d15N <- expression(paste("Raw - drift-corrected ", Delta^{15},"N (\u2030)"))
d13C <- expression(paste("Raw - drift-corrected ", Delta^{13},"C (\u2030)"))
plot(data$Ps, data$d13Cdc-data$d13CR, ylab=d13C, xlab=xtext, main="Effect of drift correction vs position")
plot(data$Ps, data$d15Ndc-data$d15NR, ylab=d15N, xlab=xtext, main="Effect of drift correction vs position")
dev.copy2pdf(file="plot3.pdf", encoding="WinAnsi")




####################################################################
# 		What effect did normalization have?						   #
####################################################################


### For carbon
measuredC <- c(mean(RM1$d13Cdc), mean(RM2$d13Cdc))
trueC <- c(RM1T.C, RM2T.C)
trueCsd <- c(RM1Tsd.C, RM2Tsd.C)
par(mfrow=c(1,1))
plot(measuredC, trueC, type="n", ylab=expression(paste("True ", delta^{13},"C (\u2030)")), 
     xlab=expression(paste("Measured ", delta^{13},"C (\u2030)")), 
     ylim=c(min(trueC-1), max(trueC+1)), 
     xlim=c(min(measuredC-1), max(measuredC+1)))
measuredCsd <- c(sd(RM1$d13Cdc), sd(RM2$d13Cdc))
arrows(measuredC-measuredCsd, trueC, measuredC+measuredCsd, trueC, angle=90, code=3, length=0.01)
arrows(measuredC, trueC-trueCsd, measuredC, trueC+trueCsd, angle=90, code=3, length=0.01)
text(measuredC, trueC, c(paste(RM1.name), paste(RM2.name)), cex=0.6, adj=c(0,1))
slope <- abs(trueC[1]-trueC[2])/abs(measuredC[1]-measuredC[2])
intercept <- RM1T.C-(mean(RM1$d13Cdc)*slope)
abline(a=intercept, b=slope, col="red")
legend("topleft", paste("slope =", round(slope,4), ", intercept =", round(intercept,4)), col="red", lty=1)
par(new=T)
plot(samples$d13Cdc, samples$normd13C, axes=F, ylab="", xlab="", ylim=c(min(trueC-1), max(trueC+1)), 
     xlim=c(min(measuredC-1), max(measuredC+1)))
mtext("Measured versus True (normalized) d13C", line=3)
mtext("How much effect does the normalization regression have? How different is the slope from 1,", line=2, cex=0.7)
mtext("how different is the intercept from zero? Are the error bars for the two reference materials as small", line=1, cex=0.7)
mtext("as they should be? Do your samples lie on the regression line between the two reference materials?", line=0, cex=0.7)
dev.copy2pdf(file="plot4.pdf", encoding="WinAnsi")

### For nitrogen

measuredN <- c(mean(RM1$d15Ndc), mean(RM2$d15Ndc))
trueN <- c(RM1T.N, RM2T.N)
trueNsd <- c(RM1Tsd.N, RM2Tsd.N)
par(mfrow=c(1,1))
plot(measuredN, trueN, type="n", ylab=expression(paste("True ", delta^{15},"N (\u2030)")), 
     xlab=expression(paste("Measured ", delta^{15},"N (\u2030)")), 
     ylim=c(min(trueN-1), max(trueN+1)), 
     xlim=c(min(measuredN-1), max(measuredN+1)))
measuredNsd <- c(sd(RM1$d15Ndc), sd(RM2$d15Ndc))
arrows(measuredN-measuredNsd, trueN, measuredN+measuredNsd, trueN, angle=90, code=3, length=0.01)
arrows(measuredN, trueN-trueNsd, measuredN, trueN+trueNsd, angle=90, code=3, length=0.01)
text(measuredN, trueN, c(paste(RM1.name), paste(RM2.name)), cex=0.6, adj=c(0,1))
slope <- abs(trueN[1]-trueN[2])/abs(measuredN[1]-measuredN[2])
intercept <- RM1T.N-(mean(RM1$d15Ndc)*slope)
abline(a=intercept, b=slope, col="red")
legend("topleft", paste("slope =", round(slope,4), ", intercept =", round(intercept,4)), col="red", lty=1)
par(new=T)
plot(samples$d15Ndc, samples$normd15N, axes=F, ylab="", xlab="", ylim=c(min(trueN-1), max(trueN+1)), 
     xlim=c(min(measuredN-1), max(measuredN+1)))
mtext("Measured versus True (normalized) d15N", line=3)
mtext("How much effect does the normalization regression have? How different is the slope from 1,", line=2, cex=0.7)
mtext("how different is the intercept from zero? Are the error bars for the two reference materials as small", line=1, cex=0.7)
mtext("as they should be? Do your samples lie on the regression line between the two reference materials?", line=0, cex=0.7)
dev.copy2pdf(file="plot5.pdf", encoding="WinAnsi")


####################################################################
# 		How are the C/N ratios and mass of C and N?				   #
####################################################################

par(mfrow=c(2,2))
par(mar=c(4,5,3,1))
CNmin <- if (min(samples$CN) < 2.9) {min(samples$CN)} else { 2.8}
CNmax <- if (max(samples$CN) > 3.6) {max(samples$CN)} else { 3.7}
d15N <- expression(paste("Normalized ", delta^{15},"N (\u2030)"))
d13C <- expression(paste("Normalized ", delta^{13},"C (\u2030)"))
plot(samples$CN, samples$normd13C, xlim=c(CNmin, CNmax), xlab="C/N ratio", ylab=d13C)
abline(v=2.9, col="red")
abline(v=3.6, col="red")
plot(samples$CN, samples$normd15N, xlim=c(CNmin, CNmax), xlab="C/N ratio", ylab=d15N)
abline(v=2.9, col="red")
abline(v=3.6, col="red")
xtext <- expression(paste("Sample weight (", mu,"g)"))
plot(samples$CugR, samples$normd13C, ylab=d13C, xlab=xtext)
plot(samples$NugR, samples$normd15N, ylab=d15N, xlab=xtext)
par(mfrow=c(1,1))
mtext("Do the C/N ratios fall within a good range?", cex=0.7, line=2)
mtext("Within the good range, is there any possibility of a trend between C/N and d13C or d15N?", cex=0.7, line=1)
par(mfrow=c(2,1))
mtext("Similarly, was the weight of C and N measured high enough?", cex=0.7, line=2)
mtext("Is there any biasing effect of weight on isotopic ratio, even for the OK samples?", cex=0.7, line=1)
dev.copy2pdf(file="plot6.pdf", onefile=T, encoding="WinAnsi")

#########################################################################
#########################################################################
##					Remove BAD data									   ##
#########################################################################
#########################################################################

# For example, but you can use your own numbers based on what the plots have shown
#samples <- samples[samples$CugR < 300,]
#samples <- samples[samples$NugR > 60,]
#samples <- samples[samples$CN > 2.9 & samples$CN < 3.6, ]


#########################################################################
#########################################################################
##					Export the accepted samples						   ##
#########################################################################
#########################################################################

write.csv(samples, "Samples.160812.csv")


