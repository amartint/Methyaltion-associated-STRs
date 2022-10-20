



cat (paste0("\n", Sys.Date(),"\n"))
cat(paste0("R version: ",getRversion()),"\n")

############################################################################################################################################
#	Loading libraries
############################################################################################################################################

cat (paste0("Loading libraries\n"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(bedr))

#############################################################
#	Set variables 
#############################################################

# create parser object and specify our desired options 
parser <- ArgumentParser()
parser$add_argument("-d", "--dir", help="your working directory")
parser$add_argument("-q", "--chr", help="chromosome")  
parser$add_argument("-s", "--str", help="file containing str genotypes and normalized methylation data for a given STR and CpGs present in its flanks")  
args <- parser$parse_args()

# set variables 
cat (paste0("Arguments:\n"))
MYDIR=args$dir; cat (paste0(" \tWorking directory): ", MYDIR, "\n"))
MYCHROM=args$chr; cat (paste0(" \tchromosome: ", MYCHROM, "\n"))
MYREGION=args$str; cat (paste0(" \tRegion: ", MYREGION, "\n"))


cat (paste0(".... Setting variables\n"))

DIROUT=paste0( MYDIR, "/LinReg/"); if (!file.exists(DIROUT)){dir.create(DIROUT)} 
LoadName <- paste0(MYREGION,".txt")
SaveName <- paste0(DIROUT,MYREGION,".LinReg.txt")


#############################################################
#	Loading Files 
#############################################################

setwd (MYDIR)
if (file.exists(SaveName)) { unlink(SaveName);print.noquote (paste0("Deleting ",SaveName))}
         
print.noquote (paste0("Loading ", LoadName))
df<- fread (LoadName, sep="\t", check.names =F, header = T) # file containing meth and genotypes ("avg_repeats") per sample and CpG:STR. required columns: regionID,sampleID, ,probeID,sampleID, methylation, genotypes)

#############################################################
#	Regression data using residuals
#############################################################				         
cat (paste0( "Running associations\n"))
#create an empty file
names = c("STR","NumberOfProbes","ProbeID","Pair.obs","intersept","Slope_repeatUnits","Rsquared", "LinRegP")
tosave <- data.frame(matrix(ncol = length(names), nrow = 0)); colnames(tosave) <- names
write.table(tosave,SaveName, sep ="\t", row.names = F, quote = F)
					
probes = unique (df$ProbeID)
for (i in 1:length (probes)) {
	cat (paste0(".... ",i, " .... ", probes[i],"\t"))
	x <- df %>% filter (ProbeID == probes[i])

	CpG <- unique(x$ProbeID)
	STR <- unique(x$regionId)
	NumberOfProbes <- length(unique(df$ProbeID))
	
	tobj = try (fit <- summary(lm(value ~ avg_repeats, data=x)),silent = TRUE)
	if(is(tobj,"try-error")) 
	{
		cat("... error ....\t")
		intersept <- NA
		Slope_repeatUnits <- NA
		Rsquared <- NA
		LinRegP <-  NA
		Pair.obs <- nrow (x %>%filter (!is.na(avg_repeats)) %>% filter (!is.na(value)))
		 
	} else {
		cat("... regressing ....\t")
		intersept <-round(coef(fit)[1],6)
		Slope_repeatUnits <- round(fit$coefficients [,"Estimate"]["avg_repeats"] ,4)
		Rsquared <- round(fit$r.squared, digits = 4)
		LinRegP <-  scientific(fit$coefficients[,"Pr(>|t|)"]["avg_repeats"], digits = 5)
		Pair.obs <- attributes(logLik(lm(value ~ avg_repeats, data=x)))$nall
		rm (fit)
	}

    dat <- c(STR,NumberOfProbes,CpG,Pair.obs,intersept, Slope_repeatUnits,Rsquared,LinRegP)
	names(dat) <- names 
	dat <- as.data.frame(t(dat))									
	write.table(dat,SAVENAME,append = T,quote = F,sep = "\t", row.names = F,col.names = F)
	
	rm (x,dat, CpG,STR,NumberOfProbes,intersept, Slope_repeatUnits,Rsquared,LinRegP,Pair.obs,tobj)
} # end of loop fer [probes]

cat ("The end\n")	

