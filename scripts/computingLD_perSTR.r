cat (paste0("\n", Sys.Date(),"\n"))
cat (paste0("R version: ",getRversion(),"\n"))

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
suppressPackageStartupMessages(library(broom))


############################################################################################################################################
#	Setting variables
############################################################################################################################################

# create parser object and specify our desired options 
parser <- ArgumentParser()
parser$add_argument("-d", "--directory", help="your working directory")
parser$add_argument("-i", "--input", help="file containing genotypes for a given STR and SNVs located <250 kb of the tested STR. Required: SampleID, STRId,SNVId,STR genotype (avg_repeats) ,SNV genotype (value)")  
args <- parser$parse_args()

# set variables 
cat (paste0("Arguments:\n"))
MYDIR=args$directory; cat (paste0(" \tWorking directory: ", MYDIR, "\n"))
MYINPUT=args$input; cat (paste0(" \tInput: ", MYINPUT, "\n"))
SAVENAME= gsub (".txt",".LD.txt",MYINPUT)

############################################################################################################################################
#	set directory and load file 														  
############################################################################################################################################

cat (paste0("... Loading file \n"))
df <- fread (MYINPUT, sep = "\t", check.names = F, header =T) %>% as.data.frame

############################################################################################################################################
#	computing LD														  
############################################################################################################################################

cat (paste0( "Running associations\n"))
#create an empty file
names = c("STRId","NumberOfSNVs","SNVId","Pair.obs","Rsquare","Slope","LinRegP")
tosave <- data.frame(matrix(ncol = length(names), nrow = 0)); colnames(tosave) <- names
write.table(tosave,SAVENAME, sep ="\t", row.names = F, quote = F)

STR = unique (df$STRId)
n = length(unique (df$SNVId))
regions = unique (df$SNVId)

for (i in 1:length(regions)) {
	# select genotypes for a given SNVs and convert genotypes into alternate allele dosages
	x <- df %>% filter (SNVId == regions[i] ) %>% 
		mutate (value = ifelse (value == "./.",NA,ifelse (
								value =="0/0",0,ifelse (
								value =="0/1",1,ifelse (
								value =="1/1",2,NA)))) )
	
	tobj = try (fit <- summary(lm(avg_repeats ~ value, data=x)))
	if (is(tobj,"try-error")) 
	{
		LinRegP <- NA
		Rsq     <- NA
		Slope   <- NA
		Pair.obs <- nrow (x %>% filter (!is.na(value)) %>% filter (!is.na(avg_repeats )))
	} else {	
		LinRegP  <- scientific (glance(fit)$p.value,digits = 5)
		Rsq      <- round(glance(fit)$r.squared,4)   #ow much variation of the dependent variable is explained by a modelcor
		Slope    <- round(fit$coefficients [,"Estimate"]["value"] ,6)
		Pair.obs <- nrow (x %>% filter (!is.na(value)) %>% filter (!is.na(avg_repeats )))
	}
	
	dat <- c(STR,n,regions[i],Pair.obs,Rsq,Slope,LinRegP)
	names(dat) <- names 
	dat <- as.data.frame(t(dat))									
	write.table(dat,SAVENAME,append = T,quote = F,sep = "\t", row.names = F,col.names = F)
	
	rm (x,LinRegP,Rsq,Slope,Pair.obs,dat)
}	

cat (paste0( "The end\n"))
