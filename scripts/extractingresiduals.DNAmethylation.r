

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

#############################################################
#	Set variables 
#############################################################

# create parser object and specify our desired options 
parser <- ArgumentParser()
parser$add_argument("-d", "--dir", help="your working directory")
parser$add_argument("-s", "--sufix", help="sufix in your file name. File Format (columns:ProbeID, sample1, sample2, sample3, ... sample4)")
parser$add_argument("-q", "--chr", help="chromosome")  
args <- parser$parse_args()

# set variables 
cat (paste0("Arguments:\n"))
MYDIR = args$directory;cat (paste0(" \tWorking directory: ", MYDIR, "\n"))
MYSUFIX = args$directory;cat (paste0(" \tSufix: ", MYSUFIX, "\n"))
MYCHROM = args$chr; cat (paste0(" \tChromosome: ", MYCHROM, "\n"))

DIROUT = paste0(MYDIR, "/Residuals/"); if (!file.exists(DIROUT)){dir.create(DIROUT)} 
MYFILE = paste(MYCHROM,MYSUFIX,"txt", sep = ".")
MYOUTPUT = paste(MYCHROM,MYSUFIX,"residuals.txt", sep = ".")


#############################################################
#	Loading file 
#############################################################

setwd(MYDIR)
cat(paste0 (" ... Loading file containing methyaltion data (ProbeIds in rows and samples in columns) ... \n"))
df <- fread (MYFILE, sep ="\t", check.names = F, header = T) %>% as.data.frame %>%		
			select (ProbeID, <COLUMNS CONTAINING DATA>) %>% 
			tibble::column_to_rownames ("ProbeID") 

cat(paste0 (" ... Loading metadata file ... \n"))
metadata <- fread ("Metadata.txt", sep ="\t", check.names = F, header = T) %>% as.data.frame # file containing metadata, i.e, PCs, age, gender, diagnosis status, ... where each row correspond to an individual sample 	 
rownames (metadata) = metadata$sampleId
#############################################################
#	Extracting residuals
#############################################################

if (file.exists(MYOUTPUT  )) { unlink (MYOUTPUT  ) } 	
for (i in 1:nrow(df)) {
	x<- df[i,]
	CpG <- rownames(x) 
	x <- t(x)
	y = merge (x, metadata , by ="row.names")  # merge by sampleId
	colnames(y)[2] <- "betas"
	tobj = try (fit <- summary(lm(betas ~ cov1 + cov2 + cov3 + ... , data = y)),silent=TRUE)	
	if(is(tobj,"try-error")) 
	{
		cat ("error\t")
		res <- rep (NA,nrow(y))
		names(res) <- y$Row.names
		res <- t(res); res <- cbind (ProbeID = CpG, res)
	} else {
		cat (paste0("getting residuals\t"))
		res = residuals(lm(betas ~ cov1 + cov2 + cov3 + ... , data = y), data = y,na.action= na.exclude))
		names(res) <-y$Row.names
		res <- t(res); res <- cbind (ProbeID= CpG, round(res,6))	
	}
	cat("... Saving data ....\n")
	if (i==1) 
	{ 
		write.table(res,MYOUTPUT , append = T,quote = F,sep = "\t", row.names = F,col.names = T)
		rm (res, CpG,x,y,tobj )
	} else {
		write.table(res,MYOUTPUT, append = T,quote = F,sep = "\t", row.names = F,col.names = F)
		rm (res, CpG,x,y,tobj )
	}
}

cat ("The end\n")
