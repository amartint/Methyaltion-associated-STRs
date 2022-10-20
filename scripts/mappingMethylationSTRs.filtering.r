


cat (paste0("\n", Sys.Date(),"\n"))
cat(paste0("R version: ",getRversion()),"\n")

############################################################################################################################################
#	Loading libraries
############################################################################################################################################

cat (paste0("Loading libraries\n"))
suppressPackageStartupMessages(library (dplyr))
suppressPackageStartupMessages(library (data.table))
suppressPackageStartupMessages(library (tidyverse))
suppressPackageStartupMessages(library (scales))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(bedr))


############################################################################################################################################
#	Setting variables
############################################################################################################################################

# create parser object and specify our desired options 
parser <- ArgumentParser()
parser$add_argument("-d", "--directory", help="your working directory")
parser$add_argument("-p", "--pattern", help="add pattern to generate a list to compile")  
parser$add_argument("-m", "--min", help="minimum number of observations per association")
args <- parser$parse_args()

# set variables 
cat (paste0("Arguments:\n"))
MYDIR=args$directory; cat (paste0(" \tWorking directory): ", MYDIR, "\n"))
MYPATTERN=args$pattern; cat (paste0(" \tPattern: ", MYPATTERN, "\n"))
MYFILTER=as.numeric(args$min); cat (paste0(" \tMin. observations: ", MYFILTER, "\n"))

SAVENAME= paste0("Allregions", gsub (".txt","",MYPATTERN),".raw.txt")
SAVENAME_FILT = paste0("Allregions", gsub (".txt","",MYPATTERN),".filt.txt")

############################################################################################################################################
#	compile resulting files from association														  
############################################################################################################################################

compile_files_per_chrom <- function (mypattern) {

	#mypattern = MYPATTERN
	cat (paste0(" \t ... compiling files: *", mypattern, "\n"))                                                                                                                                                                  

	file.list <- list.files(pattern=mypattern)
	#file.list <- file.list [c(4,5)]
	big.list.of.data.frames <- lapply(file.list, function(x){
		fread(x,header = T,sep ="\t", check.names = F)
	})
	toreturn <- do.call(rbind,big.list.of.data.frames)
	return (toreturn); 
	rm (toreturn,file.list,big.list.of.data.frames )
}

setwd (MYDIR)
df <- compile_files_per_chrom (MYPATTERN)
write.table (df, SAVENAME,sep ="\t", row.names =F, quote = F )

############################################################################################################################################
#	filtering 
############################################################################################################################################


cat (" \t ... filtering for paired observations:\t")
filt <- df %>% filter (Pair.obs >= MYFILTER)

cat (" \t ... adjusting pvalue  ...\n"); 
filt <- filt %>% filter (!is.na (LinRegP) ) %>% mutate (
		FDR = p.adjust (method ="fdr", LinRegP), 
		Bonf = p.adjust (method ="bonferroni", LinRegP)) %>%
		arrange (Bonf)

cat (paste0("Saving file\n"))	
write.table (filt,SAVENAME_FILT, sep ="\t", row.names = F, quote = F)
gzip (SAVENAME_FILT, paste0(SAVENAME_FILT, ".gz"), overwrite = TRUE)	

