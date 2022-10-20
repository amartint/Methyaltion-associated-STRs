


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


############################################################################################################################################
#	Setting variables
############################################################################################################################################

# create parser object and specify our desired options 
parser <- ArgumentParser()
parser$add_argument("-d", "--directory", help="your working directory")
parser$add_argument("-i", "--input", help="file containing allelic genotypes of STR. Required: SampleID, str allele1 , allele2, ancestry of the sample")  
args <- parser$parse_args()

# set variables 
cat (paste0("Arguments:\n"))
MYDIR=args$directory; cat (paste0(" \tWorking directory): ", MYDIR, "\n"))
MYINPUT=args$input; cat (paste0(" \tInput: ", MYINPUT, "\n"))
DIROUT = paste0(MYDIR, "/vst/"); if (!file.exists(DIROUT)){dir.create(DIROUT)} 
SAVENAME= paste0(DIROUT,gsub (".txt",".vst.txt",MYINPUT))


############################################################################################################################################
#	set directory and load file 														  
############################################################################################################################################

cat (paste0("... Loading file \n"))
df <- fread (paste0(MYINPUT , ".txt"), sep = "\t", check.names = F, header =T)

cat (paste0("... preparing file: arranging str_allele1 and str_allele2 in a single column \n"))
df <- df %>% filter (!is.na ( str_allele1)) %>%  filter (!is.na ( str_allele2)) %>%
		mutate (toremove = paste0( str_allele1,";", str_allele2)) %>% separate_rows (toremove, sep =";",convert =TRUE)  %>% 
		dplyr::rename (repeats = toremove) %>% as.data.frame  # if required: number of copies can be rounded. 


#############################################################
#	selection of columns and rounding alleles
#############################################################
	

cat (paste0( "... computing vst ...\t"))

#create an empty file
names = c("STRId","Superpop","NumberOfsamples","Va","Vt","Vb","Ct","Cb","Vst")
tosave <- data.frame(matrix(ncol = length(names), nrow = 0)); colnames(tosave) <- names
write.table(tosave,SAVENAME, sep ="\t", row.names = F, quote = F)

# define regions and ancestries
regions <- unique (df$STRId)
pop <- unique (df$Superpopulation)


for (i in 1:length (regions)) {
	cat (paste0(i ," ...", regions[i], "\n"))
	y = df %>% filter (STRId ==regions[i]) %>% filter (!is.na (repeats)) %>% as.data.frame
	attach (y)
	Va = round (var(y$repeats, na.rm=FALSE),4) 							  #calculate variance of the alleles across all samples
	for (j in 1:length (pop)) {
		if ( nrow (y %>% filter (Superpopulation==pop[j])) >= 100) 			# filter for samples
		{
			n  = length (y[Superpopulation==pop[j],]$repeats)
			Vt =  var(y[Superpopulation == pop[j],]$repeats, na.rm=FALSE)  # variance for each STR across target population
			Vb = var(y[Superpopulation != pop[j],]$repeats, na.rm=FALSE)   # variance for each STR across background population ( whole set population - background)
			Ct = round(n/N_alleles, 4)										# fraction target/whole set of population
			Cb = round(1 - Ct,4)
			VST =round((Va -((Ct *Vt )+(Cb * Vb )))/Va,4)
		} else {
			cat (paste0("\tno enough samples\n"))
			n = NA
			VST = NA
			Vt = NA
			Vb = NA
			Ct = NA
			Cb = NA
		}
		dat = c(regions[i],as.character(pop[j]),n,Va,Vt,Vb,Ct,Cb,VST)
		names(dat) <- names 
		dat <- as.data.frame(t(dat))									
		write.table(dat,SAVENAME,append = T,quote = F,sep = "\t", row.names = F,col.names = F)
		rm (n, VST,Vt,Vb,Ct,Cb, dat)	
	}
	detach(y)
	rm (Va,y)
}

cat ("The end\n")
