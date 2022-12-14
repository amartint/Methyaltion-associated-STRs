cat (paste0("\n", Sys.Date(),"\n"))
cat (paste0("R version: ",getRversion(),"\n"))


##########################################################################################################################
# Loading require libraries										 			
##########################################################################################################################

cat (paste0("Loading libraries\n"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(argparse))


##########################################################################################################################
#											 			variables 										 			#
##########################################################################################################################

# create parser object and specify our desired options 
parser <- ArgumentParser()
parser$add_argument("-d", "--directory", help="working directory")
parser$add_argument("-i", "--input", help="input file containing methylation data (columns: ProbeID, sample1, sample2, ... sample n)") 
args <- parser$parse_args()
if (length (args)==0) { stop("please provide arguments!!")}

  
cat (paste0("Arguments:\n"))
MYDIR=args$directory; cat(" \tWorking directory): ", MYDIR, "\n")
myfile=args$input; cat(" \tFilename: ", myfile, "\n")
MYDIROUT= paste0(MYDIR, "/qc"); if (!file.exists(MYDIROUT)) { dir.create(MYDIROUT)}
mytheme <- theme_bw () +
        		theme(plot.title = element_text(face="bold",size = 22,colour = "black"),
        		plot.subtitle = element_text(size = 20,colour= "black"),
        		axis.title = element_text (size = 20,face= "bold"),
        		axis.text.x= element_text(vjust = 0.5,hjust =0.5,color = "black", size = 20,angle= 0),
        		axis.text.y= element_text(colour = "black", size = 20,vjust = 0.5,hjust = 1),
        		axis.text.y.right= element_text(colour = "black", size = 20,vjust = 0.5,hjust = 1),
        		axis.line =  element_line(colour="black", size=1, linetype ="solid"),
        		panel.grid.minor = element_blank(),
        		axis.ticks.length=unit(0.1,"inch"),
        		legend.key = element_rect(colour = "white",fill ="white",size = 0.5,linetype='solid'),
        		legend.direction = "horizontal",
        		legend.position ="bottom",
        		legend.text = element_text(size =18, colour = "black", angle = 0),
        		legend.title = element_text(size =20, face="bold",colour = "black", angle = 0),
        		legend.background = element_rect(colour="black"))


##########################################################################################################################
#	custom functions								 		
##########################################################################################################################

densityplot <- function (MYFILE) {
	
	cat (paste0(" ... calculating density copy number across regions on ", myfile,"\n"))
	n <- ncol(MYFILE) -1;
	Nregions <- format(length (unique(MYFILE$ProbeID)),big.mark=",")
	plotname <- paste0(MYDIROUT, "/",gsub(".txt","", myfile),".qc_densityplot.n",n,".png")
	txtname  <- paste0(MYDIROUT, "/",gsub(".txt","", myfile),".qc_densityplot.n",n,".txt")
	YLAB="Density"; 
	XLAB="ylab"
	TITLE=paste0("Mytitle")
	SUBTITLE=paste0("regions: ",Nregions , " across ", n," samples")
	
	# get density peak per sample
	x <- melt (MYFILE) %>% rename (sampleId = variable) %>% filter (!is.na (value)) %>% as.data.frame
	labels <- x %>% group_by (sampleId) %>% summarise( 
					yPos=max(density(value)$y), 
					yPos = round (yPos,2)) %>% 
				ungroup %>% data.frame	%>%
				arrange (desc(yPos)) 
	cat (paste0(" ... saving ... \n"))
	write.table (labels,txtname, sep ="\t", row.names =F, quote =F)
	
	cat (paste0(" ... plotting ... \n"))
	g <- ggplot (x)  +  
		geom_density (aes(x = value, color= sampleId ), show.legend = FALSE)  +  
		labs (xlab = XLAB, ylab = YLAB, title = TITLE, subtitle = SUBTITLE) + 
		scale_y_continuous(breaks = pretty_breaks ()) +
		scale_x_continuous(breaks = pretty_breaks ()) + 
		mytheme 
	
	png (plotname, width = 750, height = 750);  print(g) ; dev.off()
	
	cat (paste0(" The end \n"))

}	
pca <- function(MYFILE) {
	
	cat (paste0(" ... calculating pca across regions on ", myfile,"\n"))
	Nregions=format (length(unique(MYFILE$ProbeID)), big.mark = ",")
	n <- ncol(MYFILE) -1;
	plotname <- paste0(MYDIROUT, "/",gsub(".txt","", myfile),".qc_pca.n",n,".png")
	txtname  <- paste0(MYDIROUT, "/",gsub(".txt","", myfile),".qc_pca.n",n,".txt")
	varname  <- paste0(MYDIROUT, "/",gsub(".txt","", myfile),".qc_pca_variance.n",n,".txt")
	
	rownames(MYFILE) = MYFILE[,1]; 
	temp <- MYFILE[,-1];
	temp <- na.omit(temp) 
	temp  <- t(temp )
	
	cat ("\t\t ... comuputing pcs ... \n")
 	res.pca <- prcomp(temp, scale = FALSE)  
     
 	# get principal components
 	cat ("\t\t ... extracting pcs ... \n")
 	pca = res.pca$x %>% data.frame
 	pca = cbind(sampleID = rownames(pca), round(pca,4)) %>% select (1:50)
	
 	cat ("\t\t ... getting variance explained by each PC... \n") # get and calculate proportion of variance for component
 	std_dev <- res.pca$sdev; rm (res.pca)
 	pr_var <- std_dev^2
 	prop_varex <- round((pr_var/sum(pr_var) *100 ),2)
 	a <- paste0("PC",1:length (prop_varex ));
 	prop_varex <- as.data.frame(cbind("PCs"=a,prop_varex))
 	write.table (prop_varex,varname , sep ="\t", row.names = F,quote = F) 
 	prop_varex = read.table(varname , sep ="\t", check.names = F,header =T)
 	
 	TITLE=paste0("PCA")
 	SUBTITLE=paste0("Regions:  ", format(Nregions,big.mark=","), " regions across ",format (n,big.mark = ",")," samples")
 	
 	pdf ( plotname , width = 8, height = 8)
 	for (j in 1:9) {
          	A <- pca [,grep (paste0("^PC",j,"$"),colnames(pca))]
          	B <- pca [,grep (paste0("^PC",j+1,"$"),colnames(pca))]
          	XLAB = paste0("PC",j," (",prop_varex[j,2],"%)")
          	YLAB = paste0("PC",j+1," (",prop_varex[j+1,2],"%)")
          	p<- ggplot (pca,aes(x=A,y=B)) +  geom_point (size=2,alpha = 0.9) + 
          	 	mytheme +
          	labs (xlab = XLAB, ylab = YLAB , title = TITLE, subtitle = SUBTITLE)
     	 	print (p)
     	 	rm (A,B, XLAB, YLAB,p)
 	}

	TITLE="Variance explained by PCs"
	XLAB="Principal Components";YLAB="Variance explained (%)"
 	
	p<-ggplot (data=prop_varex[1:10,]) +
	     geom_bar (aes(y=prop_varex,x=reorder(PCs,-prop_varex)),stat="identity",color = "black",fill="black") +
	     labs (xlab = XLAB, ylab= YLAB, title = TITLE,subtitle = SUBTITLE)
	     scale_y_continuous (breaks = pretty_breaks(n=10)) + mytheme
	 print (p)
	
	
	dev.off()
	write.table(pca[,1:20],txtname,sep ="\t", row.names = F, quote = F))


}


##########################################################################################################################
#	loading file	and run custom functions								 		
##########################################################################################################################

setwd(MYDIR)
cat (paste0("... loading eth data: ", MYFILE,"\n"))
df <- fread (myfile, sep = "\t", check.names =F, header = T) %>% as.data.frame
densityplot (df)
pca (df)

cat ("The end")
