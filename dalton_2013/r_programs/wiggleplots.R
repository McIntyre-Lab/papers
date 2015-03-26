#===============================================================================
#
#          FILE:  wiggleplots_example.R
# 
#         USAGE:  R CMD wiggleplots_example.R 
# 
#   DESCRIPTION:  This is an example script to create wiggle plots using R.
# 
#  REQUIREMENTS:  If you are going to output as an SVG you need to install
#  "RSvgDevice" by running from inside R:
#
#                 > install.packages("RSvgDevice")
#
#        AUTHOR:  Justin Fear
# 
#===============================================================================

args <- commandArgs(TRUE);
MCLAB <- Sys.getenv("MCLAB");

# set working direcotry to where files are
setwd('/home/jfear/tmp/wiggles')
setwd(args[1])

# Get list of files in directory
files <- list.files()

# load wiggle plot functions
source(file=paste(MCLAB,"/arbeitman_fru_network/r_programs/wiggleplot_functions.R",sep=""));

# Read in table with gene symbol and start/end coordinates. 

symbol2coord.file <- paste(MCLAB,"/useful_dmel_data/flybase530/symbol2coord.csv",sep="");
symbol2coord <- read.csv(symbol2coord.file, stringsAsFactors=F, header=T);

my.gene <- args[2]
my.gene <- 'Sxl'

my.entry <- symbol2coord[symbol2coord$symbol == my.gene,]

my.chrom <- my.entry$chrom
my.start <- my.entry$start
my.end <- my.entry$end


# Read in exon information for building gene mnodels
exon2symbol.file <- paste(MCLAB,"/useful_dmel_data/flybase530/exon2symbol.csv",sep="")
exon2symbol  <- read.csv(exon2symbol.file, stringsAsFactors=F, header=T,colClasses=c(rep("NULL",3),"integer","integer","NULL","character",rep("NULL",3),"character","NULL"))

exon2symbol.subset <- subset(exon2symbol,exon2symbol$symbol == my.gene)

uniq.fbtr <- unique(sort(exon2symbol.subset$FBtr))


## Read in track information. I needed to use mean counts across multiple lanes by using output from:
## /share/mclab/jfear/wiggles/SAS_scripts/wiggles.sas
## Would be simple to edit function for pileups. Follow steps above and 
## track.colnames <- c("chrom","coord","mean");


frua.file <- files[1];
frub.file <- files[2];
fruc.file <- files[3];

#frua.file <- "/home/jfear/tmp/tor_AH_MALE_FRUM_A.csv"
#frub.file <- "/home/jfear/tmp/tor_AH_MALE_FRUM_B.csv"
#fruc.file <- "/home/jfear/tmp/tor_AH_MALE_FRUM_C.csv"

frua <- read.csv(frua.file,header=FALSE);
frub <- read.csv(frub.file,header=FALSE);
fruc <- read.csv(fruc.file,header=FALSE);

colnames(frua) <- c("chrom", "coord", "count")
colnames(frub) <- c("chrom", "coord", "count")
colnames(fruc) <- c("chrom", "coord", "count")

max.y <- max(frua$count,frub$count,fruc$count)

png(paste("/home/jfear/mclab/arbeitman_fru_network/reports/wiggles/",my.gene,".png",sep=""))

par(mfrow=c(2,1),mar=c(3,5,5,5))
plot(x=NULL, y=NULL,axes=T,xlab=NA,ylab=NA,xlim=c(my.start,my.end),ylim=c(0,max.y));
lines(x=frua$coord,y=frua$count,type="h",col="black");
lines(x=frub$coord,y=frub$count,type="h",col="red");
lines(x=fruc$coord,y=fruc$count,type="h",col="blue");

# Plot gene model

par(mar=c(5,5,0,5))
bar.thickness=0.5
model.ylim=(4*length(uniq.fbtr))
model.color="black"

btop <- -1;
bbottom <- btop-bar.thickness;

plot(x=NULL,y=NULL,xlab=NA,ylab=NA,xlim=c(my.start,my.end),ylim=c(-model.ylim,0), axes=F);

# Loop through each transcript and create gene model #
j <- 1
while(j <= length(uniq.fbtr)){
    transcript <- subset(exon2symbol.subset,exon2symbol.subset$FBtr == uniq.fbtr[j]);
    # Loop through and create each Exon and Intron for that transcript
    i <- 1;
    while(i <= nrow(transcript)) {
        # Create Exon "Thick Box"
        rect(xleft=transcript[i,1],xright=transcript[i,2],ytop=btop,ybottom=bbottom,col=model.color,border=NA);

        # Create Intron "Thin Line"
        if(i+1 <= nrow(transcript)){
            trancenter <- btop - bar.thickness/2;
            sttop <- trancenter+0.1;
            sbottom <- trancenter-0.1;
            rect(xleft=transcript[i,2],xright=transcript[i+1,1],ytop=sttop,ybottom=sbottom,col=model.color,border=NA);

        }

        i <- i + 1;
    }

    # Move down to plot next gene model
    btop <- bbottom - 1;
    bbottom <- btop-bar.thickness;
    j <- j + 1;
}

dev.off()
