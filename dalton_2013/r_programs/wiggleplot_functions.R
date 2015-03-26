#===============================================================================
#
#          FILE:  wiggleplot_functions.R
# 
#         USAGE:  Called from another R script
# 
#   DESCRIPTION:  A set of functions used to create R based wiggle plots. 
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#===============================================================================

read.wiggle <- function(wiggle.file,my.chrom,my.start,my.end){
    con <- file(wiggle.file,"r")
    newline <- readLines(con,n=1) # Initiate header line
    data2 <- data.frame()

    while(length(newline) != 0){
        mydata <- as.vector(strsplit(newline,split=",")[[1]])

        if(mydata[[1]] == my.chrom & my.start <= mydata[[2]] & mydata[[2]]<=my.end){
            print(mydata[[1]])
            data2 <- rbind(data2,matrix(mydata,ncol=3),deparse.level=0)
        }

        newline <- readLines(con,n=1)
    }

    close(con)

    colnames(data2) <- c("chrom","pos","count")
    return(data2)

}

read.fusion.key <- function(fusion.file) {
    # Function for reading in a table containing mapping fusion/gene_id to exon locations
    # 
    # INPUT
    # TSV table need at least: fusion/gene_id, exon_id, chromosome, exon_start, exon_end
    # My need to edit "fusion2gene.colnames" to fit your table
    #
    # OUTPUT
    # a data.frame
    #
    fusion2gene <- read.table(fusion.file, stringsAsFactors=F, skip=1,fill=T);
    fusion2gene.colnames <- c("fusion_id",
                              "exon_id",
                              "chrom",
                              "start",
                              "end",
                              "exon_gene_id",
                              "exon_name",
                              "chromosome",
                              "start_exon",
                              "end_exon",
                              "strand",
                              "fbtr",
                              "fbtrs_per_exon",
                              "fbgn",
                              "fbpp",
                              "symbol",
                              "sequence_loc");
    colnames(fusion2gene) <- fusion2gene.colnames;
    return(fusion2gene)
}

read.track.data <- function(track.file) {
    # Function to read in a table containing track data, could be a simple counts
    # from a pileup or mean counts from SAS output
    #
    # INPUT
    # TSV table with: chromosome, coordinate, mean/count
    #
    # OUTPUT
    # a data.frame
    #
    track.colnames <- c("chrom","coord","mean");
    track.data <- read.table(track.file,stringsAsFactors=F);
    colnames(track.data) <- track.colnames;
    return(track.data);
}

subset.fusion.table <- function(fusion.table, gene){
    # Function that subsets fusion2gene table based upon the given gene_id
    #
    # INPUT
    # fusion2gene
    #
    # OUTPUT
    # Subset of fusion2gene for only one gene_id
    #
    gene.subset <- subset(fusion.table, fusion.table[16] == gene);
    gene.chrom <- gene.subset$chromosome[1];
    gstart.pos  <- min(gene.subset$start_exon);
    gend.pos  <- max(gene.subset$end_exon);
    return(list(gsubset=gene.subset,chrom=gene.chrom,gstart=gstart.pos,gend=gend.pos));
}

track.subset.fun <- function(track,chromosome,gene.start.position,gene.end.position){
    # Function that subsets track/pilupe table based upon the given gene_id
    #
    # INPUT
    # track.data
    #
    # OUTPUT
    # Subset of track.data for only one gene_id
    #
    track.subset <- subset(track,track$chrom==chromosome & track$coord>gene.start.position & track$coord<gene.end.position);
    return(track.subset)
}

track.combine.fun <- function(track.1, track.2,col1.name,col2.name){
    # Function that merges multiple tracks into one data set by coordinate
    # location. Most likely you can pull different data from the separate
    # tracks, but I was worried about how coordinates with information only in
    # one tracks would be handled. This functions sets all NA data points to 0.0
    #
    # INPUT
    # track.1.subset track.2.subset
    #
    # OUTPUT
    # Combined track data for given gene_id subset
    #
    track.1 <- track.1[-1];
    track.2 <- track.2[-1];
    colnames(track.1) <- sub("mean",col1.name,colnames(track.1));
    colnames(track.2) <- sub("mean",col2.name,colnames(track.2));
    track.combined <- merge(x=track.1,y=track.2,by="coord",all=TRUE);
    track.combined[,c(col1.name,col2.name)] <- apply(track.combined[,c(col1.name,col2.name)],2,function(x){replace(x,is.na(x),0)});
    track.combined[,"minmean"] <- apply(track.combined[,c(col1.name,col2.name)],1,min);
    return(track.combined);
}

plot.wiggle <- function(track1,main.title=NULL,color="black"){
    # Function that plots a single wiggle track. Note enclose this function in
    # pdf/png(<OPTIONS>) <plot.wiggles> dev.off() to get saved output.
    #
    # INPUT
    # track.subset
    #
    # OUTPUT
    # wiggle track
    #
    gstart.pos <- min(track1$coord);
    gend.pos <- max(track1$coord);
    plot(x=NULL,
         y=NULL,
         axes=F,
         xlab=NA,
         ylab=NA,
         xlim=c(gstart.pos,gend.pos),
         ylim=c(0,100),
         main=main.title);
    axis(1,pos=0,xaxp=c(gstart.pos,gend.pos,20))
    axis(2,pos=gstart.pos)

    # Plot Reads #
    lines(x=track1$coord,y=track1$count,type="h",col=color);
}

plot.overlay.wiggle <- function(track.combine,col1.name,col2.name,main.title=NULL,color1="red",color2="blue",color3="maroon",gstart.pos=NULL,gend.pos=NULL){
    # Function that plots multiple wiggle tracks in an overlay fashion. Note enclose this function in
    # pdf/png(<OPTIONS>) <plot.wiggles> dev.off() to get saved output.
    #
    # INPUT
    # track.combined.subset
    #
    # OUTPUT
    # wiggle track
    #
    if(is.null(gstart.pos)){gstart.pos <- min(track.combine$coord)-100}
    if(is.null(gend.pos)){gend.pos <- max(track.combine$coord)+100}
    ## PLOT OVERLAY ##
    plot(x=NULL,
         y=NULL,
         axes=F,
         xlab=NA,
         ylab=NA,
         xlim=c(gstart.pos,gend.pos),
         ylim=c(0,ceiling(max(c(max(track.combine[,col1.name]),max(track.combine[,col2.name]))))+20),
         main=main.title);

    axis(1,pos=0,xaxp=c(gstart.pos,gend.pos,20))
    axis(2,pos=gstart.pos)

    # Plot All Mean Reads #
    lines(x=track.combine$coord,y=track.combine[,col1.name],type="h",col=color1);
    lines(x=track.combine$coord,y=track.combine[,col2.name],type="h",col=color2);
    lines(x=track.combine$coord,y=track.combine$minmean,type="h",col=color3);
}

plot.gene.model <- function(gene.subset.lst, uniq.fbtr, num.trans, bar.thickness=3, model.ylim=20, model.color="black",gstart.pos=NULL, gend.pos=NULL){
    # Function that plots gene models for a given gene_id. Note enclose this function in
    # pdf/png(<OPTIONS>) <plot.gene.model> dev.off() to get saved output.
    #
    # INPUT
    # fusion2gene.subset
    #
    # OUTPUT
    # Plot of models
    #
    btop <- -1;
    bbottom <- btop-bar.thickness;

    if(is.null(gstart.pos)){gstart.pos <- min(gene.subset.lst$gsubset$start_exon)}
    if(is.null(gend.pos)){gend.pos <- max(gene.subset.lst$gsubset$end_exon)}
    
    #plot(x=NULL,y=NULL,xlab=NA,ylab=NA,xlim=c(gstart.pos,gend.pos),ylim=c(-(7*num.trans),0),axes=F);
    plot(x=NULL,y=NULL,xlab=NA,ylab=NA,xlim=c(gstart.pos,gend.pos),ylim=c(-model.ylim,0), axes=F);

    # Loop through each transcript and create gene model #
    for (j in uniq.fbtr){
        transcript <- subset(gene.subset.lst$gsubset,gene.subset.lst$gsubset$fbtr == j);
        # Loop through and create each Exon and Intron for that transcript
        i <- 1;
        while(i <= nrow(transcript)) {
            # Create Exon "Thick Box"
            rect(xleft=transcript[i,9],xright=transcript[i,10],ytop=btop,ybottom=bbottom,col=model.color,border=NA);

            # Create Intron "Thin Line"
            if(i+1 <= nrow(transcript)){
                trancenter <- btop - bar.thickness/2;
                sttop <- trancenter+0.1;
                sbottom <- trancenter-0.1;
                rect(xleft=transcript[i,10],xright=transcript[i+1,9],ytop=sttop,ybottom=sbottom,col=model.color,border=NA);

            }

            i <- i + 1;
        }

        # Move down to plot next gene model
        btop <- bbottom - bar.thickness;
        bbottom  <- bbottom - (2* bar.thickness);
    }
}

plot.scale.bar <- function(gene.subset.lst, bar.thickness=3, bar.length=1000, bar.ylim=10, bar.color="black",gstart.pos=NULL, gend.pos=NULL){
    btop <- 3;
    bbottom <- btop-bar.thickness;

    if(is.null(gstart.pos)){gstart.pos <- min(gene.subset.lst$gsubset$start_exon)-100}
    if(is.null(gend.pos)){gend.pos <- max(gene.subset.lst$gsubset$end_exon)+100}
    
    plot(x=NULL,y=NULL,xlab=NA,ylab=NA,xlim=c(gstart.pos,gend.pos),ylim=c(-bar.ylim,0), axes=F);

    # Create Exon "Thick Box"
    rect(xleft=gstart.pos ,xright=gstart.pos + bar.length, ytop=btop, ybottom=bbottom, col=bar.color);

}
