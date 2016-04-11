plot.overlay.wiggle <- function(mydat,main.title=NULL,gstart.pos=NULL,gend.pos=NULL){
    # Function that plots multiple wiggle tracks in an overlay fashion. Note enclose this function in
    # pdf/png(<OPTIONS>) <plot.wiggles> dev.off() to get saved output.
    #
    # INPUT
    # data frame with (chrom, pos, count1, count2, ...)
    #
    # OUTPUT
    # wiggle track
    #
    library(RColorBrewer)
	library(scales)
    
    if(is.null(gstart.pos)){gstart.pos <- min(mydat[2])-100}
    if(is.null(gend.pos)){gend.pos <- max(mydat[2])+100}

    # How many columns are there?
    col.num <- ncol(mydat)

    # Make some colors
    my.color <- brewer.pal(col.num-2, "Set1")
	

	# Get a list of names
	name=names(mydat)[3:col.num]
	
	
    ## PLOT OVERLAY ##
    plot(x=NULL,
         y=NULL,
         axes=F,
         xlab=NA,
         ylab=NA,
         xlim=c(gstart.pos,gend.pos),
         ylim=c(0,ceiling(max(mydat[3:col.num]))+20),
         main=main.title);

    axis(1,pos=0,xaxp=c(gstart.pos,gend.pos,20))
    axis(2,pos=gstart.pos)
	
	
	
    # Plot All Mean Reads #
    for (i in 3:col.num){ 
        lines(x=mydat[,2],y=mydat[,i],type="h",col=alpha(my.color[i-2], 0.25));
       
    }
    legend("topright", name, lwd=1, col=my.color)
    	
}


plot.gene.model <- function(mydat, gene, bar.thickness=3, model.ylim=20, model.color="black",gstart.pos=NULL, gend.pos=NULL){
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


    # Create Subset
    names(mydat) <- c("gene_symbol","chrom","estart","eend","transid")    # renaming incase it is different
    mysub <- subset(mydat,gene_symbol == gene)
    uniq_trans <- unique(mysub[5])

    if(is.null(gstart.pos)){gstart.pos <- min(mysub[3])-100}
    if(is.null(gend.pos)){gend.pos <- max(mysub[4])+100}
    
    plot(x=NULL,y=NULL,xlab=NA,ylab=NA,xlim=c(gstart.pos,gend.pos),ylim=c(-model.ylim,0), axes=F);

    # Loop through each transcript and create gene model #
    for (j in uniq_trans[,1]){
        transcript <- subset(mysub,mysub$transid == j);
	    transcript.sort <- transcript[order(transcript$estart),]

        # Loop through and create each Exon and Intron for that transcript
        i <- 1;
        while(i <= nrow(transcript.sort)) {

            # Create Exon "Thick Box"
            rect(xleft=transcript.sort[i,3],xright=transcript.sort[i,4],ytop=btop,ybottom=bbottom,col=model.color);

            # Create Intron "Thin Line"
            if(i+1 <= nrow(transcript.sort)){
                trancenter <- btop - bar.thickness/2;
                sttop <- trancenter+0.002;
                sbottom <- trancenter-0.002;
                rect(xleft=transcript.sort[i,4],xright=transcript.sort[i+1,3],ytop=sttop,ybottom=sbottom,col=model.color);

            }

            i <- i + 1;
        }

        # Move down to plot next gene model
        btop <- bbottom - bar.thickness;
        bbottom  <- bbottom - (2* bar.thickness);
    }
}

