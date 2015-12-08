/* Assemble splicing event catalogue */
/* libraries */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';
libname splice2 '/media/jrbnewman/SAS_WRK1/';


/* add in missing variables and rename existing ones */
/* for intron retention events */

data intron_retention_formatting;
    set splice2.intron_retention_events_hg19;
    length event_type $20.;
    length xscripts_cat $3000.;
    length feature1_id $65.;
    length feature1_type $25.;
    length feature2_id $65.;
    length feature2_type $25.;
    length gene_cat $120.;
    event_type="intron_retention";
    rename gene_cat=genes_cat;
    if strand="+" then do; *split intron retention on strand;
        feature1_id=stop_exon;
        feature1_type="exon_donor";
        feature1_start=intron_ret_start;
        feature1_stop=event_stop;
        feature2_id="intron";
        feature2_type="intron_start";
        feature2_start=event_stop+1;
        feature2_stop=intron_ret_stop;
	end;
    else do;
        feature1_id="intron";
        feature1_type="intron_start";
        feature1_start=intron_ret_start;
        feature1_stop=event_start-1;
	feature2_id=start_exon;
        feature2_type="fusion_donor";
        feature2_start=event_start;
        feature2_stop=intron_ret_stop;
        end;
    xscripts_cat="";
    flag_junction_annotated=0;
    flag_firstexon=0;
    flag_altfirstexon=0;
    flag_lastexon=0;
    flag_altlastexon=0;
    flag_exonskip=0;
    flag_alt_donor=0;
    flag_alt_acceptor=0;
    num_xscripts=0;
    num_genes=1;
    flag_multigene=0;
run;

/* for exon-exon junctions */

data exon_junction_formatting;
    length event_type $20.;
    length chr $4.;
    length exonA $65.;
    length exonB $65.;
    length feature1_type $25.;
    length feature2_type $25.;
    length junction_id $100.;
    set splice2.junctions_w_flags_fixed_hg19;
    rename junction_id=event_id;
    event_type="exon_junction";
    rename exonA=feature1_id;
    feature1_type="exon_donor";
    rename donor_start=feature1_start;
    rename donor_stop=feature1_stop;
    rename exonB=feature2_id;
    feature2_type="exon_acceptor";
    rename acceptor_start=feature2_start;
    rename acceptor_stop=feature2_stop;
    flag_intron_retention=0;
    rename flag_exonBskip=flag_exonskip;
run;

/* fixing variable lengths */


data intron_retention_formatting2;
    length genes_cat $120.;
    length xscripts_cat $3000.;
    set intron_retention_formatting;
    format genes_cat $120.;
    informat genes_cat $120.;
run;


data exon_junction_formatting2;
    length feature1_id $65.;
    length feature2_id $65.;
    set exon_junction_formatting;
    format feature1_id $65.;
    format feature2_id $65.;
    format chr $4.;
    informat feature1_id $65.;
    informat feature2_id $65.;
    informat chr $4.;
run;

/* reorder variables */

data intron_retention_formatting3;
    retain
        event_id
        event_type
        genes_cat
        xscripts_cat
        num_genes
        num_xscripts
        flag_multigene
        chr
        feature1_id
        feature1_type
        feature1_start
        feature1_stop
        feature2_id
        feature2_type
        feature2_start
        feature2_stop
        strand
        flag_intron_retention
        flag_junction_annotated
        flag_firstexon
        flag_altfirstexon
        flag_lastexon
        flag_altlastexon
        flag_exonskip
        flag_alt_donor
        flag_alt_acceptor
        ;
    set intron_retention_formatting2;
    drop
        event_start
        event_stop
        intron_ret_stop
        intron_ret_start
        start_exon
        stop_exon
        ;
run;


data exon_junction_formatting3;
    retain
        event_id
        event_type
        genes_cat
        xscripts_cat
        num_genes
        num_xscripts
        flag_multigene
        chr
        feature1_id
        feature1_type
        feature1_start
        feature1_stop
        feature2_id
        feature2_type
        feature2_start
        feature2_stop
        strand
        flag_intron_retention
        flag_junction_annotated
        flag_firstexon
        flag_altfirstexon
        flag_lastexon
        flag_altlastexon
        flag_exonskip
        flag_alt_donor
        flag_alt_acceptor
        ;
    set exon_junction_formatting2;
    drop gene_id;
run;

/* put datasets together */

data splicing_event_catalogue;
    set intron_retention_formatting3 exon_junction_formatting3;
run;

/*
data info_for_bed;
   length score $1.;
   length color $7.;
   set splice.splicing_event_catalogue;
   score='.';
   color='255,0,0';
   block1_start=0;
   other_start=feature1_start;
   other_stop=feature2_stop;
   if flag_intron_retention=1 then do;
      num_blocks=1;
      if strand="+" then block1_length=feature2_stop-feature1_start;
      else block1_length=feature1_start-feature2_stop;
      end;
   else do;
      num_blocks=2;
      block1_length=feature1_stop-feature1_start;
      block2_length=feature2_stop-feature2_start;
      block2_start=feature2_start-feature1_start;
      end;
   keep event_id chr feature1_start feature1_stop feature2_start feature2_stop strand flag_intron_retention score color num_blocks block1_length block2_length block1_start block2_start other_start other_stop;
run;

/* BED12 format 


chrom		chr
totalStart	feature1_start
totalStop	feature2_stop
name		event_id
score		"."
strand		strand
totalStart	feature1_start
totalStop	feature2_stop
color		255,0,0
num		"1" if intron retention, "2" if junction
lengths		block lengths (total if intron_retention, feature1_stop-feature1_start,feature2_stop-feature2_start
starts		if retention then 0, if junction then (0,feature2_start-feature1_start)

*/

data assemble_bed;
    set info_for_bed;
    retain chr;
    retain feature1_start;
    retain feature2_stop;
    retain event_id;
    retain score;
    retain other_start;
    retain other_stop;
    retain color;
    retain num;
    length lengths_cat $10.;
    length starts_cat $10.;
    lengths_cat=catx(',', block1_length, block2_length);
    starts_cat=catx(',', block1_start, block2_start);
    drop flag_intron_retention block1_length block2_length block1_start block2_start;
run;



chrom		chr
totalStart	feature1_start
totalStop	feature2_stop
name		event_id
score		"."
strand		strand
totalStart	feature1_start
totalStop	feature2_stop
color		255,0,0
num		"1" if intron retention, "2" if junction
lengths		block lengths (total if intron_retention, feature1_stop-feature1_start,feature2_stop-feature2_start
starts		if retention then 0, if junction then (0,feature2_start-feature1_start)


run;

proc sort data=splicing_event_catalogue;
    by event_id chr feature1_start feature1_stop feature2_start feature2_stop;
run;



/* Make permenant */

data splice2.splicing_event_catalogue_hg19;
   set splicing_event_catalogue;
run;

proc export data=splice2.splicing_event_catalogue_hg19
	outfile='/home/jrbnewman/McLab/junction_annotations/generated_files/hg19_splicing_catalogue.csv'
	dbms=csv replace;
	run;

/* Clean up */

proc datasets nolist;
     delete intron_retention_formatting exon_junction_formatting
     intron_retention_formatting2 exon_junction_formatting2
     splicing_event_catalogue
     ;
     run;
     quit;









