libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname anno "!MCLAB/useful_maize_info/RefGenV4/sasdata";


/*

link geneName to genes and transcripts 

*/

/* import subset */
data fsm_ism_iso;
retain isoform geneID ZMgn ;
format ZMgn $ 44.;
set pacbio.fsm_ism_isoform_zmtr ;
ZMgn = scan(associated_transcript, 1, '_T') ;
var1 = scan(isoform,1, '.') ;
var2 = scan(isoform,2, '.') ; 
var3 = scan(isoform,3, '.') ;  
geneID = compress((var1||'.'||var2));
drop associated_transcript var: ;
run ;

data nic_nnc_iso ;
retain isoform geneID zmgn ;
set pacbio.nic_nnc_isoform_zmgn ;
var1 = scan(isoform,1, '.') ;
var2 = scan(isoform,2, '.') ; 
var3 = scan(isoform,3, '.') ;  
geneID = compress((var1||'.'||var2));
drop var: ;
run ;

/* merge in geneNames */

data anno ;
set anno.B73_V4_geneID_geneName;
rename gene = ZMgn ;
run ;

proc sort data = anno ;
by zmgn ;
proc sort data = fsm_ism_iso ;
by zmgn ;
proc sort data = nic_nnc_iso ;
by zmgn ;
run ;

data pacbio.fsm_ism_wanno ;
merge fsm_ism_iso (in=in1) anno (in=in2) ;
by zmgn  ;
if in1 ;
run;

data pacbio.nic_nnc_wanno ;
merge nic_nnc_iso (in=in1) anno (in=in2) ;
by zmgn  ;
if in1 ;
run;


proc sort data = pacbio.fsm_ism_wanno ;
by geneID ;
proc sort data = pacbio.nic_nnc_wanno ;
by geneID ;
run ;

data fsm_ism_nic_nnc ;
set pacbio.fsm_ism_wanno pacbio.nic_nnc_wanno ;
run ;

proc sort data = fsm_ism_nic_nnc nodups ;
by _all_ ;
run; /* nothing deleted */

proc freq data = fsm_ism_nic_nnc ;
tables isoform / out = cnts ;
run;
data ck ;
set cnts ;
where count ne 1 ;
run;  /* uniq */

data pacbio.fsm_ism_nic_nnc_wanno ;
set fsm_ism_nic_nnc ;
run ;


