

libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;

/*
https://wviechtb.github.io/metafor/reference/rma.mv.html



prep pdiff t values (=effect size) for meta analysis and create design file
metafor R:   
   output = rma(yi, vi, measure="MD", method="FE", data=datain)
        yi = ES estimate (tvalue from pdiff)
        vi = estimated sampling variance (stdErr from pdiff? 
        
        CD4 female = 39+15 = 54
        CD4 male = 46+13 = 59
        

input:
seq.de_pd_frag_&cell._FbyCase

output:

    group               ES          EV
    CD4_F_CaseVCon      tvalue      stderr
    CD4_M_CaseVCon      tvalue      stderr      
    CD8_F_CaseVCon      tvalue      stderr
    CD8_M_CaseVCon      tvalue      stderr      
    

*/



/* file linking frag to gene  */
proc import datafile = "/TB14/TB14/t1d_case_control_cellType/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv"
out = frag_anno 
dbms = csv replace ;
guessingrows = MAX;
run;
proc contents data = frag_anno ; run;

data frag_anno2 ;
set frag_anno ;
keep ef_id gene_id ;
rename gene_id = geneID ;
rename ef_id = featureID ;
run ; 

proc sort data = frag_anno2 nodups ;
by _all_ ;
proc sort data = frag_anno2 nodups ;
by geneID ;
run;

/* t1d candidate gene list */
proc import datafile = "!MCLAB/t1d_case_control_cellType/candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt"
out = list 
dbms = tab replace ;
run;  /* 147 obs */

data list2 ;
length ensembl_geneID $235.;
set list ;
rename ensembl_geneID = geneID ;
run;

proc sort data = list2 ;
by geneID ;
run;

data frag2gene ;
merge list2 (in=in1) frag_anno2 (in=in2) ;
by geneID ;
if in1 and in2 ;
run;


%macro prepping (cell) ;

/* in DE pdiff output, only keep rows for female caseVScontrol and male caseVScontrol */
data pd_&cell. ;
set seq.de_pd_frag_&cell._FbyCase ;
if flag_female = 0 and flag_case = 0 and _flag_female = 0 and _flag_case = 1 then group = 'M_CaseVCon' ; 
else if flag_female = 1 and flag_case = 0 and _flag_female = 1 and _flag_case = 1 then group = 'F_CaseVCon' ; 
else group = 'A' ;
run ;

data pd2_&cell. ;
set pd_&cell. ;
if find(group, 'M_') ge 1 or find(group, 'F_') ge 1 ;
run ;

proc sort data = pd2_&cell. ;
by featureID ;
run;

/* Effect size */
data ES_&cell. ;
retain featureID group tvalue ;
set pd2_&cell. ;
label tvalue = "effect_size" ;
keep featureID tvalue group ;
rename tvalue = effect_size ;
run;

/* variance */
data SV_&cell. ;
retain featureID group variance;
set pd2_&cell. ;
variance = StdErr ;
keep featureID variance group ;
run ;

proc sort data = ES_&cell. ;
by featureID group ;
proc sort data = SV_&cell. ;
by featureID group ;
run;

/* merge ES and SV */
data ES_SV_&cell. ;
merge ES_&cell. SV_&cell. ;
by featureID group ;
run ;

data ES_SV_&cell. ;
retain geneID featureID group sex ;
set  ES_SV_&cell. ;
if find(group, "F_") ge 1 then sex = "F" ;
else if find(group, "M_") ge 1 then sex = "M" ;
else sex = "A" ;
run;

/* merge in t1d gene list by featureID  */
proc sort data = frag2gene ;
by featureID ;
proc sort data = ES_SV_&cell.;
by featureID ;
run ;

data ES_SV_&cell._t1d ;
merge frag2gene (in=in1) ES_SV_&cell. (in=in2) ;
by featureID ;
if in1 and in2 then output ES_SV_&cell._t1d ;
run ;

data ES_SV_&cell._t1d_4_ma ;
retain geneID featureID group sex cellType effect_size variance ;
set  ES_SV_&cell._t1d ;
cellType = "&cell.";
run ;


proc export data = ES_SV_&cell._t1d_4_ma 
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/ES_SV_&cell._t1d_pdiff_4_ma.csv"
dbms = csv replace ;
run;

/*
 make testsets 
data test_1_gene_&cell.  ;
set ES_SV_&cell._t1d_4_ma ;
where geneID = "ENSG00000003056" ;
run ;

proc export data = test_1_gene_&cell. 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs/test_1_gene_&cell..csv"
dbms = csv replace ;
run;

data test_2_gene_&cell.  ;
set ES_SV_&cell._t1d_4_ma ;
where geneID = geneID = "ENSG00000003056" or geneID = "ENSG00000013725" ;
run ;

proc export data = test_2_gene_&cell. 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs/test_2_gene_&cell..csv"
dbms = csv replace ;
run;
*/
%mend ;

%prepping (CD4);
%prepping (CD8);
%prepping (CD19);


/* combine cell types into single dataset for meta-analysis with cellType as moderator */
data both ;
set  ES_SV_CD4_t1d_4_ma ES_SV_CD8_t1d_4_ma ;
run ;

proc sort data = both ;
by geneID featureID group sex cellType ;
run;

/* check each gene has both cell types - dropping ENSG00000117560, only in CD8 */
proc freq data = both ;
by geneID ;
tables cellType / out = aaa ;
run;
proc transpose data = aaa out = aaa_flip ;
by geneID ;
var count ;
id cellType ;
run;

data aa_drop ;
set aaa_flip ;
if CD4 = . or CD8 = . then flag_drop_gene = 1 ;
else flag_drop_gene = 0 ;
keep geneID flag_drop_gene ;
run ;

proc sort data = aa_drop ;
by geneID ;
proc sort data = both ;
by geneID ;
run ;

data ES_SV_CD4_CD8 ;
merge aa_drop both ;
by geneID ;
run ;

data ES_SV_CD4_CD8_t1d ;
retain geneID exonID featureID group sex cellType cellType_sex ;
set ES_SV_CD4_CD8 ;
var1 = scan(featureID, 1, ":");
var2 = scan(featureID, 2, ":");
exonID = compress(var1||":"||var2);
where flag_drop_gene = 0;
drop flag_drop_gene var1 var2 ;
cellType_sex = compress(cellType||'_'||sex);
run;

proc export data = ES_SV_CD4_CD8_t1d 
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/ES_SV_CD4_CD8_t1d_pdiff_4_ma.csv"
dbms = csv replace ;
run;

proc freq data = ES_SV_CD4_CD8_t1d ;
by exonID ;
tables cellType_sex / out = aa_fet_cnts ;
run;
proc freq data = aa_fet_cnts ;
tables exonID / out = aa_mor_cnts ;
run;


data contrast ;
set aa_fet_cnts ;
if count > 2 then flag_run = 1 ;
else flag_run = 0 ;
keep exonID flag_run ;
run ;

proc sort data = ES_SV_CD4_CD8_t1d ;
by exonID ;
proc sort data = contrast ;
by exonID ;
run ;

data added ;
merge contrast (in=in1) ES_SV_CD4_CD8_t1d (in=in2) ;
by exonID ;
run ;

data ES_SV_CD4_CD8_t1d_exon ;
set added ;
where flag_run = 1 ;
drop flag_run ;
run;

proc export data = ES_SV_CD4_CD8_t1d_exon 
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/ES_SV_CD4_CD8_t1d_pdiff_exon.csv"
dbms = csv replace ;
run;

proc freq data = aa_fet_cnts ;
tables exonID / out = aa_mre_cnts ;
run;

/* all 3 */

/* combine cell types into single dataset for meta-analysis with cellType as moderator */
data three ;
set ES_SV_CD19_t1d_4_ma ES_SV_CD4_t1d_4_ma ES_SV_CD8_t1d_4_ma ;
run ;

proc sort data = three ;
by geneID featureID group sex cellType ;
run;

/* check each gene has both cell types - dropping ENSG00000117560, only in CD8 */
proc freq data = three ;
by geneID ;
tables cellType / out = aaa ;
run;
proc transpose data = aaa out = aaa_flip ;
by geneID ;
var count ;
id cellType ;
run;

data aa_drop ;
set aaa_flip ;
if CD4 = . or CD8 = . or CD19 = . then flag_drop_gene = 1 ;
else flag_drop_gene = 0 ;
keep geneID flag_drop_gene ;
run ;

proc sort data = aa_drop ;
by geneID ;
proc sort data = three ;
by geneID ;
run ;

data ES_SV_CD4_CD8 ;
merge aa_drop three ;
by geneID ;
run ;

data ES_SV_CD4_CD8_CD19_t1d ;
retain geneID exonID featureID group sex cellType cellType_sex ;
set ES_SV_CD4_CD8 ;
var1 = scan(featureID, 1, ":");
var2 = scan(featureID, 2, ":");
exonID = compress(var1||":"||var2);
where flag_drop_gene = 0;
drop flag_drop_gene var1 var2 ;
cellType_sex = compress(cellType||'_'||sex);
run;

proc export data = ES_SV_CD4_CD8_CD19_t1d 
outfile = "!MCLAB/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs/ES_SV_CD4_CD8_CD19_t1d_pdiff.csv"
dbms = csv replace ;
run;

proc freq data = ES_SV_CD4_CD8_CD19_t1d ;
by exonID ;
tables cellType_sex / out = aa_fet_cnts ;
run;
proc freq data = aa_fet_cnts ;
tables exonID / out = aa_mor_cnts ;
run;


data contrast ;
set aa_fet_cnts ;
if count > 2 then flag_run = 1 ;
else flag_run = 0 ;
keep exonID flag_run ;
run ;

proc sort data = ES_SV_CD4_CD8_CD19_t1d ;
by exonID ;
proc sort data = contrast ;
by exonID ;
run ;

data added ;
merge contrast (in=in1) ES_SV_CD4_CD8_CD19_t1d (in=in2) ;
by exonID ;
run ;

data ES_SV_CD4_CD8_CD19_t1d_exon ;
set added ;
where flag_run = 1 ;
drop flag_run ;
run;

proc export data = ES_SV_CD4_CD8_CD19_t1d_exon 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs/ES_SV_CD4_CD8_t1d_pdiff_exon.csv"
dbms = csv replace ;
run;

proc freq data = aa_fet_cnts ;
tables exonID / out = aa_mre_cnts ;
run;

