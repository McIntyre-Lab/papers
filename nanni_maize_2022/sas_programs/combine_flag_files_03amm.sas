
libname make "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/sasdata";


/*

merge   
        make.flag_detect_cvrg_cnts_shortRead
        make.flag_detect_cvrg_cnts_ccs
        make.onCalls_shrt_gene_tpm5
        make.flag_assembled_transcripts 

        pav and ortho from hoopes:  
            /mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/Hoopes_2018/Zea_mays_Hoopes_2018_orthogroup_para_flag.csv
            /mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/Hoopes_2018/Zea_mays.PAV.freq.txt


*/

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/Hoopes_2018/Zea_mays.PAV.freq.txt"
out = pav 
dbms = tab replace ;
guessingrows = MAX ;
run;

data pav2 ;
set pav ;
rename gene = geneID ;
rename frequency = freq_pav_hoopes ;
flag_pav = 1 ;
run;


proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/gene_lists/Hoopes_2018/Zea_mays_Hoopes_2018_orthogroup_para_flag.csv"
out = para
dbms = csv replace ;
run;

data para2 ;
set para ;
rename gene_id = geneID ;
rename orthogroup = orthogroup_hoopes ;
rename flag_zea_mays_paralog = flag_zea_mays_paralog_hoopes ;
run ;  /* 24458 genes */

proc freq data = para2 ;
tables geneID / out = ck_para2 ;
run ;
data cking ;
set ck_para2 ;
where count ne 1 ;
run;


data assembled ;
set make.flag_assembled_transcripts ;
rename gene = geneID ;
run ;

data cvrg_shrt ;
set make.flag_detect_cvrg_cnts_shortRead ;
rename primary_FBgn = geneID ;
run ;   /* 43,523 genes */

data cvrg_ccs ;
set make.flag_detect_cvrg_cnts_ccs ;
rename primary_FBgn = geneID ;
run ;  /* 43,523 genes */

data ons_shrt ;
set make.onCalls_shrt_gene_tpm5 ;
rename primary_FBgn = geneID ;
run ;


proc sort data = assembled ;
by geneID ;
proc sort data = cvrg_shrt ;
by geneID ;
proc sort data = cvrg_ccs ;
by geneID ;
proc sort data = ons_shrt ;
by geneID ;
proc sort data = pav2 ;
by geneID ;
proc sort data = para2 ;
by geneID ;
run ;

data combine other ;
merge assembled (in=in1) ons_shrt  (in=in2) cvrg_shrt (in=in3) cvrg_ccs  (in=in4) pav2 (in=in5) para2 (in=in6) ;
by geneID ;
if in1 then output combine ;
else output other ;
run;  /* other = 0 */



proc contents data = combine ; run;

/*

flag_detect_shrtRd_geno_gt0
flag_detect_ccs_geno_gt0
flag_assembled_transcript 

*/

data combo2 ;
retain geneID

flag_assembled_fusion_transcript
flag_assembled_novel_transcript
flag_assembled_transcript
flag_b73_reference_gene
flag_detect_ccs_b73_amb_gt0
flag_detect_ccs_b73_ele_gt0
flag_detect_ccs_b73_gt0
flag_detect_ccs_c123_amb_gt0
flag_detect_ccs_c123_ele_gt0
flag_detect_ccs_c123_gt0
flag_detect_ccs_hp301_amb_gt0
flag_detect_ccs_hp301_ele_gt0
flag_detect_ccs_hp301_gt0
flag_detect_ccs_mo17_amb_gt0
flag_detect_ccs_mo17_ele_gt0
flag_detect_ccs_mo17_gt0
flag_detect_ccs_nc338_amb_gt0
flag_detect_ccs_nc338_ele_gt0
flag_detect_ccs_nc338_gt0
flag_detect_shrtRd_b73_amb_gt0
flag_detect_shrtRd_b73_ele_gt0
flag_detect_shrtRd_b73_gt0
flag_detect_shrtRd_c123_amb_gt0
flag_detect_shrtRd_c123_ele_gt0
flag_detect_shrtRd_c123_gt0
flag_detect_shrtRd_hp301_amb_gt0
flag_detect_shrtRd_hp301_ele_gt0
flag_detect_shrtRd_hp301_gt0
flag_detect_shrtRd_mo17_amb_gt0
flag_detect_shrtRd_mo17_ele_gt0
flag_detect_shrtRd_mo17_gt0
flag_detect_shrtRd_nc338_amb_gt0
flag_detect_shrtRd_nc338_ele_gt0
flag_detect_shrtRd_nc338_gt0
flag_in_cvrg_cnts_bed
flag_multigene 
flag_zea_mays_paralog_hoopes 
flag_pav 
freq_pav_hoopes
orthogroup_hoopes
;

set combine ;
if flag_detect_shrtRd_b73_gt0 = . then flag_detect_shrtRd_b73_gt0 = 0 ;
if flag_detect_shrtRd_mo17_gt0 = . then flag_detect_shrtRd_mo17_gt0 = 0 ;
if flag_detect_shrtRd_c123_gt0 = . then flag_detect_shrtRd_c123_gt0 = 0 ;
if flag_detect_shrtRd_hp301_gt0 = . then flag_detect_shrtRd_hp301_gt0 = 0 ;
if flag_detect_shrtRd_nc338_gt0 = . then flag_detect_shrtRd_nc338_gt0 = 0 ;

if flag_detect_ccs_b73_gt0 = .  then flag_detect_ccs_b73_gt0 = 0 ;
if flag_detect_ccs_mo17_gt0 = .  then flag_detect_ccs_mo17_gt0 = 0 ;
if flag_detect_ccs_c123_gt0 = .  then flag_detect_ccs_c123_gt0 = 0 ;
if flag_detect_ccs_hp301_gt0 = .  then flag_detect_ccs_hp301_gt0 = 0 ;
if flag_detect_ccs_nc338_gt0 = .  then flag_detect_ccs_nc338_gt0 = 0 ;

if flag_pav = . then flag_pav = 0;
run ;

data make.combo_flags_shrtRd_ccs ;
set combo2 ;
run;

proc freq data = make.combo_flags_shrtRd_ccs ;
tables flag_detect_ccs_b73_gt0 * flag_detect_shrtRd_b73_gt0 ;
run;  /* files are different - good */

proc export data = make.combo_flags_shrtRd_ccs 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/combination_flag_file_shrtRd_ccs_hoopes.csv" 
dbms = csv replace ;
run;

options orientation=landscape ;
ods pdf file = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/crossTab_genes_ccs_shrtRd_assembled.pdf" ;
%macro freqs (genotype) ;

proc freq data = make.combo_flags_shrtRd_ccs  noprint;
tables flag_b73_reference_gene *  flag_assembled_transcript * flag_in_cvrg_cnts_bed * flag_multigene * flag_detect_ccs_&genotype._gt0 * flag_detect_shrtRd_&genotype._gt0 / out = freq_&genotype.;
run ;

title "crossTab for &genotype." ;
proc print data = freq_&genotype. ;
run;
title "";

%mend ;

%freqs (b73) ;
%freqs (mo17) ;
%freqs (c123) ;
%freqs (hp301) ;
%freqs (nc338)  ;
ods pdf close ;

proc contents data = make.combo_flags_shrtRd_ccs ; run;



    





