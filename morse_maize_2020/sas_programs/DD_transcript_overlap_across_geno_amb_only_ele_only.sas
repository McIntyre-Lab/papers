libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";


/* transcripts */
%macro trans (genotype) ;

data iso_&genotype.;
set pacbio.sub_geno_trt_isoform_onCall_tpm0;
keep flag_&genotype._: transcriptID ;
run;

data iso2_&genotype. ;
set iso_&genotype. ;
length &genotype. $25.;
if flag_&genotype._amb_on0 = 0 and flag_&genotype._ele_on0 = 0  then &genotype. = "Not Detected";
else if flag_&genotype._amb_on0 = 0 and flag_&genotype._ele_on0 = 1  then &genotype. = "Elevated Only";
else if flag_&genotype._amb_on0 = 1 and flag_&genotype._ele_on0 = 0  then &genotype. = "Ambient Only";
else if flag_&genotype._amb_on0 = 1 and flag_&genotype._ele_on0 = 1  then &genotype. = "Both Conditions";
else &genotype. = "oops" ;
drop flag_: ;
run ;

proc freq data = iso2_&genotype. ;
tables &genotype. ;
run ;  /* no oops */

proc sort data = iso2_&genotype. ;
by transcriptID ;
run;

data ele_&genotype ;
set iso2_&genotype ;
where &genotype = "Elevated Only" ;
run ;

data amb_&genotype ;
set iso2_&genotype ;
where &genotype = "Ambient Only" ;
run ;

%mend ;

%trans (B73)  ;
%trans (C123)  ;
%trans (Hp301)  ;
%trans (Mo17)  ;
%trans (NC338)  ;

data amb_all ;
merge amb_: ;
by transcriptID ;
run;


data amb_all_cnt ;
set amb_all ;
if find(B73, "Amb") ge 1 then flag_amb_B73 = "B73" ; else flag_amb_B73 = 0 ;
if find(C123, "Amb") ge 1 then flag_amb_C123 = "C123" ; else flag_amb_C123 = 0 ; 
if find(Hp301, "Amb") ge 1 then flag_amb_Hp301 = "Hp301" ; else flag_amb_Hp301 = 0 ; 
if find(Mo17, "Amb") ge 1 then flag_amb_Mo17 = "Mo17" ; else flag_amb_Mo17 = 0 ; 
if find(NC338, "Amb") ge 1 then flag_amb_NC338 = "NC338" ; else flag_amb_NC338 = 0 ; 
pattern_amb = compress(flag_amb_B73||'_'||flag_amb_C123||'_'||flag_amb_Hp301||'_'||flag_amb_Mo17||'_'||flag_amb_NC338) ;

if find(B73, "Amb") ge 1 then flag2_amb_B73 = "B73" ; else flag_amb_B73 = 0 ;
if find(C123, "Amb") ge 1 then flag2_amb_C123 = "C123" ; else flag_amb_C123 = 0 ; 
if find(Hp301, "Amb") ge 1 then flag2_amb_Hp301 = "Hp301" ; else flag_amb_Hp301 = 0 ; 
if find(Mo17, "Amb") ge 1 then flag2_amb_Mo17 = "Mo17" ; else flag_amb_Mo17 = 0 ; 
if find(NC338, "Amb") ge 1 then flag_amb_NC338 = "NC338" ; else flag_amb_NC338 = 0 ; 
*sum_across_geno = sum(flag_amb_B73, flag_amb_C123, flag_amb_Hp301, flag_amb_Mo17,flag_amb_NC338) ; 

keep transcriptID flag_: pattern_amb ;
run ;

data ele_all ;
merge ele_: ;
by transcriptID ;
run ;

data ele_all_cnt ;
set ele_all ;
if find(B73, "Ele") ge 1 then flag_ele_B73 = "B73" ; else flag_ele_B73 = 0 ;
if find(C123, "Ele") ge 1 then flag_ele_C123 = "C123" ; else flag_ele_C123 = 0 ; 
if find(Hp301, "Ele") ge 1 then flag_ele_Hp301 = "Hp301" ; else flag_ele_Hp301 = 0 ; 
if find(Mo17, "Ele") ge 1 then flag_ele_Mo17 = "Mo17" ; else flag_ele_Mo17 = 0 ; 
if find(NC338, "Ele") ge 1 then flag_ele_NC338 = "NC338" ; else flag_ele_NC338 = 0 ; 

if find(B73, "Ele") ge 1 then flag2_ele_B73 = 1 ; else flag2_ele_B73 = 0 ;
if find(C123, "Ele") ge 1 then flag2_ele_C123 = 1 ; else flag2_ele_C123 = 0 ; 
if find(Hp301, "Ele") ge 1 then flag2_ele_Hp301 = 1 ; else flag2_ele_Hp301 = 0 ; 
if find(Mo17, "Ele") ge 1 then flag2_ele_Mo17 = 1 ; else flag2_ele_Mo17 = 0 ; 
if find(NC338, "Ele") ge 1 then flag2_ele_NC338 = 1 ; else flag2_ele_NC338 = 0 ; 


sum_across_geno = sum(flag2_ele_B73, flag2_ele_C123, flag2_ele_Hp301, flag2_ele_Mo17,flag2_ele_NC338) ;
pattern_ele = compress(flag_ele_B73||'_'||flag_ele_C123||'_'||flag_ele_Hp301||'_'||flag_ele_Mo17||'_'||flag_ele_NC338) ; 
*keep transcriptID flag_:  pattern_ele sum_across_geno;
run ;

proc freq data = ele_all_cnt ;
tables sum_across_geno ;
run ;
    /*  
                                            Cumulative
sum_across_geno    Frequency     Percent     Frequency
-------------------------------------------------------
              1        4568       63.09          4568
              2        1796       24.80          6364
              3         645        8.91          7009
              4         194        2.68          7203
              5          38        0.52          7241
*/







ods pdf file = "!MCLAB/maize_ozone_FINAL/pacbio_paper/penultimate_version/draft_figs_tables/DD_transcript_overlap_across_geno_for_amb_only_and_ele_only.pdf" ;

proc freq data = amb_all_cnt noprint ;
tables pattern_amb  / out = cnts_amb_all;
run ;

proc sort data = cnts_amb_all ;
by descending count ;
run ;

title "Diff detected transcript overlap across genotypes for transcripts detected in ambient only" ;
proc print data = cnts_amb_all ; run;

proc freq data = ele_all_cnt noprint;
tables pattern_ele / out = cnts_ele_all;
run ;
proc sort data = cnts_ele_all ;
by descending count ;
run ;

title "Diff detected transcript overlap across genotypes for transcripts detected in ozone only" ;
proc print data = cnts_ele_all ; run;
ods pdf close ;




