
ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname event '!MCLAB/event_analysis/sas_data';


/* Convert Ensembl IDs to Aceview IDs
   Since there is no simple Ensembl-to-Aceview ID conversion, I will need to first convert Ensembl IDs to Refseq
   then RefSeq to Aceview */

* Import ENS-to-Refseq;

proc import datafile="!MCLAB/event_analysis/references/ensemb2refseq_hg19.txt" out=ens2refseq
    dbms=tab replace;
guessingrows=67000;
run;

data ens2refseq2;
   set ens2refseq;
   rename gene_stable_id=ens_gene_id entrezgene_id=refseq_gene_id;
run;

proc import datafile="!PATCON/useful_human_data/aceview_hg19/downloaded_files/Aceview2Entrez.txt" out=refseq2av
    dbms=tab replace;
guessingrows=600000;
run;

data refseq2av2;
  set refseq2av;
  if entrez_gene_id__if_any_="" then delete;
  do i=1 by 1 while(scan(entrez_gene_id__if_any_,i,"|") ^= '' );
        refseq_gene_id=scan(entrez_gene_id__if_any_,i,"|")+0;
        output;
        end;
  keep aceview_gene refseq_gene_id;
  rename aceview_gene=gene_id ;
run;

proc sort data=refseq2av2 nodup;
   by refseq_gene_id gene_id;
proc sort data=ens2refseq2 nodup;
  by refseq_gene_id ens_gene_id;
run;

data ens2refseq2av;
   merge ens2refseq2 (in=in1) refseq2av2 (in=in2);
   by refseq_gene_id;
   if in1 and in2;
run;

/* Make this permenant so I can use it later */

data event.hg19_ens2refseq2av;
   set ens2refseq2av;
run;


data ens2av;
 set ens2refseq2av;
 drop refseq_gene_id;
run;

proc sort data=ens2av nodup;
   by ens_gene_id gene_id;
run;

proc freq noprint data=ens2av;
   tables ens_gene_id / out=gene_count;
proc sort data=gene_count;
  by descending count;
proc print data=gene_count(obs=10);
run;
*max=14 aceview genes per ENS max;

data cat_av; 
  array av[14] $38.;
  retain av1-av14;
  set ens2av;
  by ens_gene_id;
  if first.ens_gene_id then do;
     call missing(of av1-av14);
     records = 0;
  end;
  records + 1;
  av[records]=gene_id;
  if last.ens_gene_id then output;
run;



  *clean up the output file;
data cat_av2;
  set cat_av;
  length cat_aceview_id $ 550.;
  cat_aceview_id= catx("|", OF av1-av14);
  keep ens_gene_id cat_aceview_id;
  run;


* Merge in MISO results;

data miso_by_gene;
   set event.t1d_flag_miso_results_by_gene2; 
run;

proc sort data=miso_by_gene;
   by ens_gene_id;
proc sort data=cat_av2;
   by ens_gene_id;
run;

data miso_by_gene_comparable;
   merge cat_av2 (in=in1) miso_by_gene (in=in2);
   by ens_gene_id;
   if in1 and in2;
run;
*23238 ENS genes with Aceview IDs;
*15185 analyzable ENS genes;
*12324 analyzed ENS genes with Aceview IDs;

/* Uncat aceview IDs and summarize to Aceview ID */

data miso_by_gene_av;
   length gene_id $38.;
   set miso_by_gene_comparable;
   do i=1 by 1 while(scan(cat_aceview_id,i,"|") ^= "" );
        gene_id=scan(cat_aceview_id,i,"|");
        output; end;
run;

proc sort data=miso_by_gene_av;
   by gene_id;
proc means data=miso_by_gene_av noprint;
   by gene_id;
   var flag_cd4_on flag_cd8_on flag_Cd19_on flag_all_off flag_any_on flag_miso_testable
       flag_1cell_event flag_2cells_event
       flag_cd4cd8_testable flag_cd4cd19_testable flag_cd8cd19_testable
       flag_sig_cd4cd8_bf5 flag_sig_cd4cd19_bf5 flag_sig_cd8cd19_bf5 flag_sig_bf5
       flag_sig_cd4cd8_bf10 flag_sig_cd4cd19_bf10 flag_sig_cd8cd19_bf10 flag_sig_bf10;
   output out=miso_by_av max=;
run;

/* Merge with T1D splicing results and count overlap */

data t1d_qs;
   set event.t1d_flag_gene_dd_ds_exons;
   if flag_gene_cell_specific=1 then delete;
   if flag_gene_monoexon=0 and flag_gene_exon_dd=1 and flag_cell_by_fus_fdr05=0 then output;
   if flag_gene_monoexon=0 and flag_gene_exon_dd=0 and flag_cell_by_fus_fdr05=1 then output;
   if flag_gene_monoexon=1 and flag_gene_exon_dd=0 and flag_cell_by_fus_fdr05=1 then output;
   if flag_gene_monoexon=0 and flag_gene_exon_dd=1 and flag_cell_by_fus_fdr05=1 then output;
   if flag_gene_monoexon=0 and flag_gene_exon_dd=0 and flag_cell_by_fus_fdr05=0 then output;
   if flag_gene_monoexon=1 and flag_gene_exon_dd=0 and flag_cell_by_fus_fdr05=0 then output;
run;


proc sort data=t1d_qs;
  by gene_id;
proc sort data=miso_by_av;
  by gene_id;
run;

data av_miso_w_qs;
   merge miso_by_av (in=in1) t1d_qs (in=in2);
   by gene_id;
   if in1 then flag_in_miso=1; else flag_in_miso=0;
   if in2 then flag_in_qds=1; else flag_in_qds=0;
   if flag_cd4_on=1 or flag_cd8_on=1 or flag_Cd19_on=1 then flag_gene_any_on=1; else flag_gene_any_on=0;
   if flag_cd4_on=0 or flag_cd8_on=0 or flag_Cd19_on=0 then flag_gene_all_off=1; else flag_gene_all_off=0;
   if flag_cd4cd8_testable=1 or flag_cd4cd19_testable=1 or flag_cd8cd19_testable=1 then flag_miso_gene_testable=1;
   else flag_miso_gene_testable=0;
   if flag_1cell_event=1 or flag_2cells_event=1 then  flag_miso_gene_event_dd=1; else flag_miso_gene_event_dd=0;
   if in2 then output;
run;

/* Now count everything! */

proc freq data=av_miso_w_qs noprint;
   where flag_in_qds=1 and flag_in_miso=1;
   tables flag_cell_by_fus_fdr05*flag_gene_exon_dd*flag_miso_gene_event_dd*
   flag_sig_bf5 / out=gene_count;
run;

proc print data=gene_count;
run;

/*
  flag_cell_     flag_      flag_miso_
    by_fus_      gene_     gene_event_     flag_
     fdr05      exon_dd         dd        sig_bf5    COUNT    PERCENT

       0           0            0            0        473     18.2625
       0           0            1            0        301     11.6216
       0           1            0            0        151      5.8301
       0           1            1            0        129      4.9807
       0           1            1            1          1      0.0386
       1           0            0            0        475     18.3398
       1           0            0            1          3      0.1158
       1           0            1            0        424     16.3707
       1           0            1            1          5      0.1931
       1           1            0            0        285     11.0039
       1           1            0            1          2      0.0772
       1           1            1            0        340     13.1274
       1           1            1            1          1      0.0386


*/


proc freq data=av_miso_w_qs noprint;
      where flag_in_qds=1 and flag_in_miso=1;
   tables  flag_cell_by_fus_fdr05*flag_gene_exon_dd*flag_miso_gene_event_dd*flag_sig_bf10 / out=gene_count;
run;



proc print data=gene_count;
run;

/*
 flag_cell_     flag_      flag_miso_
   by_fus_      gene_     gene_event_      flag_
    fdr05      exon_dd         dd        sig_bf10    COUNT    PERCENT

      0           0            0             0        473     18.2625
      0           0            1             0        301     11.6216
      0           1            0             0        151      5.8301
      0           1            1             0        130      5.0193.

      1           0            0             0        476     18.3784
      1           0            0             1          2      0.0772
      1           0            1             0        429     16.5637

      1           1            0             0        286     11.0425
      1           1            0             1          1      0.0386
      1           1            1             0        341     13.1660

OLD:


                        flag_cell_     flag_       flag_
flag_in_    flag_in_      by_fus_      gene_       miso_       flag_
   qds        miso         fdr05      exon_dd    testable    sig_bf10    COUNT

    1           0            0           0           0           .        601
    1           0            0           1           0           .         79
    1           0            1           0           0           .        384
    1           0            1           1           0           .        114
    1           1            0           0           0           0        255
    1           1            0           0           1           0        519
    1           1            0           1           0           0        119
    1           1            0           1           1           0        162
    1           1            1           0           0           0        239
    1           1            1           0           1           0        666
    1           1            1           0           1           1          2
    1           1            1           1           0           0        261
    1           1            1           1           1           0        366
    1           1            1           1           1           1          1



*/




proc freq data=av_miso_w_qs ;
      where flag_in_qds=1 and flag_in_miso=1;
   tables flag_gene_exon_dd*flag_miso_gene_event_dd
      flag_cell_by_fus_fdr05*flag_sig_bf10 flag_cell_by_fus_fdr05*flag_sig_bf5 ;
run;

proc print data=av_miso_w_qs (keep=gene_id flag_sig_bf5 flag_Cell_by_fus_fdr05);
   where flag_sig_bf5=1 and flag_Cell_by_fus_fdr05=0;
run;

*HK1 gene is not QDS in events, but is DS in MISO;

proc print data=ens2refseq2av (keep=ens_gene_id gene_id);
   where gene_id="HK1";
run;

*ENSG00000156515;

/* Make permenant so I can look through the gene lists more */

data event.t1d_miso2qds_comparison2;
  set av_miso_w_qs;
run;









