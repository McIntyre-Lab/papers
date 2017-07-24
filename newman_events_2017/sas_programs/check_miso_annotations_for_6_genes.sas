ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* 6 MISO genes with exon skipping, but not DD or QDS
   (1) What are they?
   (2) Why no DD or QDS? Is this because of how exon detection is being called?
            (ie, these are cases where there is DD of exons at APN>0, but the exons are ambig with APN>5
   (3) Or are there additional exons not seen in RefSeq annotations? */

/* Get genes */

data gene_check;
  set event.miso_refseq_exonskip_cmpr_dd_qds ;
  where flag_miso_se_diff_bf5=1 and flag_gene_dd=0 and flag_gene_qds=0;
run;

proc print data=gene_check(keep=ens_gene_id gene_id);
run;

/*
  Obs    ens_gene_id           gene_id

   1     ENSMUSG00000005374    27368
   2     ENSMUSG00000031511    54126
   3     ENSMUSG00000038286    68021
   4     ENSMUSG00000022641    70508
   5     ENSMUSG00000038084    74143
   6     ENSMUSG00000024293    77805

Extracting out the set of splicing events tested in miso for these genes:


chr18:10596337:10596439:-@chr18:10595674:10595873:-@chr18:10588437:10588548:-	ENSMUSG00000024293
chr18:10610024:10610351:-@chr18:10596337:10596439:-@chr18:10588437:10588548:-	ENSMUSG00000024293
chr18:10596337:10596439:-@chr18:10588437:10588548:-@chr18:10585976:10586039:-	ENSMUSG00000024293
chr18:10596337:10596439:-@chr18:10593744:10595873:-@chr18:10588437:10588548:-	ENSMUSG00000024293
chr18:10582053:10582184:-@chr18:10577636:10577757:-@chr18:10576999:10577088:-	ENSMUSG00000024293
chr8:11821960:11822049:+@chr8:11824505:11824681:+@chr8:11831493:11831592:+	ENSMUSG00000031511
chr8:11815161:11815300:+@chr8:11817660:11817884:+@chr8:11819638:11819731:+	ENSMUSG00000031511
chr16:29589749:29589834:+@chr16:29597672:29597782:+@chr16:29602298:29602351:+	ENSMUSG00000038084
chr16:29588338:29588445:+@chr16:29588909:29588962:+@chr16:29589767:29589834:+	ENSMUSG00000038084
chr13:34254833:34255118:+@chr13:34257930:34258005:+@chr13:34269850:34270045:+	ENSMUSG00000038286
chr16:50331142:50331215:-@chr16:50320413:50320526:-@chr16:50280482:50280652:-	ENSMUSG00000022641
chr16:50220558:50220706:-@chr16:50209198:50209257:-@chr16:50202491:50202688:-	ENSMUSG00000022641
chr16:50331142:50331215:-@chr16:50280872:50280918:-@chr16:50280482:50280652:-	ENSMUSG00000022641
chr5:136108475:136108788:+@chr5:136110974:136111093:+@chr5:136111865:136115110:+	ENSMUSG00000005374


What exons are these?
1. Are they sig diff in MISO?
2. Do these exons exist in RefSeq?


*/

data miso_bf5;
   set event.miso_bin_diffs_and_bayes;
   where NSC_OLD_bayes_factor ge 5 and NSC_OLD_diff ge abs(0.2);
   keep event_name;
run;

data miso2event;
  set event.miso_se2gene;
  where ens_gene_id in ("ENSMUSG00000005374","ENSMUSG00000031511","ENSMUSG00000038286",
                    "ENSMUSG00000022641","ENSMUSG00000038084","ENSMUSG00000024293");
run;


proc sort data=miso_bf5;
  by event_name;
proc sort data=miso2event;
  by event_name;
run;

data miso_bf5_w_gene;
  merge miso_bf5 (in=in1) miso2event (in=in2);
  by event_name;
  if in1 and in2;
run;


/* Sig MISO:

chr13:34254833:34255118:+@chr13:34257930:34258005:+@chr13:34269850:34270045:+	ENSMUSG00000038286
chr16:29589749:29589834:+@chr16:29597672:29597782:+@chr16:29602298:29602351:+	ENSMUSG00000038084
chr16:50220558:50220706:-@chr16:50209198:50209257:-@chr16:50202491:50202688:-	ENSMUSG00000022641
chr18:10596337:10596439:-@chr18:10593744:10595873:-@chr18:10588437:10588548:-	ENSMUSG00000024293
chr5:136108475:136108788:+@chr5:136110974:136111093:+@chr5:136111865:136115110:+	ENSMUSG00000005374
chr8:11821960:11822049:+@chr8:11824505:11824681:+@chr8:11831493:11831592:+	ENSMUSG00000031511

Exons for these events:

ENSMUSG00000038286:					In RefSeq?
chr13	SE	exon	34162964	34163249	Y
chr13	SE	exon	34166061	34166136	Y
chr13	SE	exon	34177981	34178176	N (shorter version is, end is 34178172)

ENSMUSG00000038084:					In RefSeq?
chr16	SE	exon	29589663	29589748	N (shorter version is, start position is 29589748)
chr16	SE	exon	29597586	29597696	Y
chr16	SE	exon	29602212	29602265	Y


ENSMUSG00000022641:					In RefSeq?
chr16	SE	exon	50202378	50202575	Y
chr16	SE	exon	50209085	50209144	Y
chr16	SE	exon	50220445	50220593	Y


ENSMUSG00000024293:					In RefSeq?
chr18	SE	exon	10588439	10588550	Y
chr18	SE	exon	10593746	10595875	Y
chr18	SE	exon	10596339	10596441	Y


ENSMUSG00000005374:					In RefSeq?
chr5	SE	exon	135632605	135632918	N (longer version is, start position is 135632654)
chr5	SE	exon	135635104	135635223	Y
chr5	SE	exon	135635995	135639240	N (shorter version is, stop position is 135636402)

ENSMUSG00000031511:					In RefSeq?
chr8	SE	exon	11821960	11822049	Y
chr8	SE	exon	11824505	11824681	Y
chr8	SE	exon	11831493	11831592	Y

Okay. Effectively all these exons are in RefSeq (some annotation differences, but shouldn't affect outcome.
Now, for these genes, what does fusion detection look like? Add in start and stop positions so I can match these */

data miso_only;
  set event.miso_refseq_exonskip_cmpr_dd_qds;
  where flag_miso_se_diff_bf5=1 and flag_gene_dd=0 and flag_gene_qds=0;
  keep ens_gene_id gene_id;
run;

data fus2gene;
   set mm10.mm10_refseq_fusion_si_info_v2;
   keep primary_gene_id fusion_id;
   rename primary_gene_id=gene_id;
run;

data fus_info;
   set mm10.mm10_refseq_fusion_si_bed_v2;
   keep chr fusion_start fusion_stop fusion_id;
run;

data fus_dtct;
   set event.flag_fusion_on_apn5;
run;


proc sort data=miso_only;
   by gene_id;
proc sort data=fus2gene nodup;
  by gene_id;
run;

data fus2gene2;
  merge fus2gene (in=in1) miso_only (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=fus2gene2;
  by fusion_id;
proc sort data=fus_info;
  by fusion_id;
proc sort data=fus_dtct;
  by fusion_id;
run;

data fus_check;
  merge fus2gene2 (in=in1) fus_info fus_dtct;
  by fusion_id;
  if in1;
run;

/* Export so I can view and manually check the data */
proc export data=fus_check outfile="!MCLAB/event_analysis/analysis_output/miso_6_genes_to_check.csv"
   dbms=csv replace;
run;

/*


proc sort data=mm10.mm10_refseq_fusion_si_bed_v2;
  by chr fusion_start fusion_stop;
run;

ENSMUSG00000038286:					Detected?
chr13	SE	exon	34162964	34163249	S63891_SI	69666	68021
chr13	SE	exon	34166061	34166136	S63892_SI	69666
chr13	SE	exon	34177981	34178176	S63893_SI	69666
These are all on, but in RefSeq these are annotated to a different gene

ENSMUSG00000038084:					Detected?
chr16	SE	exon	29589663	29589748	Yes
chr16	SE	exon	29597586	29597696	Excluded -- ambiguous detection
chr16	SE	exon	29602212	29602265	Y


ENSMUSG00000022641:					Detected?
chr16	SE	exon	50202378	50202575	Excluded -- ambiguous detection
chr16	SE	exon	50209085	50209144	Excluded -- ambiguous detection
chr16	SE	exon	50220445	50220593	Excluded -- ambiguous detection


ENSMUSG00000024293:					Detected?
chr18	SE	exon	10588439	10588550	Y
chr18	SE	exon	10593746	10595875	Y
chr18	SE	exon	10596339	10596441	Excluded -- ambiguous detection


ENSMUSG00000005374:					Detected?
chr5	SE	exon	135632605	135632918	S197664_SI		215160	27368
chr5	SE	exon	135635104	135635223	S197665_SI		215160
chr5	SE	exon	135635995	135639240	S197667_SI/S197668_SI	215160
These are all on, but in RefSeq these are annotated to a different gene

ENSMUSG00000031511:					Detected?
chr8	SE	exon	11821960	11822049	Ambiguous/Off
chr8	SE	exon	11824505	11824681	N
chr8	SE	exon	11831493	11831592	Ambiguous

Okay so most instances are due to ambiguous detection of exons
