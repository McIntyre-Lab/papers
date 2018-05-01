ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Subset Cav1 RefSeq transcripts, and how much of their events are detected

RefSeq		PacBio2
NM_001243064			
NM_007616   	PB.5785.1 (FSM)
XM_006504974	PB.5785.2 (FSM)
                PB.5785.3 (NIC)
XM_006504975	PB.5785.4 (FSM)
                PB.5785.5 (NIC)
XM_006504976	PB.5785.6 (FSM)

*/

data xs2gene;
   set event.feature2xs2gene_exp_only_nomulti;
   where gene_id="12389";
   keep gene_id transcript_id;
run;

proc sort data=xs2gene nodup;
  by transcript_id;
run;

data xs_info;
  set event.xscripts_w_unique_by_bin;
run;

proc sort data=xs_info;
   by transcript_id;
run;

data cav1_xs_info;
  merge xs2gene (in=in1) xs_info (in=in2);
  by transcript_id;
  if in1 and in2;
run;

ods listing;
proc print data=cav1_xs_info(keep=transcript_id perc_features_dtct perc_unique_features_dtct);
run;

/*
                             perc_      perc_unique_
           transcript_     features_      features_
    Obs    id                 dtct          dtct

     1     NM_001243064        1              1
     2     NM_007616           1              .
     3     XM_006504974        1              1
     4     XM_006504975        1              1
     5     XM_006504976        1              1
 So, we should see all of these, with the possible exception of NM_007616 (it has no unique pieces so we can't really tell for certain if it is expressed).

*/

/* Since NM_001243064 is eliminated in PacBio, let's look at the coverage of its unique fragments */

data frag2xs;
   set mm10.mm10_exon_fragment_flagged;
   where gene_id ? "12389";
   keep fragment_id transcript_id;
run;


data frag_cov;
   set event.mm10_refseq_fragment_counts;
   where sample_id ? "NSC";
   keep fragment_id apn;
run;

proc sort data=frag2xs nodup;
   by fragment_id;
proc sort data=frag_cov;
   by fragment_id;
run;

data frag_check;
  merge frag_cov (in=in1)  frag2xs (in=in2);
  by fragment_id;
  if in1 and in2;
run;

proc sort data=frag_check;
   by event_id transcript_id;
proc means data=frag_check noprint;
  by fragment_id transcript_id;
  var apn;
  output out=mean_apn_per_frag mean=;
run;

proc print data=mean_apn_per_frag(keep = fragment_id transcript_id apn);
run;

/*
 fragment_id     transcript_id                                                            apn

 F201665_SI:1    NM_007616                                                       0.4166666667
 F201665_SI:2    NM_007616|XM_006504974                                          32.712765957
 S201666_SI:1    XM_006504975                                                    0.4203539823
 S201667_SI:1    XM_006504976                                                     0.472972973
 F201668_SI:1    NM_001243064                                                    16.586776859
 F201668_SI:2    NM_001243064|NM_007616|XM_006504974|XM_006504975|XM_006504976   145.22891566
 F201669_SI:1    NM_001243064|NM_007616|XM_006504974|XM_006504975|XM_006504976   94.406593406
 F201669_SI:2    NM_001243064|NM_007616|XM_006504975|XM_006504976                53.108346709
 F201669_SI:3    NM_001243064|NM_007616|XM_006504974|XM_006504975|XM_006504976   66.758116883

F201668_SI:1 is the only unique piece for NM_001243064 (which is not in PB)
Data would suggest that XM_006504974 should be the MET here

*/


*all fragments are on;

/* check coverage of fusions */

data event2xs;
  set event.feature2xs2gene_exp_only_nomulti;
   where gene_id="12389";
run;


data fus;
   set event2xs;
   where feature_id ? "_SI";
   if count(feature_id,":") > 0 then delete;
   keep feature_id;
run;

data coverage_check;
   set event.mm10_refseq_fusion_counts;
   where sample_id ? "NSC";
   keep fusion_id apn;
   rename fusion_id=feature_id;
run;

proc sort data=fus nodup;
   by feature_id;
proc sort data=coverage_check;
   by feature_id;
run;

data fus_check;
  merge coverage_check (in=in1)  fus(in=in2);
  by feature_id;
  if in1 and in2;
run;

proc means data=fus_check noprint;
  by feature_id;
  var apn;
  output out=mean_apn_per_fus mean=;
run;

data fus2xs;
   set mm10.mm10_si_fusions_unique_flagged;
   keep fusion_id transcript_id;
   rename fusion_id=feature_id;
run;

proc sort data=fus2xs nodup;
   by feature_id;
proc sort data=mean_apn_per_fus;
   by feature_id;
run;

data fus_apn_w_xs;
  merge mean_apn_per_fus (in=in1) fus2xs (in=in2);
   by feature_id;
   if in1 and in2;
run;



proc print data=fus_apn_w_xs(keep=feature_id transcript_id apn);
run;

/*
 Obs    feature_id                apn    transcript_id

Obs   feature_id               apn   transcript_id

 1    F201665_SI      31.613013699   NM_007616|XM_006504974
 2    F201668_SI      69.095823096   NM_001243064|NM_007616|XM_006504974|XM_006504975|XM_006504976
 3    F201669_SI      67.558186739   NM_001243064|NM_007616|XM_006504974|XM_006504975|XM_006504976
 4    S201666_SI      0.4203539823   XM_006504975
 5    S201667_SI       0.472972973   XM_006504976


*/

/* Now for these transcripts, get the mean TPM from RSEM, from:
   all RefSeq transcripts
   all EA transcripts
   EA transcripts with 100% events detected (APN>0)
   all pacbio transcripts
*/

data tpm_all;
   set event.rsem_refseq_all;
   where transcript_id in ("NM_001243064","NM_007616","XM_006504974","XM_006504975","XM_006504976");
   mean_TPM_refseq_all=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_TPM_refseq_all;
run;

data tpm_any;
   set event.rsem_events_exp_any;
   where transcript_id in ("NM_001243064","NM_007616","XM_006504974","XM_006504975","XM_006504976");
   mean_TPM_event_all=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_TPM_event_all;
run;

data tpm_100;
   set event.rsem_events_exp_100perc;
   where transcript_id in ("NM_001243064","NM_007616","XM_006504974","XM_006504975","XM_006504976");
   mean_TPM_event_100=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_TPM_event_100;
run;

data tpm_pb;
   length transcript_id $20.;
   set event.rsem_pacbio_all;
   where transcript_id ? "PB.5785.";
   if transcript_id="PB.5785.1" then transcript_id="NM_007616";
   if transcript_id="PB.5785.2" then transcript_id="XM_006504974";
   if transcript_id="PB.5785.3" then transcript_id="Novel_NNC";
   if transcript_id="PB.5785.4" then transcript_id="XM_006504975";
   if transcript_id="PB.5785.5" then transcript_id="Novel_NIC";
   if transcript_id="PB.5785.6" then transcript_id="XM_006504976";
   mean_TPM_pacbio_all=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_TPM_pacbio_all;
run;


proc sort data=tpm_all;
   by transcript_id;
proc sort data=tpm_any;
   by transcript_id;
proc sort data=tpm_100;
   by transcript_id;
proc sort data=tpm_pb;
   by transcript_id;
run;

data all_mean_tpm;
  merge tpm_all tpm_any tpm_100 tpm_pb;
  by transcript_id;
run;

proc print data=all_mean_tpm;
run;


/*
       transcript_      mean_TPM_    mean_TPM_    mean_TPM_     mean_TPM_
Obs    id              refseq_all    event_all    event_100    pacbio_all

 1     NM_001243064      10.600        12.64        12.930         .
 2     NM_007616          8.390        10.53        11.125       29.765
 3     Novel_NIC           .             .            .          14.260
 4     Novel_NNC           .             .            .           0.865
 5     XM_006504974      13.970        16.84        18.090        4.320
 6     XM_006504975       0.415         0.49         0.540        0.340
 7     XM_006504976       0.345         0.45         0.545        0.215


Okay, so NM_001243064 and XM_006504974 could be syphoning off reads from
    Novel_NIC coverage and NM_007616, so now to check the sequence/structure
   of these transcripts 
*/




/* 

For Cav1:

Cav1 has 5 Refseq transcripts, of which 4 have a matching PacBio transcript. There also two novel PB transcripts:

RefSeq		PacBio?	Match	MET in PB
NM_001243064	No		
NM_007616   	Yes	FSM	Yes
XM_006504974	Yes	FSM
XM_006504975	Yes	FSM
XM_006504976	Yes	FSM
Novel_NIC	Novel	NIC
Novel_NNC	Novel	NNC			

Looking at the proportion of events detected. All RefSeq transcripts have all their events detected.
As NM_001243064 was not found in the PacBio data, I checked to see the quantification of the unique pieces of NM_001243064. This transcript has a unique fragment at the 5'end of its first exon, however it is clearly present in the data (APN is >5 in both reps).

Looking at RSEM estimates shows that the MET for Event transcripts is XM_006504974, but for PacBio transcripts it is NM_007616:

Event Analysis transcripts:
NM_001243064: mean TPM=12.930
NM_007616: mean TPM=11.125
XM_006504974: mean TPM=18.090
XM_006504975: mean TPM=0.540
XM_006504976: mean TPM=0.545

PacBio transcripts:
NM_007616: mean TPM=29.765
XM_006504974: mean TPM=4.320
XM_006504975: mean TPM=0.340
XM_006504976: mean TPM=0.215
Novel_NIC: mean TPM=14.260
Novel_NNC: mean TPM=0.865




Here, the MET for Event transcripts is XM_006504974, but for PacBio transcripts it is NM_007616. I suspect NM_001243064 and XM_006504974 could be siphoning off reads from NM_007616 and Novel_NIC, possibly  so I then decided to check the exon structure of these transcripts to check this, specifically where these transcripts are 

NM_001243064 comprised two exons, of which only the 5'end of the first is unique to this transcript (its junction and remaining exonic sequence is shared by other RefSeq transcripts. Part of the 5'-end of NNC overlaps with NM_007616



RefSeq		PacBio2
NM_001243064			
NM_007616   	PB.5785.1 (FSM)
XM_006504974	PB.5785.2 (FSM)
                PB.5785.3 (NIC)
XM_006504975	PB.5785.4 (FSM)
                PB.5785.5 (NIC)
XM_006504976	PB.5785.6 (FSM)



chr6    Aceview exon    17306335        17306479        .       +       .       Name=12389:1;Parent=NM_007616;parent_type=mRNA
chr6    Aceview gene    17306335        17341328        .       +       .       ID=12389;Name=12389
chr6    Aceview mRNA    17306335        17341328        .       +       .       ID=NM_007616;Name=NM_007616;Parent=12389
chr6    Aceview exon    17306340        17306479        .       +       .       Name=12389:2;Parent=XM_006504974;parent_type=mRNA
chr6    Aceview mRNA    17306340        17341328        .       +       .       ID=XM_006504974;Name=XM_006504974;Parent=12389
chr6    Aceview exon    17307053        17307164        .       +       .       Name=12389:3;Parent=XM_006504975;parent_type=mRNA
chr6    Aceview mRNA    17307053        17341328        .       +       .       ID=XM_006504975;Name=XM_006504975;Parent=12389
chr6    Aceview exon    17307283        17307318        .       +       .       Name=12389:4;Parent=XM_006504976;parent_type=mRNA
chr6    Aceview mRNA    17307283        17341328        .       +       .       ID=XM_006504976;Name=XM_006504976;Parent=12389
chr6    Aceview exon    17307640        17308045        .       +       .       Name=12389:5;Parent=NM_001243064;parent_type=mRNA
chr6    Aceview mRNA    17307640        17341328        .       +       .       ID=NM_001243064;Name=NM_001243064;Parent=12389
chr6    Aceview exon    17307881        17308045        .       +       .       Name=12389:6;Parent=XM_006504974,XM_006504975,XM_006504976,NM_007616;parent_type=mRNA
chr6    Aceview exon    17339113        17339475        .       +       .       Name=12389:7;Parent=XM_006504974;parent_type=mRNA
chr6    Aceview exon    17339113        17341328        .       +       .       Name=12389:8;Parent=XM_006504975,XM_006504976,NM_001243064,NM_007616;parent_type=mRNA
chr6    Aceview exon    17340098        17341328        .       +       .       Name=12389:9;Parent=XM_006504974;parent_type=mRNA



chr6    gffutils_derived_gff3_converted exon    17306392        17306479        .       +       .       Name=PB.5785:1;Parent=PB.5785.1;parent_type=mRNA
chr6    gffutils_derived_gff3_converted gene    17306392        17341450        .       +       .       ID=PB.5785;Name=PB.5785
chr6    gffutils_derived_gff3_converted transcript      17306392        17341321        .       +       .       ID=PB.5785.1;Name=PB.5785.1;Parent=PB.5785
chr6    gffutils_derived_gff3_converted exon    17306402        17306479        .       +       .       Name=PB.5785:2;Parent=PB.5785.2;parent_type=mRNA
chr6    gffutils_derived_gff3_converted transcript      17306402        17341309        .       +       .       ID=PB.5785.2;Name=PB.5785.2;Parent=PB.5785
chr6    gffutils_derived_gff3_converted exon    17306418        17306458        .       +       .       Name=PB.5785:3;Parent=PB.5785.3;parent_type=mRNA
chr6    gffutils_derived_gff3_converted transcript      17306418        17341450        .       +       .       ID=PB.5785.3;Name=PB.5785.3;Parent=PB.5785
chr6    gffutils_derived_gff3_converted exon    17306941        17307164        .       +       .       Name=PB.5785:4;Parent=PB.5785.4;parent_type=mRNA
chr6    gffutils_derived_gff3_converted transcript      17306941        17341315        .       +       .       ID=PB.5785.4;Name=PB.5785.4;Parent=PB.5785
chr6    gffutils_derived_gff3_converted exon    17307723        17308045        .       +       .       Name=PB.5785:5;Parent=PB.5785.5,PB.5785.6;parent_type=mRNA
chr6    gffutils_derived_gff3_converted transcript      17307723        17341315        .       +       .       ID=PB.5785.6;Name=PB.5785.6;Parent=PB.5785
chr6    gffutils_derived_gff3_converted transcript      17307723        17341326        .       +       .       ID=PB.5785.5;Name=PB.5785.5;Parent=PB.5785
chr6    gffutils_derived_gff3_converted exon    17307860        17308045        .       +       .       Name=PB.5785:6;Parent=PB.5785.3;parent_type=mRNA
chr6    gffutils_derived_gff3_converted exon    17307881        17308045        .       +       .       Name=PB.5785:7;Parent=PB.5785.1,PB.5785.2,PB.5785.4;parent_type=mRNA
chr6    gffutils_derived_gff3_converted exon    17339113        17339475        .       +       .       Name=PB.5785:8;Parent=PB.5785.2,PB.5785.5;parent_type=mRNA
chr6    gffutils_derived_gff3_converted exon    17339113        17341315        .       +       .       Name=PB.5785:9;Parent=PB.5785.4,PB.5785.6;parent_type=mRNA
chr6    gffutils_derived_gff3_converted exon    17339113        17341321        .       +       .       Name=PB.5785:10;Parent=PB.5785.1;parent_type=mRNA
chr6    gffutils_derived_gff3_converted exon    17339113        17341450        .       +       .       Name=PB.5785:11;Parent=PB.5785.3;parent_type=mRNA
chr6    gffutils_derived_gff3_converted exon    17340098        17341309        .       +       .       Name=PB.5785:12;Parent=PB.5785.2;parent_type=mRNA
chr6    gffutils_derived_gff3_converted exon    17340098        17341326        .       +       .       Name=PB.5785:13;Parent=PB.5785.5;parent_type=mRNA












Here's what I have:

We intially looked at 2 genes for Ana : Rbm7 and Dhrs7b. In the SQANTI paper, only Dhrs7b is really discussed in any depth and has the most expressed transcript identified, so I started with only Dhrs7b.

Dhrs7b has 3 RefSeq transcripts, of which two have matching PacBio IDs:

RefSeqID	PacBioID	Match-type	MET PB?
NM_145428	PB.1035.1	FSM		Yes
NM_001172117	PB.1035.2	FSM		No
XM_006532869	(none)

I then looked at the proportion of events detected in each. The two transcripts with PacBio matches have all their events detected. The remaining non-PB-matching transcript (XM_006532869) is missing a unique junction between its first and second exons. Its first exon (S37379_SI, also the only other unique event for this transcript) is also inconsistently AND lowly-expressed (ie, only in one of the two reps, and its coverage is much lower compared to the other exons of Dhrs7b: mean APN od=1.32 vs >20 for all other remaining exons). So we can eliminate XM_006532869, as its unique events are essentially not detected, leaving only the two PB-matching transcripts.

When we look at the RSEM estimates of these two transcripts, we see that the most expressed transcript from SQANTI (NM_145428) is also the most expressed transcript of Dhrs7b (after filtering), although it is not that much different from the other remaining transcript.

                  perc_      perc_unique_   mean_TPM_    mean_TPM_    mean_TPM_
transcript_      events_       events_         no_        after_        pacbio
id                dtct          dtct         filtering  filtering
						
NM_001172112     1.00000          1.0	       7.79         27.98	33.70
NM_145428        1.00000          1.0	       3.04         28.99       46.76
XM_006532869     0.92308          0.5	      40.70

