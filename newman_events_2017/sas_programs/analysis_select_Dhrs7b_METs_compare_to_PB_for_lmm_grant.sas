ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Most expressed RefSeq transcript = Most expressed PacBio transcript

For Dhrs7b, we need to show:
(1) We select the same set of transcripts as SQANTI
	- OR, we eliminate the same transcripts that aren't detected by PB
	- this is easier
(2) For the remaining set of transcripts, we need to show that reducing transcriptomic complexity
    results in the same transcript classified as the "most expressed transcript" as in SQANTI

So, we have 3 Dhrs7b transcripts:
RefSeqID	PacBioID	Match-type
NM_145428	PB.1035.1	FSM
NM_001172117	PB.1035.2	FSM
XM_006532869	(none)

Which do we eliminate and which do we keep?
How does their quantification compare?
(Also: are there any additional novel junctions??) */

/* Subset Dhrs7b RefSeq transcripts, and how much of their events are detected */

data xs2gene;
   set event.feature2xs2gene_exp_only_nomulti;
   where gene_id="216820";
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

data dhrs7b_xs_info;
  merge xs2gene (in=in1) xs_info (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc print data=dhrs7b_xs_info(keep=transcript_id perc_features_dtct perc_unique_features_dtct);
run;

/*

                  perc_      perc_unique_
transcript_     features_      features_
id                 dtct          dtct
					
NM_001172112     1.00000          1.0
NM_145428        1.00000          1.0
XM_006532869     0.92308          0.5		

*/


/* Missing pieces of XM_006532869 */

data frags_off;
  set event.flag_fragment_on;
  where flag_fragment_nsc_on=0;
  keep fragment_id;
  rename fragment_id=feature_id;
run;

data fus_off;
  set event.flag_fusion_on;
  where flag_fusion_nsc_on=0;
  keep fusion_id;
  rename fusion_id=feature_id;
run;

data juncs_off;
   set event.flag_splicing_on;
   where flag_event_nsc_on=0;
   keep event_id;
   rename event_id=feature_id;
run;

data events_off;
   set juncs_off frags_off fus_off;
run;

data event2xs;
  set event.feature2xs2gene_exp_only_nomulti;
   where gene_id="216820";
run;

proc sort data=events_off;
   by feature_id;
proc sort data=event2xs;
   by feature_id;
run;

data events_off_by_xs;
  merge event2xs (in=in1) events_off (in=in2);
  by feature_id;
  if in1 and in2;
run;

/* 1 event off and it is the junction 216820:3_4, belonging to XM_006532869 */

/* check coverage of fusions */

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

  1     F37378_SI        32.209923664    NM_001172112|NM_145428
  2     F37385_SI        62.879066478    NM_001172112|NM_145428|XM_006532869
  3     S37379_SI        1.3235294118    XM_006532869
  4     S37380_SI                84.7    NM_001172112|NM_145428|XM_006532869
  5     S37381_SI        46.147619048    NM_001172112|NM_145428|XM_006532869
  6     S37382_SI        130.05733945    NM_001172112|NM_145428|XM_006532869
  7     S37383_SI        107.41489362    NM_001172112|NM_145428|XM_006532869
  8     S37384_SI        131.25649351    NM_001172112|NM_145428|XM_006532869

*/

/* Now for these transcripts, get the mean TPM from RSEM, from:
   all RefSeq transcripts
   all EA transcripts
   EA transcripts with 100% events detected (APN>0)
   all pacbio transcripts
*/

data tpm_all;
   set event.rsem_refseq_all;
   where transcript_id in ("XM_006532869","NM_145428","NM_001172112");
   mean_TPM_refseq_all=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_TPM_refseq_all;
run;

data tpm_any;
   set event.rsem_events_exp_any;
   where transcript_id in ("XM_006532869","NM_145428","NM_001172112");
   mean_TPM_event_all=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_TPM_event_all;
run;

data tpm_100;
   set event.rsem_events_exp_100perc;
   where transcript_id in ("XM_006532869","NM_145428","NM_001172112");
   mean_TPM_event_100=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_TPM_event_100;
run;

data tpm_pb;
   length transcript_id $20.;
   set event.rsem_pacbio_all;
   where transcript_id in ("PB.1035.1","PB.1035.2");
   if transcript_id="PB.1035.1" then transcript_id="NM_145428";
   if transcript_id="PB.1035.2" then transcript_id="NM_001172112";
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

  1     NM_001172112       7.790         9.570       27.980        33.70
  2     NM_145428          3.035         3.745       28.985        46.76
  3     XM_006532869      40.695        49.240         .             .
*/




/* 

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

