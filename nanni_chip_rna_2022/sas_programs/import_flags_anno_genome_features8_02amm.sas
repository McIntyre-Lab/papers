/*libname RNA 'Z:\SHARE\McIntyre_Lab\Dros_PB_ChIP\RNAseq\group_summary';
libname chip 'Z:\SHARE\McIntyre_Lab\Dros_PB_ChIP\ChIPseq\group_summary';
libname anno 'Z:\SHARE\McIntyre_Lab\Dros_PB_ChIP\annotation\sequenced_transcripts';
libname fb 'Z:\SHARE\McIntyre_Lab\useful_dmel_data\flybase617\sas_data';

libname local 'C:\a1stuff\dros';
*/

libname chip "!MCLAB/Dros_PB_ChIP/ChIPseq/group_summary";
libname anno "!MCLAB/Dros_PB_ChIP/annotation/sequenced_transcripts";
libname fb "!MCLAB/useful_dmel_data/flybase617/sas_data"; 


 %macro import(id,in,out,path);
proc import out=work.&out
 datafile="&path.&in..csv"
 dbms=csv replace;
 getnames=yes;
 guessingrows=MAX;
 datarow=2;
 run;

 proc sort data=&out;
 by &id;
 run;

 %mend;

/*chip data detection*/

%import(featureid,mel_chip_5U_3U_TSS_frag_intr_inter_flag,mel_chip_frag_flags, !MCLAB/Dros_PB_ChIP/ChIPseq/detection_above_background/features/);
%import(featureid,sim_chip_5U_3U_TSS_frag_intr_inter_flag,sim_chip_frag_flags, !MCLAB/Dros_PB_ChIP/ChIPseq/detection_above_background/features/);

 /*RNA data bias*/

 %import(featureid,mel_5U_3U_TSS_frag_intr_inter_RNA_sex_bias,mel_RNA_bias,!MCLAB/Dros_PB_ChIP/RNAseq/group_summary/);
 %import(featureid,sim_5U_3U_TSS_frag_intr_inter_RNA_sex_bias,sim_RNA_bias,!MCLAB/Dros_PB_ChIP/RNAseq/group_summary/);

/* following does not exist..... 
proc contents data=sim_frag_anno;
run;
*/

 /*RNA data detection */

%import(featureid, mel_5U_3U_TSS_frag_intr_inter_on_off_uq_ff,mel_RNA_frag_detect_flags,!MCLAB/Dros_PB_ChIP/RNAseq/sample_summary/);
%import(featureid, sim_5U_3U_TSS_frag_intr_inter_on_off_uq_ff,sim_RNA_frag_detect_flags,!MCLAB/Dros_PB_ChIP/RNAseq/sample_summary/);

 /*fragment annotations*/

%import(featureid,mel_5U_3U_TSS_frag_intr_inter_uniq,mel_frag_anno,!MCLAB/useful_dmel_data/flybase617/fb_features/);
%import(featureid,sim_5U_3U_TSS_frag_intr_inter_uniq,sim_frag_anno,!MCLAB/useful_dsim_data/flybase202/fb_features/);


/*need fbgn2go*/

%import(fbgn,gene2go_with_terms,fbgn_2go,!MCLAB/useful_dmel_data/flybase617/dmel_annotation/);

/*need symbol for both species... note that while fbgn_2go has mel symbol no sim data avaialble and no secondary fbgn*/
%import(primary_fbgn,fbgn_annotation_ID,fbgn_anno,!MCLAB/useful_dmel_data/flybase617/dmel_annotation/);

%import(primary_fbgn,fbgn2coord,fbgn2coord,!MCLAB/useful_dmel_data/flybase617/dmel_annotation/);

%import(primary_fbgn,fbgn2coord,sim_fbgn2coord,!MCLAB/useful_dsim_data/flybase202/dsim_annotation/);

%import(FlyBase_FBgn,fbgn_fbtr_fbpp,fbgn_fbtr_fbpp, !MCLAB/useful_dmel_data/flybase617/dmel_annotation/);


data fbgn_anno_all;
set fbgn_anno;
rename primary_fbgn=fbgn;
run;

proc sort data=fbgn_anno_all;
by fbgn;
run;

data fbgn_2_coord;
set fbgn2coord;
rename primary_fbgn=fbgn;
run;

proc sort data=fbgn_2_coord;
by fbgn;
run;


data sim_fbgn_2_coord;
set sim_fbgn2coord;
rename primary_fbgn=fbgn;
run;

proc sort data=sim_fbgn_2_coord;
by fbgn;
run;

data fbgn_fbtr;
set fbgn_fbtr_fbpp;
rename flybase_fbgn=fbgn;
drop flybase_fbpp;
run;

proc sort data=fbgn_fbtr nodupkey;
by fbgn flybase_fbtr;
run;


proc freq data=fbgn_fbtr;
tables fbgn/out=count_transcripts;
run;

data count_transcripts1;
set count_transcripts;
rename count=num_transcripts;
if count>1 then flag_single_transcript=0;
else if count=1 then flag_single_transcript=1;
run;

/*191883 fbgn values*/

/* need ortholist */
*MCLAB/useful_dmel_dsim_compare/dsim_GO_from_ortho/gene2go_with_terms_from_mel_ortho.csv;

%import(mel_geneId,dmel_orthologs_dsim_fb_2017_04,mel_2_sim,!MCLAB/useful_dmel_dsim_compare/);




