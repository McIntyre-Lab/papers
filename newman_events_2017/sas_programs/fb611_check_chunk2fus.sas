/* Libraries */

ods listing; ods html close;
libname dmel611 '!MCLAB/useful_dmel_data/flybase611/sas_data';


/* Import exon chunks and check that starts match with fusions */

proc import datafile='!MCLAB/conesa_pacbio/created_files/fb611_exon_chunks.csv' out=chunks
     dbms=csv replace; guessingrows=88000;
run;


/* get fusions */

data fus_info;
   set dmel611.fb611_si_fusion_2_gene_id;
   keep chrom start end fusion_id exon_name FBgn;
run;

data chunks2;
   set chunks;
   start=group_start;
   rename group_stop=end;
run;

proc sort data=fus_info nodup;
   by fusion_id exon_name;
run;

proc sort data=chunks2 nodup;
   by exon_id;
run;

proc sort data=fus_info;
   by chrom start;
proc sort data=chunks2;
   by chrom start;
run;

data fus2chunk nochunk nofus;
   merge fus_info (in=in1) chunks2 (in=in2);
   by chrom start;
   if in1 and in2 then output fus2chunk;
   else if in1 then output nochunk;
   else output nofus;
run;

data chunks3;
length exon_name $35.; 
set chunks2; 
do i=1 by 1 while(scan(exon_id,i,'|') ^=' '); 
exon_name =scan(exon_id,i,'|'); 
output; 
end; 
run;


proc sort data=fus_info;
   by exon_name;
proc sort data=chunks3;
   by exon_name;
run;


data fus2chunk nochunk nofus;
   merge fus_info (in=in1) chunks3 (in=in2);
   by exon_name;
   if in1 and in2 then output fus2chunk;
   else if in1 then output nochunk;
   else output nofus;
run;



