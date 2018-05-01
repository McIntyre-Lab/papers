
/* compare MEI and MEE between genes */

data mei;
  set event.t1d_mei_comparison_by_gene;
run;

data mee;
  set eventloc.t1d_mee_comparison_by_gene;
run;

proc sort data=mei;
  by gene_id;
proc sort data=mee;
  by gene_id;
run;

data mee_vs_mei;
   merge mee (in=in1) mei (in=in2);
   by gene_id;
   if in1 and in2;
run;

proc freq data=mee_vs_mei noprint;
  tables flag_cd4cd8_mee_diff*flag_cd4cd19_mee_diff*flag_cd8cd19_mee_diff*
         flag_mono_isoform_gene / out=exon_check;
proc print data=exon_check;
run;
quit;

/*
                                                  flag_mono_
flag_cd4cd8_    flag_cd4cd19_    flag_cd8cd19_     isoform_
  mee_diff         mee_diff         mee_diff         gene       COUNT

      0               0                0               1         2766
      0               1                1               1           31
      1               0                1               1           16
      1               1                0               1           13
      1               1                1               1            1
      0               0                0               0         2725
      0               1                1               0          151
      1               0                1               0           52
      1               1                0               0           47
      1               1                1               0            8

*/

proc freq data=mee_vs_mei noprint;
  tables flag_cd4cd8_mei_diff*flag_cd4cd19_mei_diff*flag_cd8cd19_mei_diff*
         flag_mono_isoform_gene / out=iso_check;
proc print data=iso_check;
run;
quit;


/*
                                                  flag_mono_
flag_cd4cd8_    flag_cd4cd19_    flag_cd8cd19_     isoform_
  mei_diff         mei_diff         mei_diff         gene       COUNT    PERCENT

      0               0                0               0         2722    46.8503
      0               0                0               1         2827    48.6575
      0               1                1               0          169     2.9088
      1               0                1               0           41     0.7057
      1               1                0               0           44     0.7573
      1               1                1               0            7     0.1205

*/

proc freq data=mee_vs_mei noprint;
  where flag_mono_isoform_gene=0;
  tables flag_cd4cd8_mee_diff*flag_cd4cd19_mee_diff*flag_cd8cd19_mee_diff*
         flag_cd4cd8_mei_diff*flag_cd4cd19_mei_diff*flag_cd8cd19_mei_diff/ out=gene_check;
proc print data=gene_check;
run;
quit;

/*
Only genes with multiple isoforms in the reduced transcriptome set (APN>5, 75% dtct):


flag_cd4cd8_ flag_cd4cd19_ flag_cd8cd19_ flag_cd4cd8_ flag_cd4cd19_ flag_cd8cd19_
  mee_diff      mee_diff      mee_diff     mei_diff      mei_diff      mei_diff   COUNT

No differences:
      0            0             0             0            0             0        2502

Only isoform differences:
      0            0             0             0            1             1         145
      0            0             0             1            0             1          36
      0            0             0             1            1             0          36
      0            0             0             1            1             1           6

Only exon differences:
      0            1             1             0            0             0         126
      1            0             1             0            0             0          49
      1            1             0             0            0             0          38
      1            1             1             0            0             0           7

Isoform and exon differences:
      0            1             1             0            1             1          16
      0            1             1             1            0             1           3
      0            1             1             1            1             0           6
      1            0             1             0            1             1           3
      1            1             0             0            1             1           5
      1            1             0             1            0             1           2
      1            1             0             1            1             0           2
      1            1             1             1            1             1           1

Of genes with multiple isoforms in the reduced set:
2502 genes with no differences in MEE or MEI
223 genes with only MEI differences
220 genes with only MEE differences
38 genes with both MEI and MEE differences
*/


