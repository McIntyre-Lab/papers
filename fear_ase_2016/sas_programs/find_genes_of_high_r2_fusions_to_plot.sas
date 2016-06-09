libname dmel "Z:\useful_dmel_data\flybase551\sasdata";

data gene_id;
set dmel.fb551_si_fusion_2_gene_id;
keep fusion_id chrom start end exon_gene_id symbol;
run;

%macro find (fusion);
data &fusion;
set gene_id;
if fusion_id = "&fusion";
run;
%mend;
%find(F11767_SI);
%find(F13922_SI);
%find(F15445_SI);
%find(F10005_SI);
%find(F15876_SI);

%find(F15765_SI);
%find(F37624_SI);
%find(S682_SI);
%find(S52856_SI);

%find(F587_SI);
%find(F1558_SI);
%find(F21007_SI);
%find(S18508_SI);
%find(S38687_SI);

%find(F32473_SI);
%find(S10318_SI);
%find(S13419_SI);
%find(S2713_SI);

*do these genes show up on the splice list? ;
libname cegs "Z:/cegs_ase_paper/sas_data";

data splice;
set cegs.genes_qith_spice_ai;
run;

%macro gene (symbol);
data &symbol;
set splice;
if symbol_cat = "&symbol";
run;

%mend;
%gene(Acn); *no;
%gene(Appl); *yes;
%gene(CG15765); *no;
%gene(His3.3B); *no;
%gene(l(1)G0230); *no;

%gene(Idgf4); *no;
%gene(brp); *no;
%gene(IA-2); *yes;
%gene(Men); *no;

%gene(CG3662); *no;
%gene(Syt1); *yes;
%gene(Opb19d); *no;
%gene(CG8974); *no;
%gene(Listericin); *no;

%gene(Pka-R1); *no;
%gene(Tep4); *no;
%gene(CHOp24); *no;
%gene(CG3036); *no;
