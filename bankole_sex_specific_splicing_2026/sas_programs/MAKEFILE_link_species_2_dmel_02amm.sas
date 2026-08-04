/* sas MAKEFILE ---> link dsim2/dsan/dyak/dser fivespecies to dmel  */


/* link each species to dmel */
%include "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sas_programs/link_dsim2_fivespecies_2_dmel_04amm.sas" ;
%include "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sas_programs/link_dsan_fivespecies_2_dmel_04amm.sas" ;
%include "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sas_programs/link_dser_fivespecies_2_dmel_05amm.sas" ;
%include "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sas_programs/link_dyak_fivespecies_2_dmel_04amm.sas" ;

/* add geneids */
%include "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sas_programs/add_geneids_fivespecies_03ksb.sas" ;

/* check for sex det genes */
%include "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sas_programs/sex_det_genes_by_mapping_03amm.sas" ;

/* link dmel 5 species to other species */
%include "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sas_programs/link_dmel_fivespecies_2_others.sas" ;
