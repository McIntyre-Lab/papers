

/*edges are the link between the components and the genes.  This is in the component_map_by_node.csv that is output from the initial network*/

PROC IMPORT OUT= component_map
            DATAFILE= "!MCLAB/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/component_map_by_node.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
run;

     data WORK.COMPONENT_MAP    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile  '!MCLAB/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/component_map_by_node.csv' 
     delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat jxnhash $64. ;
        informat geneID $12. ;
        informat source $5. ;
        informat component_id $ 12. ;
        format jxnhash $64. ;
        format geneID $12. ;
        format source $5. ;
        format component_id $ 12. ;
     input
                 jxnhash  $
                 geneID  $
                 source  $
                 component_id $
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;
     
data component_map1;
set component_map;
drop jxnhash;
run;

proc contents data=component_map1;
run;

/*keep the unique links between the components and the genes*/

proc sort data= component_map1 nodupkey out=geneset_edges;
by component_id geneid;
run;


proc sort data=geneset_edges nodupkey out=component_list;
by component_ID;
run;

/*56370*/


proc sort data=geneset_edges nodupkey out=gene_list;
by geneID;
run;


/*88259*/

data gene_list1;
set gene_list;
rename geneID=nodeID;
drop component_id source;


data component_list1;
set component_list;
rename component_id=nodeID;
drop geneid source;
run;


data geneset_nodes;
set gene_list1 component_list1;
run;

/*144629*/


proc export data = geneset_nodes
outfile = "!MCLAB/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/geneset_nodes.csv"
dbms = csv replace ;
run;

proc export data = geneset_edges
outfile = "!MCLAB/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/geneset_edges.csv"
dbms = csv replace ;
run;

