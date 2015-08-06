# Set up major variables
PROJ=$MCLAB/cegs_sem_sd_paper
WORK=$PROJ/analysis_output/adding_genes/cegsV_add_sex_members
MNAME=cegsv_by_gene_sbs
ORIG=$PROJ/sasdata/${MNAME}.sas7bdat
TMPDIR=/home/jfear/tmp

# Make a copy of the orignal data
cd $TMPDIR
cp $ORIG $TMPDIR/

for GENE in dsx fl_2_d fru her ix snf Spf45 Sxl tra2 tra vir yp2
do
    if [ $GENE == 'yp2' ]
    then
        PYGEN=$PROJ/scripts/cegsV_adding_yp2.py
    else
        PYGEN=$PROJ/scripts/cegsV_adding_${GENE}_yp2.py
    fi

    # Generate SAS models for GENE
    if [[ ! -e $WORK/sas_programs/$GENE ]]; then mkdir -p $WORK/sas_programs/$GENE; fi;
    if [[ ! -e $WORK/model_logs/$GENE ]]; then mkdir -p $WORK/model_logs/$GENE; fi;
    if [[ ! -e $WORK/gen_logs ]]; then mkdir -p $WORK/gen_logs; fi;
    if [[ ! -e $WORK/comb_logs ]]; then mkdir -p $WORK/comb_logs; fi;
    if [[ ! -e $WORK/sasdata ]]; then mkdir -p $WORK/sasdata; fi;

    echo "[`date`] Generating SEM models for $GENE"
    python $PYGEN -o $WORK/sas_programs/$GENE/$GENE.sas -l $TMPDIR -m $MNAME -g $GENE -n $GENE --log $WORK/gen_logs/${GENE}_model_gen.log
    echo "[`date`] Finished generating SEM models for $GENE"

    # Run SAS Models
    echo "[`date`] Running SEM models for $GENE"
    for MODEL in $WORK/sas_programs/$GENE/$GENE*.sas
    do
        NAME=`basename $MODEL .sas`
        sas -work $TMPDIR -MEMSIZE 3g -nonews -sysin $MODEL -log $WORK/model_logs/$GENE/${NAME}.log -print /dev/null
    done
    echo "[`date`] Finished running SEM models for $GENE"

    echo "[`date`] Combining SEM models for $GENE"
    sas -work $TMPDIR -MEMSIZE 3g -nonews -sysin $PROJ/sas_programs/dsrp_combine_genome_wide_sem_models.sas -log $WORK/comb_logs/${GENE}_combine.log -print /dev/null -sysparm "lib1=$TMPDIR,lib2=$WORK/sasdata,gene=$GENE"
    echo "[`date`] Finished combining SEM models for $GENE"

    # Delete current gene models to prepare for next iteration
    rm gene_*.sas7bdat
done
