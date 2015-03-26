libname fru '/mclab/Fru_network/sasdata';
libname fru '/home/Justin/fru';

proc import out=work.alignment_counts
            datafile="/mclab/Fru_network/pipeline_output/fru_network_alignments.csv"
            dbms=csv replace;
            getnames=yes;
            datarow=2;
            run;

data fru.alignment_counts;
    set alignment_counts;
    total_aln = spikes_aln + dmel_aln;
    percent_spike = (spike_aln / total_aln)*100;
    run;
