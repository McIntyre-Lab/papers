library(ComplexHeatmap)
library(circlize)
library(seriation)

#groups <- list("TCA", "UGT", "NI", "NInoN2")
groups <- list("ALL")
types <- list("rp_pos", "rp_neg", "hilic_pos")
#types <- list("rp_pos")

for (group in groups) {

    for (type in types) {
#group = 'TCA'

# Get input data
        data = read.csv(paste("~/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/MA_effectSizes_", type, "_", group, ".csv", sep=""))

# Set rowID as index
        rownames(data) <- data$rowID

# Reset second columns name
        colnames(data)[2] <- 'Meta_pathway'

# Get column index of the last mutant (first mutant starts on index 3)
        lastcol = 16
# Get mutant names
        for (i in 3:lastcol) {
                colnames(data)[i] <- strsplit(colnames(data)[i], "_")[[1]][2]
        }

# Subset main dataset into
# 1. other info: 'flag_same_direction', 'sig_counts', 'direction', 'p_meta_pathway'
# 2. effect sizes only
#        info = data[,c((lastcol+1):(lastcol+4))]
        effectS = as.matrix(data[,c(2:lastcol)])

# Set columns for heatmap: one "pathway model" and "within-set-model" repeated for each mutant
        colsplit = c("pathway\nmodel", rep('within-set-model', (lastcol-2)))

# Set direction AMM changes neg_val to le_PD1074_val and pos_val to ge_PD1074
        #info$direction[info$direction == 0] <- 'le_PD1074_val'
        #info$direction[info$direction == 1] <- 'ge_PD1074_val'
#        info$direction[info$direction == 0] <- 'neg_val'
#        info$direction[info$direction == 1] <- 'pos_val'

# Set colors vector
        colors = c('#ff4d4d','#1a75ff','#ff8c1a', '#cc66ff','#00e6ac','black')

# Make sig_in_N_model and positive/negative heatmaps (far right)
#        row_ha = rowAnnotation(sig_in_N_model = as.factor(info$sig_counts), 
#                               direct = as.factor(info$direction),
#                               col = list(direct = c("neg_val" = "#99b3e6", 
#                                                     "pos_val" = "#ff3399"),
#                                          sig_in_N_model = c('0' = colors[1],
#                                                  '1' = colors[2],
#                                                  '2' = colors[3],
#                                                  '3' = colors[4],
#                                                  '4' = colors[5],
#                                                  '5' = colors[6])), show_annotation_name = FALSE)

# Set p-value colors
        col_fun = colorRamp2(c(0, 0.05, 1), c("#263300", '#dfff80',"#ffffff"))

# Make p-value heatmap (far left)
#        row_ha_left = rowAnnotation(p_pathway = info$p_meta_pathway, show_annotation_name = FALSE, col = list(p_pathway = col_fun))

        print(paste("Number of rows is", nrow(effectS)))

## clustered heatmap
# Set up output PDF file
#        pdf(paste('/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/MA_heatmap_', type, "_", group, '.pdf', sep = ""))
        svg(paste('/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/MA_heatmap_', type, "_", group, '_complete_clustered.svg', sep = ""))

# Plot full heatmap
        print(Heatmap(effectS, name = "Effect",
#                right_annotation = row_ha,
#                left_annotation = row_ha_left,
                column_names_rot = 45,
#                row_split = info$sig_counts,
                column_split = colsplit,
                cluster_rows = TRUE, 
                cluster_columns = TRUE, 
                show_row_names = FALSE,
                row_title = c("metabolites"),
                column_title_gp = gpar(fontsize = c(9,rep(12,5))),
                column_names_gp = gpar(col = c("purple", 'orange'), 
                                        fontsize = c(10, 12))))
    
        dev.off()

## not-clustered heatmap 
# Set up output PDF file
#        pdf(paste('/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/MA_heatmap_', type, "_", group, '.pdf', sep = ""))
        svg(paste('/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/MA_heatmap_', type, "_", group, '_complete.svg', sep = ""))

# Plot full heatmap
        print(Heatmap(effectS, name = "Effect",
#                right_annotation = row_ha,
#                left_annotation = row_ha_left,
                column_names_rot = 45,
#                row_split = info$sig_counts,
                column_split = colsplit,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                row_title = c("metabolites"),
                column_title_gp = gpar(fontsize = c(9,rep(12,5))),
                column_names_gp = gpar(col = c("purple", 'orange'),
                                        fontsize = c(10, 12))))

        dev.off()


    }
}
