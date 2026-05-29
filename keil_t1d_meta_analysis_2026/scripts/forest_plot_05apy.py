import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import matplotlib.patches as patches

# Path to the CSV files
path = '/nfshome/p.yang/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/data_4_forest_plots'
graphs_path = '/nfshome/p.yang/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/graphs/forest_plots'

# Create directory for graphs if it doesn't exist
os.makedirs(graphs_path, exist_ok=True)

# Get a list of all files in the directory
all_files = os.listdir(path)

# Filter the list for files that end with '_study.csv'
study_files = [f for f in all_files if f.endswith('_study.csv')]

# Plotting function
def plot_data(df_female, df_male, gene_name, cell_type, title, ylim, filename):
    num_fragments = len(df_female)
    fig_width = max(12, num_fragments * 2)  # Dynamically set figure width
    base_height = 1  # Base height for rectangles and ovals
    plt.figure(figsize=(fig_width, 8))  # Adjust width based on number of fragments
    ax = plt.gca()  # Get current axis
    x_positions = range(num_fragments)
    x_offset = 0.15  # Offset for males

    # Fixed width for the markers
    fixed_width = 0.2

    # Add sample plot elements for the legend
    plt.plot([], [], 's', color='red', alpha=1.0, label='Female (exon)')
    plt.plot([], [], 's', color='red', alpha=0.5, linestyle='dashed', label='Female (intron)')
    plt.plot([], [], 'o', color='blue', alpha=1.0, label='Male (exon)')
    plt.plot([], [], 'o', color='blue', alpha=0.5, linestyle='dashed', label='Male (intron)')

    # Plot for females
    for x, row in zip(x_positions, df_female.itertuples()):
        alpha = 0.5 if row.ef_ir_flag == 1 else 1.0
        linestyle = 'dashed' if row.ef_ir_flag == 1 else 'solid'
        # Adjust height based on weight proportion
        height = base_height * (row.weight / 100)  # Scale rectangle height based on weight
        rect = patches.Rectangle((x - x_offset - fixed_width / 2, row.effect_size - height / 2), fixed_width, height, 
                                 linewidth=1, edgecolor='red', facecolor='red', alpha=alpha, linestyle=linestyle)
        ax.add_patch(rect)

        # Add error bars
        plt.errorbar(x - x_offset, row.effect_size, 
                     yerr=[[row.effect_size - row.ci_lower], [row.ci_upper - row.effect_size]],
                     fmt='none', ecolor='red', alpha=alpha, capsize=5)

    # Plot for males
    for x, row in zip(x_positions, df_male.itertuples()):
        alpha = 0.5 if row.ef_ir_flag == 1 else 1.0
        linestyle = 'dashed' if row.ef_ir_flag == 1 else 'solid'
        # Adjust height based on weight proportion
        height = base_height * (row.weight / 100)  # Scale oval height based on weight
        oval = patches.Ellipse((x + x_offset, row.effect_size), fixed_width, height,
                               linewidth=1, edgecolor='blue', facecolor='blue', alpha=alpha, linestyle=linestyle)
        ax.add_patch(oval)

        # Add error bars
        plt.errorbar(x + x_offset, row.effect_size, 
                     yerr=[[row.effect_size - row.ci_lower], [row.ci_upper - row.effect_size]],
                     fmt='none', ecolor='blue', alpha=alpha, capsize=5)

    # Plot horizontal lines for moderator effect
    moderator_effect_female = df_female['moderator_effect'].iloc[0]
    moderator_effect_male = df_male['moderator_effect'].iloc[0]
    plt.axhline(moderator_effect_female, color='red', linestyle='-', linewidth=2, alpha=0.2, label='Female Moderator Effect')
    plt.axhline(moderator_effect_male, color='blue', linestyle='-', linewidth=2, alpha=0.2, label='Male Moderator Effect')

    plt.axhline(0, color='grey', linestyle='--', alpha=0.5)  # Add grey transparent dotted line at y=0
    plt.xlabel('Exon Fragments')
    plt.ylabel('Effect Size', labelpad=20)  # Ensure y-axis title is fully displayed
    plt.title(f'Effect Size for {gene_name} ({cell_type}) by Fragment')
    plt.grid(True, alpha=0.5)
    plt.xlim(-0.5, num_fragments - 0.5)  # Ensure even spacing
    plt.ylim(ylim)
    plt.xticks(ticks=x_positions, labels=df_female['featureID'], rotation=90)  # Label x-axis with featureID
    plt.legend()
    plt.tight_layout(pad=3.0, w_pad=3.0, h_pad=3.0)  # Adjust layout to ensure the y-axis title is fully displayed
    # plt.savefig(os.path.join(graphs_path, filename))  # Comment out this line to avoid saving the plot
    plt.show()
    #plt.close()

# Loop through each study file
for file in study_files:
    # Extract the gene name and cell type from the filename
    gene_name = re.search(r'ENSG\d+', file).group(0)
    cell_type_match = re.search(r'CD[48]', file)
    cell_type = cell_type_match.group(0) if cell_type_match else 'Unknown'
    
    # Read the CSV file
    df = pd.read_csv(os.path.join(path, file))
    
    # Create a new column 'featureID' excluding the gene name for the x-axis
    df['featureID'] = df['featureID'].apply(lambda x: ':'.join(x.split(':')[1:]))
    
    # Filter by sex
    df_female = df[df['sex'] == 'F'].sort_values(by=['ef_start'], ascending=True)
    df_male = df[df['sex'] == 'M'].sort_values(by=['ef_start'], ascending=True)
    
    # Determine y-axis range based on the actual bounds
    min_y = min(df_female['ci_lower'].min(), df_male['ci_lower'].min())
    max_y = max(df_female['ci_upper'].max(), df_male['ci_upper'].max())
    y_range = (min_y - 0.25, max_y + 0.25)
    
    # Combined plot for females and males
    plot_data(df_female, df_male, gene_name, cell_type, 'Females and Males by Fragment', ylim=y_range, filename=f'ef_sz_{gene_name}_{cell_type}_female_male_ef.png')
