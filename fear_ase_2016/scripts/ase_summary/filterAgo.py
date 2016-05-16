import pandas as pd
import matplotlib.pyplot as plt
from munging import mergeAndDrop


# Function that takes various cutoffs and makes a panel plot
def filterAgo(df, apn_cutoff=50, per_geno_cutoff=25, plot=True):
    # Only look at regions with APN >= apn_cutoff
    dfI = df.reset_index()
    dfI.set_index(['line', 'mating_status', 'fusion_id'], inplace=True)
    flag_apn = pd.DataFrame((dfI['mean_apn'] < apn_cutoff).astype(int))
    flag_apn.columns = ['flag_apn_lt_{}'.format(apn_cutoff)]
    flag_apn.reset_index(inplace=True)

    dfAPN = df[df['mean_apn'] >= apn_cutoff]

    # Num genotypes
    nGeno = float(dfAPN.groupby('line').ngroups)

    # Only look at exonic regions with at least per_geno_cutoff genotypes
    grp = dfAPN.groupby(['fusion_id', 'mating_status'])
    numGenoPerFusbyMS = grp['line'].count().to_frame()
    flag_geno = ((numGenoPerFusbyMS['line'] / nGeno * 100 >= per_geno_cutoff).astype(int)).reset_index()
    flag_geno.rename(columns={'line': 'flag_geno'}, inplace=True)

    dfGeno = mergeAndDrop(dfAPN, flag_geno, left_on=['fusion_id', 'mating_status'], right_on=['fusion_id', 'mating_status'], flagName='flag_geno', keep_logic='> 0')

    # Only look at exonic regions that are present in mated and virgin
    grp = dfGeno.groupby(['fusion_id', 'line'])
    numMsPerFusbyLine = grp['mating_status'].count().to_frame()
    flag_bothMS = ((numMsPerFusbyLine['mating_status'] == 2).astype(int)).reset_index()
    flag_bothMS.rename(columns={'mating_status': 'flag_bothMS'}, inplace=True)

    dfBothMS = mergeAndDrop(dfGeno, flag_bothMS, left_on=['fusion_id', 'line'], right_on=['fusion_id', 'line'], flagName='flag_bothMS', keep_logic='> 0')

    # Make mated and virgin into side-by-side dataset for plotting
    mated = dfBothMS[dfBothMS['mating_status'] == 'M']
    virgin = dfBothMS[dfBothMS['mating_status'] == 'V']
    clean = mated.merge(virgin, on=['line', 'fusion_id'], suffixes=['_m', '_v'])

    # Get information at the genotype level
    grp = clean.groupby('line')

    ## Number of exonic regions remaining in the dataset
    fusCnt = grp['fusion_id'].count()
    print "Min num exons: {}\nMax num exons: {}".format(fusCnt.min(), fusCnt.max())

    ## Number of exonic regions with AI
    AIM = grp['flag_AI_combined_m'].sum()
    AIV = grp['flag_AI_combined_v'].sum()

    ## Make Figure
    if plot:
        fig, axs = plt.subplots(8, 9, figsize=(30, 30))
        axs = axs.ravel()

        ## Iterate over genotypes and plot
        cnt = 0
        for name, group in grp:
            # Percent of exonic regions with AI
            aim = AIM[name] / float(fusCnt[name]) * 100
            aiv = AIV[name] / float(fusCnt[name]) * 100

            # Color plot by AI
            ## 0: blue -> No AI
            ## 1: cyan -> Mated with AI
            ## 2: yellow -> Virgin with AI
            ## 3: red -> Mated and Virgin with AI
            group.loc[:, 'colors'] = 0
            group.loc[group['flag_AI_combined_m'] == 1, 'colors'] = 1
            group.loc[group['flag_AI_combined_v'] == 1, 'colors'] = 2
            group.loc[(group['flag_AI_combined_m'] == 1) & (group['flag_AI_combined_v'] == 1), 'colors'] = 3

            # Plot
            group.plot('q5_mean_theta_m', 'q5_mean_theta_v', kind='scatter', ax=axs[cnt], c='colors', colormap='rainbow', colorbar=True)

            # Tweak plot
            axs[cnt].set_title('{}\nNumber of Exonic Regions: {}\nPercent AI M: {}\nPercent AI V: {}'.format(name, fusCnt[name], aim, aiv), fontsize=10)
            axs[cnt].set_xlabel('mated')
            axs[cnt].set_ylabel('virgin')
            axs[cnt].axhline(0.5, color='r', lw=2)
            axs[cnt].axvline(0.5, color='r', lw=2)
            axs[cnt].set_xlim(0, 1)
            axs[cnt].set_xticks([])
            axs[cnt].set_ylim(0, 1)
            axs[cnt].set_yticks([])
            cnt += 1
        # Remove unsued axes from panel
        [ax.axis('off') for ax in axs[-4:]]
        plt.tight_layout()

    return fusCnt, clean, (flag_apn, flag_geno, flag_bothMS)

if __name__ == '__main__':
    pass
