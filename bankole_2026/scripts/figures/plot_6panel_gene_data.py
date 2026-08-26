#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 15:38:23 2026

@author: mgaran
"""

import os
import gc
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

from plot_6panel_utils import (
    ALL_SPECIES_DEFAULT,
    SPECIES_NAMES,
    load_jxn_datafile,
    load_gtf_coords_by_jxnhash,
    auto_flip_by_strand,
    flip_coords_to_display,
    load_components_from_list,
    load_read_props_by_jxnhash,
    sort_components_by_erp,
    create_component_alignment,
    load_unannotated_jxns_in_region,
    get_gene_region_from_raw_coords,
    get_jxnhash_to_erp_map,
    get_jxnhash_to_erp_map_for_jxns,
    plot_annotated_panel,
    plot_novel_panel,
    plot_component_panel,
    plot_novel_label_panel,
    get_colormap,
    compute_shared_xlim,
    parse_species_paths,
    build_index_csv,
)

# Species with analyzable UJC effect sizes.
# dyak2 and dsan1 lack effect size data — sig filters are skipped for them.
SIG_SPECIES = {'dmel6', 'dsim2', 'dser1'}


# ---------------------------------------------------------------------------
# Novel UJC filtering
# ---------------------------------------------------------------------------

def filter_novel(candidates, cand_es, cand_total, cand_erp_anno,
                  erp_anno_flag, min_reads, sig_only, top_n, sp, label):
    """Filter and cap novel candidate jxnHashes.

    Parameters
    ----------
    erp_anno_flag : True  → keep flag_ERPanno==1 only
                    False → keep flag_ERPanno==0 only
                    None  → keep all
    """
    if erp_anno_flag is True:
        pool = [j for j in candidates if cand_erp_anno.get(j, 0) == 1]
    elif erp_anno_flag is False:
        pool = [j for j in candidates if cand_erp_anno.get(j, 0) != 1]
    else:
        pool = list(candidates)

    meets_reads = [j for j in pool if cand_total.get(j, 0) >= min_reads]
    n_below = len(pool) - len(meets_reads)
    kept = meets_reads
    print(f'[{label}] {sp}: {len(kept)} UJC(s) kept '
          f'({n_below} below min_reads={min_reads}).')

    if sig_only and sp in SIG_SPECIES:
        n_before = len(kept)
        kept = [j for j in kept
                if abs(cand_es.get(j, float('nan'))) > 2.78]
        print(f'[{label}] {sp}: {len(kept)} significant biased UJC(s) '
              f'({n_before - len(kept)} dropped by sig_only).')

    if top_n is not None and len(kept) > top_n:
        kept = sorted(kept,
                      key=lambda j: cand_total.get(j, 0),
                      reverse=True)[:top_n]
        print(f'[{label}] {sp}: trimmed to {len(kept)} UJCs '
              f'(top-{top_n} cap by expression).')

    return kept


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def plot_gene_data(
        component_file, component_list_file,
        anno_gtf_files, jxn_datafiles, data_gtf_files,
        gene_name, output_dir,
        panels='UJCanno',
        anno_min_reads=1,
        erp_anno_min_reads=1,
        erp_novel_min_reads=1,
        erp_anno_top_n=None,
        erp_novel_top_n=None,
        erp_anno_sig_only=False,
        erp_novel_sig_only=False,
        erp_datafiles=None,
        output_format='png',
        clamp_xlim=False):
    """Create the multi-row multi-species UJC component plot with data files."""

    all_species       = ALL_SPECIES_DEFAULT
    include_erp_anno  = panels in ('ERPanno', 'ERPnovel')
    include_erp_novel = panels == 'ERPnovel'

    if (include_erp_anno or include_erp_novel) and not data_gtf_files:
        raise ValueError(
            '--data_gtf_files is required when --panels is ERPanno or ERPnovel.')

    # ------------------------------------------------------------------
    # Load component list
    # ------------------------------------------------------------------
    comp_list_df = pd.read_csv(component_list_file)
    if 'component_id' not in comp_list_df.columns:
        raise ValueError("component_list file must contain 'component_id' column")

    if 'geneID' in comp_list_df.columns and 'source_species' in comp_list_df.columns:
        geneID         = comp_list_df['geneID'].iloc[0]
        source_species = comp_list_df['source_species'].iloc[0]
        component_panel_title = f'{geneID}:{source_species}'
    else:
        component_panel_title = 'Component ID'

    component_ids = comp_list_df['component_id'].tolist()
    components, comp_gene_ids = load_components_from_list(
        component_file, component_ids, all_species)

    # ------------------------------------------------------------------
    # Load annotated coordinates + read data per species
    # ------------------------------------------------------------------
    all_coords           = {}
    all_dominant_strands = {}
    all_raw_coords       = {}
    all_chrs             = {}
    all_strands          = {}
    all_props            = {}
    all_effect_sizes     = {}
    all_total_reads      = {}
    all_sumReads_f       = {}
    all_sumReads_m       = {}
    all_modal_genes      = {}
    all_anno_jxn_to_erp  = {}

    for sp in all_species:
        sp_jxns = [j for cd in components.values() for j in cd[sp]]

        if not sp_jxns:
            for d in (all_coords, all_raw_coords, all_chrs, all_strands,
                      all_props, all_effect_sizes, all_total_reads,
                      all_sumReads_f, all_sumReads_m):
                d[sp] = {}
            all_dominant_strands[sp] = None
            all_modal_genes[sp]      = None
            continue

        gtf = anno_gtf_files.get(sp)
        if not gtf:
            for d in (all_coords, all_raw_coords, all_chrs, all_strands,
                      all_props, all_effect_sizes, all_total_reads,
                      all_sumReads_f, all_sumReads_m):
                d[sp] = {}
            all_dominant_strands[sp]  = None
            all_modal_genes[sp]       = None
            all_anno_jxn_to_erp[sp]   = {}
            continue

        raw_coords, strands, chrs = load_gtf_coords_by_jxnhash(
            gtf, sp_jxns, return_chr=True)

        all_raw_coords[sp] = raw_coords
        all_chrs[sp]       = chrs
        all_strands[sp]    = strands

        sp_comp_genes = comp_gene_ids.get(sp, {})
        if sp_comp_genes:
            modal_gene     = max(set(sp_comp_genes.values()),
                                 key=list(sp_comp_genes.values()).count)
            target_strands = [strands[j] for j in sp_jxns
                              if j in strands
                              and sp_comp_genes.get(j) == modal_gene]
        else:
            modal_gene     = None
            target_strands = []

        if target_strands:
            dominant_strand = max(set(target_strands), key=target_strands.count)
        else:
            strand_vals     = list(strands.values())
            dominant_strand = (max(set(strand_vals), key=strand_vals.count)
                               if strand_vals else '+')

        all_coords[sp]           = flip_coords_to_display(raw_coords, dominant_strand)
        all_dominant_strands[sp] = dominant_strand
        all_modal_genes[sp]      = modal_gene

        datafile   = jxn_datafiles[sp]
        anno_jxn_df = load_jxn_datafile(datafile, sp, sp_jxns)
        props, es, total_reads, _, _, sumReads_f, sumReads_m = load_read_props_by_jxnhash(
            anno_jxn_df, sp, sp_jxns)
        all_props[sp]          = props
        all_effect_sizes[sp]   = es
        all_total_reads[sp]    = total_reads
        all_sumReads_f[sp]     = sumReads_f
        all_sumReads_m[sp]     = sumReads_m
        all_anno_jxn_to_erp[sp] = get_jxnhash_to_erp_map(anno_jxn_df, modal_gene) \
            if modal_gene else {}
        del anno_jxn_df

    # ------------------------------------------------------------------
    # Apply anno_min_reads filter
    # ------------------------------------------------------------------
    if anno_min_reads > 0:
        for sp in all_species:
            keep = {
                jxn: exons for jxn, exons in all_coords[sp].items()
                if all_total_reads[sp].get(jxn, 0) >= anno_min_reads
            }
            dropped = len(all_coords[sp]) - len(keep)
            if dropped:
                print(f'[anno] {sp}: dropped {dropped} annotated UJC(s) '
                      f'with < {anno_min_reads} reads.')
            all_coords[sp] = keep

    # ------------------------------------------------------------------
    # Component ordering
    # ------------------------------------------------------------------
    component_order = sort_components_by_erp(comp_list_df)

    master_positions, annotated_max_pos, component_positions = (
        create_component_alignment(
            components, all_species, all_coords,
            component_order=component_order))

    annotated_max_pos = max(annotated_max_pos, 5)

    # ------------------------------------------------------------------
    # Build ERP circle-marker cache per species
    # ------------------------------------------------------------------
    erp_es_by_jxn = {sp: {} for sp in all_species}
    erp_es_all    = {}

    if erp_datafiles:
        for sp in all_species:
            erp_datafile = erp_datafiles.get(sp)
            if erp_datafile and os.path.exists(erp_datafile):
                df_erp = pd.read_csv(
                    erp_datafile, low_memory=False,
                    usecols=lambda c: c in {'ERP_plus', 'effect_size_Equal'})
                if 'ERP_plus' in df_erp.columns and 'effect_size_Equal' in df_erp.columns:
                    erp_es_all[sp] = (
                        df_erp.dropna(subset=['effect_size_Equal'])
                               .set_index('ERP_plus')['effect_size_Equal']
                               .to_dict())
                del df_erp

        for sp in all_species:
            if not all_modal_genes.get(sp):
                continue
            if sp not in erp_es_all:
                print(f'[erp_marker] Warning: ERP datafile not found for {sp}')
                continue
            jxn_to_erp   = all_anno_jxn_to_erp.get(sp, {})
            erp_es_cache = erp_es_all[sp]
            erp_es_by_jxn[sp] = {
                jxn: erp_es_cache[erp]
                for jxn, erp in jxn_to_erp.items()
                if erp in erp_es_cache
            }
            n_sig = sum(1 for es in erp_es_by_jxn[sp].values() if abs(es) > 2.78)
            print(f'[erp_marker] {sp}: {len(erp_es_by_jxn[sp])} jxnHashes '
                  f'with ERP effect size, {n_sig} in significant ERP_plus.')

    # ------------------------------------------------------------------
    # Discover novel UJCs per species
    # ------------------------------------------------------------------
    erp_anno_coords     = {sp: {} for sp in all_species}
    erp_anno_props      = {sp: {} for sp in all_species}
    erp_anno_es         = {sp: {} for sp in all_species}
    erp_anno_erp_group  = {sp: {} for sp in all_species}
    erp_anno_jxn_lists  = {sp: [] for sp in all_species}
    erp_anno_sumReads_f = {sp: {} for sp in all_species}
    erp_anno_sumReads_m = {sp: {} for sp in all_species}
    erp_anno_sumReads   = {sp: {} for sp in all_species}

    erp_novel_coords    = {sp: {} for sp in all_species}
    erp_novel_props     = {sp: {} for sp in all_species}
    erp_novel_es        = {sp: {} for sp in all_species}
    erp_novel_erp_group = {sp: {} for sp in all_species}
    erp_novel_jxn_lists = {sp: [] for sp in all_species}
    erp_novel_sumReads_f = {sp: {} for sp in all_species}
    erp_novel_sumReads_m = {sp: {} for sp in all_species}
    erp_novel_sumReads   = {sp: {} for sp in all_species}

    if include_erp_anno or include_erp_novel:
        for sp in all_species:
            data_gtf = data_gtf_files.get(sp) if data_gtf_files else None

            if not data_gtf or not os.path.exists(data_gtf):
                print(f'[novel] Warning: data GTF not found for {sp}: {data_gtf}')
                continue

            raw_sp = all_raw_coords.get(sp, {})
            if not raw_sp:
                continue

            gene_strand = all_dominant_strands.get(sp)
            if not gene_strand:
                continue

            strand_filtered_raw = {
                jxn: exons for jxn, exons in raw_sp.items()
                if all_strands[sp].get(jxn) == gene_strand
            }
            gene_min, gene_max = get_gene_region_from_raw_coords(strand_filtered_raw)
            if gene_min is None:
                continue

            chr_vals = list(all_chrs[sp].values())
            if not chr_vals:
                continue
            gene_chr = max(set(chr_vals), key=chr_vals.count)

            known      = list(all_coords[sp].keys())
            candidates = load_unannotated_jxns_in_region(
                data_gtf, gene_chr, gene_min, gene_max, known,
                gene_strand=gene_strand)

            if not candidates:
                continue

            datafile     = jxn_datafiles[sp]
            novel_jxn_df = load_jxn_datafile(datafile, sp, candidates)
            cand_props, cand_es, cand_total, cand_erp_anno_flag, _, cand_f, cand_m = \
                load_read_props_by_jxnhash(novel_jxn_df, sp, candidates)

            # ---- ERPanno row ----
            if include_erp_anno:
                kept_anno = filter_novel(
                    candidates, cand_es, cand_total, cand_erp_anno_flag,
                    erp_anno_flag=True,
                    min_reads=erp_anno_min_reads,
                    sig_only=erp_anno_sig_only,
                    top_n=erp_anno_top_n,
                    sp=sp, label='ERPanno')

                if kept_anno:
                    raw_a, strands_a     = load_gtf_coords_by_jxnhash(data_gtf, kept_anno)
                    erp_anno_coords[sp]  = auto_flip_by_strand(raw_a, strands_a)
                    erp_anno_props[sp]   = {j: cand_props[j] for j in kept_anno if j in cand_props}
                    erp_anno_es[sp]      = {j: cand_es[j]    for j in kept_anno if j in cand_es}
                    erp_anno_jxn_lists[sp] = kept_anno
                    erp_anno_sumReads_f[sp] = {j: cand_f[j] for j in kept_anno if j in cand_f}
                    erp_anno_sumReads_m[sp] = {j: cand_m[j] for j in kept_anno if j in cand_m}
                    erp_anno_sumReads[sp] = {j: cand_total[j] for j in kept_anno if j in cand_total}
                    jxn_to_erp_a         = get_jxnhash_to_erp_map_for_jxns(novel_jxn_df, kept_anno)
                    erp_anno_erp_group[sp] = jxn_to_erp_a
                    if erp_datafiles and sp in erp_es_all:
                        erp_es_cache = erp_es_all[sp]
                        for jxn, erp in jxn_to_erp_a.items():
                            if erp in erp_es_cache:
                                erp_es_by_jxn[sp][jxn] = erp_es_cache[erp]

            # ---- ERPnovel row ----
            if include_erp_novel:
                kept_novel = filter_novel(
                    candidates, cand_es, cand_total, cand_erp_anno_flag,
                    erp_anno_flag=False,
                    min_reads=erp_novel_min_reads,
                    sig_only=erp_novel_sig_only,
                    top_n=erp_novel_top_n,
                    sp=sp, label='ERPnovel')

                if kept_novel:
                    raw_n, strands_n      = load_gtf_coords_by_jxnhash(data_gtf, kept_novel)
                    erp_novel_coords[sp]  = auto_flip_by_strand(raw_n, strands_n)
                    erp_novel_props[sp]   = {j: cand_props[j] for j in kept_novel if j in cand_props}
                    erp_novel_es[sp]      = {j: cand_es[j]    for j in kept_novel if j in cand_es}
                    erp_novel_jxn_lists[sp] = kept_novel
                    erp_novel_sumReads_f[sp] = {j: cand_f[j] for j in kept_novel if j in cand_f}
                    erp_novel_sumReads_m[sp] = {j: cand_m[j] for j in kept_novel if j in cand_m}
                    erp_novel_sumReads[sp] = {j: cand_total[j] for j in kept_novel if j in cand_total}
                    jxn_to_erp_n          = get_jxnhash_to_erp_map_for_jxns(novel_jxn_df, kept_novel)
                    erp_novel_erp_group[sp] = jxn_to_erp_n
                    if erp_datafiles and sp in erp_es_all:
                        erp_es_cache = erp_es_all[sp]
                        for jxn, erp in jxn_to_erp_n.items():
                            if erp in erp_es_cache:
                                erp_es_by_jxn[sp][jxn] = erp_es_cache[erp]

            del novel_jxn_df

    # ------------------------------------------------------------------
    # Compute shared x-limits per species (annotated → ERPanno → ERPnovel)
    # ------------------------------------------------------------------
    shared_xlim = {}
    for sp in all_species:
        if not all_dominant_strands.get(sp):
            shared_xlim[sp] = None
            continue
        if clamp_xlim:
            shared_xlim[sp] = compute_shared_xlim([all_coords[sp]])
        else:
            combined = {**all_coords[sp], **erp_anno_coords[sp], **erp_novel_coords[sp]}
            shared_xlim[sp] = compute_shared_xlim([combined])

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    n_erp_anno_rows  = max((len(c) for c in erp_anno_coords.values()),  default=0)
    n_erp_novel_rows = max((len(c) for c in erp_novel_coords.values()), default=0)
    n_erp_anno_rows  = max(n_erp_anno_rows,  5) if include_erp_anno  else 0
    n_erp_novel_rows = max(n_erp_novel_rows, 5) if include_erp_novel else 0

    height_ratios = [annotated_max_pos]
    fig_height    = annotated_max_pos * 0.3
    n_rows        = 1

    if include_erp_anno:
        height_ratios.append(n_erp_anno_rows)
        fig_height += n_erp_anno_rows * 0.3
        n_rows += 1

    if include_erp_novel:
        height_ratios.append(n_erp_novel_rows)
        fig_height += n_erp_novel_rows * 0.3
        n_rows += 1

    fig_height = max(fig_height, 4)

    widths = [0.4, 1, 1, 1, 1, 1]
    fig    = plt.figure(figsize=(24, fig_height))

    gs_kw = dict(width_ratios=widths, left=0.08, right=0.95, wspace=0.25)
    if n_rows > 1:
        gs_kw.update(height_ratios=height_ratios,
                     top=0.92, bottom=0.08, hspace=0.3)
    else:
        gs_kw.update(top=0.85, bottom=0.15)

    gs = fig.add_gridspec(n_rows, 6, **gs_kw)

    ann_axes       = [fig.add_subplot(gs[0, i]) for i in range(6)]
    erp_anno_axes  = [fig.add_subplot(gs[1, i]) for i in range(6)] \
                     if include_erp_anno else []
    erp_novel_axes = [fig.add_subplot(gs[2, i]) for i in range(6)] \
                     if include_erp_novel else []

    colormap, normalize = get_colormap()

    # ------------------------------------------------------------------
    # Draw annotated row
    # ------------------------------------------------------------------
    plot_component_panel(
        ann_axes[0], component_positions, annotated_max_pos,
        title=component_panel_title)

    for idx, sp in enumerate(all_species, start=1):
        show_x = not (include_erp_anno or include_erp_novel)
        plot_annotated_panel(
            ann_axes[idx],
            all_coords[sp], master_positions[sp],
            all_props[sp], all_effect_sizes[sp],
            SPECIES_NAMES[sp], annotated_max_pos,
            colormap, normalize,
            show_xaxis=show_x,
            forced_xlim=shared_xlim.get(sp),
        )

    # ------------------------------------------------------------------
    # Draw ERPanno row
    # ------------------------------------------------------------------
    if include_erp_anno and erp_anno_axes:
        plot_novel_label_panel(
            erp_anno_axes[0], n_erp_anno_rows, label='ERPanno UJCs')

        for idx, sp in enumerate(all_species, start=1):
            show_x = not include_erp_novel
            plot_novel_panel(
                erp_anno_axes[idx],
                erp_anno_coords[sp],
                erp_anno_props[sp],
                erp_anno_es[sp],
                SPECIES_NAMES[sp],
                colormap, normalize,
                erp_es_by_jxn=erp_es_by_jxn[sp],
                jxn_erp_group=erp_anno_erp_group[sp],
                panel_label='ERPanno',
                show_xaxis=show_x,
                forced_xlim=shared_xlim.get(sp),
                max_rows=n_erp_anno_rows,
            )

    # ------------------------------------------------------------------
    # Draw ERPnovel row
    # ------------------------------------------------------------------
    if include_erp_novel and erp_novel_axes:
        plot_novel_label_panel(
            erp_novel_axes[0], n_erp_novel_rows, label='ERPnovel UJCs')

        for idx, sp in enumerate(all_species, start=1):
            plot_novel_panel(
                erp_novel_axes[idx],
                erp_novel_coords[sp],
                erp_novel_props[sp],
                erp_novel_es[sp],
                SPECIES_NAMES[sp],
                colormap, normalize,
                erp_es_by_jxn=erp_es_by_jxn[sp],
                jxn_erp_group=erp_novel_erp_group[sp],
                panel_label='ERPnovel',
                show_xaxis=True,
                forced_xlim=shared_xlim.get(sp),
                max_rows=n_erp_novel_rows,
            )

    # ------------------------------------------------------------------
    # Colour bar
    # ------------------------------------------------------------------
    cbar_ax = fig.add_axes([0.96, 0.15, 0.01, 0.7])
    mappable = ScalarMappable(norm=normalize, cmap=colormap)
    mappable.set_array([])
    cbar = fig.colorbar(mappable, cax=cbar_ax)
    cbar.set_label(
        'Effect Size (Equal)\n'
        'Fill: −9 (Male) → +9 (Female)\n'
        '★ = significant UJC\n'
        '● = significant ERP_plus',
        fontsize=9)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    anno_n_str  = str(erp_anno_top_n)  if erp_anno_top_n  is not None else 'all'
    novel_n_str = str(erp_novel_top_n) if erp_novel_top_n is not None else 'all'

    if panels == 'UJCanno':
        out_file = (f'{output_dir}/{gene_name}'
                    f'_UJCanno_min_{anno_min_reads}_reads'
                    f'.{output_format}')
    elif panels == 'ERPanno':
        out_file = (f'{output_dir}/{gene_name}'
                    f'_ERPanno_annoMin_{anno_min_reads}'
                    f'_erpAnnoTop{anno_n_str}_min_{erp_anno_min_reads}_reads'
                    f'.{output_format}')
    else:  # ERPnovel
        out_file = (f'{output_dir}/{gene_name}'
                    f'_ERPnovel_annoMin_{anno_min_reads}'
                    f'_erpAnnoTop{anno_n_str}_min_{erp_anno_min_reads}_reads'
                    f'_erpNovelTop{novel_n_str}_min_{erp_novel_min_reads}_reads'
                    f'.{output_format}')

    os.makedirs(output_dir, exist_ok=True)
    save_kw = {} if output_format == 'svg' else {'dpi': 150}
    fig.savefig(out_file, bbox_inches='tight', **save_kw)
    print(f'saved: {out_file}')
    plt.close(fig)
    gc.collect()

    # ------------------------------------------------------------------
    # Save CSV lookup table
    # ------------------------------------------------------------------
    csv_rows = []
    
    # UJCanno row
    for sp in all_species:
        if not all_coords[sp]:
            continue
        jxn_list_ordered = sorted(
            master_positions[sp].keys(),
            key=lambda j: master_positions[sp][j]
        )
        sp_csv = build_index_csv(
            'UJCanno', sp, jxn_list_ordered,
            all_sumReads_f[sp], all_sumReads_m[sp], all_total_reads[sp]
        )
        csv_rows.append(sp_csv)
    
    # ERPanno row
    if include_erp_anno:
        for sp in all_species:
            if erp_anno_jxn_lists[sp]:
                sp_csv = build_index_csv(
                    'ERPanno', sp, erp_anno_jxn_lists[sp],
                    erp_anno_sumReads_f[sp], erp_anno_sumReads_m[sp], 
                    erp_anno_sumReads[sp]
                )
                csv_rows.append(sp_csv)
    
    # ERPnovel row
    if include_erp_novel:
        for sp in all_species:
            if erp_novel_jxn_lists[sp]:
                sp_csv = build_index_csv(
                    'ERPnovel', sp, erp_novel_jxn_lists[sp],
                    erp_novel_sumReads_f[sp], erp_novel_sumReads_m[sp], 
                    erp_novel_sumReads[sp]
                )
                csv_rows.append(sp_csv)
    
    if csv_rows:
        csv_df = pd.concat(csv_rows, ignore_index=True)
        csv_file = out_file.replace(f'.{output_format}', '_index_lookup.csv')
        csv_df.to_csv(csv_file, index=False)
        print(f'saved: {csv_file}')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Multi-row multi-species UJC component plot with data files — gene mode'
    )

    # ---- required ----
    parser.add_argument('--component_file', required=True)
    parser.add_argument('--component_list', required=True)
    parser.add_argument('--gene_name',      required=True)
    parser.add_argument('--output_dir',     required=True)

    parser.add_argument('--anno_gtf_files', nargs='+', required=True,
                        metavar='SPECIES:PATH')
    parser.add_argument('--jxn_datafiles',  nargs='+', required=True,
                        metavar='SPECIES:PATH')
    parser.add_argument('--data_gtf_files', nargs='+', default=None,
                        metavar='SPECIES:PATH',
                        help='Required for --panels ERPanno or ERPnovel.')

    # ---- panel selection ----
    parser.add_argument('--panels', choices=['UJCanno', 'ERPanno', 'ERPnovel'],
                        default='UJCanno')

    # ---- annotated row ----
    parser.add_argument('--anno_min_reads', type=int, default=1)

    # ---- ERPanno row ----
    parser.add_argument('--erp_anno_min_reads', type=int, default=1)
    parser.add_argument('--erp_anno_top_n',     type=int, default=None)
    parser.add_argument('--erp_anno_sig_only',  action='store_true')

    # ---- ERPnovel row ----
    parser.add_argument('--erp_novel_min_reads', type=int, default=1)
    parser.add_argument('--erp_novel_top_n',     type=int, default=None)
    parser.add_argument('--erp_novel_sig_only',  action='store_true')

    # ---- ERP circle markers ----
    parser.add_argument('--erp_datafiles', nargs='+', default=None,
                        metavar='SPECIES:PATH')

    # ---- output ----
    parser.add_argument('--output_format', choices=['svg', 'png', 'pdf'], default='png')
    parser.add_argument('--clamp_xlim', action='store_true',
                        help='Set xlim from annotated UJCs only; default uses all panels.')

    args = parser.parse_args()

    plot_gene_data(
        component_file      = args.component_file,
        component_list_file = args.component_list,
        anno_gtf_files      = parse_species_paths(args.anno_gtf_files),
        jxn_datafiles       = parse_species_paths(args.jxn_datafiles),
        data_gtf_files      = parse_species_paths(args.data_gtf_files)
                              if args.data_gtf_files else None,
        gene_name           = args.gene_name,
        output_dir          = args.output_dir,
        panels              = args.panels,
        anno_min_reads      = args.anno_min_reads,
        erp_anno_min_reads  = args.erp_anno_min_reads,
        erp_novel_min_reads = args.erp_novel_min_reads,
        erp_anno_top_n      = args.erp_anno_top_n,
        erp_novel_top_n     = args.erp_novel_top_n,
        erp_anno_sig_only   = args.erp_anno_sig_only,
        erp_novel_sig_only  = args.erp_novel_sig_only,
        erp_datafiles       = parse_species_paths(args.erp_datafiles)
                              if args.erp_datafiles else None,
        output_format       = args.output_format,
        clamp_xlim          = args.clamp_xlim,
    )