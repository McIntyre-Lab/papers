#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 13:47:40 2025

@author: mgaran
"""

import pandas as pd
import argparse

def find_components_with_misassignments(component_file, target_gene, source_species, output_file, 
                                       add_erp_from_species=None, anno_file=None):
    
    all_species = ['dmel6', 'dsim2', 'dyak2', 'dsan1', 'dser1']
    
    if source_species not in all_species:
        print(f"ERROR: Invalid source species {source_species}. Must be one of {all_species}")
        return
    
    print(f"\n{'='*70}")
    print(f"Finding components for gene: {target_gene} (source: {source_species})")
    print(f"{'='*70}")
    
    df = pd.read_csv(component_file, low_memory=False)
    df = df.rename(columns={'jxnhash':'jxnHash'})
    
    print(f"\n=== STEP 1: Finding primary components ===")
    
    # Find all components containing the target gene in the source species
    primary_components = df[
        (df['source'] == source_species) & (df['geneID'] == target_gene)
    ]['component_id'].unique()
    
    if len(primary_components) == 0:
        print(f"ERROR: No components found for {source_species} gene {target_gene}")
        return
    
    print(f"Found {len(primary_components)} primary components")
    
    # Get all geneIDs from other 4 species in these primary components
    other_species = [sp for sp in all_species if sp != source_species]
    primary_df = df[df['component_id'].isin(primary_components)]
    
    ortholog_genes = {}
    for sp in other_species:
        sp_genes = primary_df[primary_df['source'] == sp]['geneID'].unique()
        sp_genes = [g for g in sp_genes if pd.notna(g)]
        ortholog_genes[sp] = set(sp_genes)
        if sp_genes:
            print(f"  {sp}: {len(sp_genes)} ortholog gene(s)")
    
    # Find additional components with misassigned source species UJCs
    print(f"\n=== STEP 2: Searching for misassigned {source_species} UJCs ===")
    
    misassigned_components = set()
    misassigned_details = []
    
    for sp in other_species:
        if not ortholog_genes[sp]:
            continue
        
        # Find components containing these ortholog genes
        sp_comps = df[
            (df['source'] == sp) & (df['geneID'].isin(ortholog_genes[sp]))
        ]['component_id'].unique()
        
        # Check which have source species with DIFFERENT geneID
        for comp_id in sp_comps:
            if comp_id in primary_components:
                continue
            
            comp_data = df[df['component_id'] == comp_id]
            source_genes = comp_data[comp_data['source'] == source_species]['geneID'].unique()
            source_genes = [g for g in source_genes if pd.notna(g)]
            
            if source_genes and target_gene not in source_genes:
                if comp_id not in misassigned_components:
                    misassigned_components.add(comp_id)
                    sp_genes_in_comp = comp_data[comp_data['source'] == sp]['geneID'].unique()
                    misassigned_details.append({
                        'component_id': comp_id,
                        f'{source_species}_misassigned_to': ', '.join(source_genes),
                        'evidence_species': sp,
                        'evidence_gene': ', '.join(sp_genes_in_comp)
                    })
    
    if misassigned_components:
        print(f"Found {len(misassigned_components)} component(s) with misassigned {source_species} UJCs:")
        for detail in misassigned_details:
            print(f"  Component {detail['component_id']}: {source_species} assigned to {detail[f'{source_species}_misassigned_to']}")
    else:
        print(f"No misassigned components detected")
    
    all_components = sorted(list(primary_components) + list(misassigned_components))
    
    print(f"\n=== SUMMARY ===")
    print(f"Primary components: {len(primary_components)}")
    print(f"Misassigned components: {len(misassigned_components)}")
    print(f"Total components: {len(all_components)}")
    
    output_df = pd.DataFrame({
        'geneID': [target_gene] * len(all_components),
        'source_species': [source_species] * len(all_components),
        'component_id': all_components,
        'is_primary': [c in primary_components for c in all_components]
    })
    
    # Add ERP column if requested
    if add_erp_from_species and anno_file:
        print(f"\n=== STEP 3: Adding ERP from {add_erp_from_species} ===")
        
        if add_erp_from_species not in all_species:
            print(f"ERROR: Invalid species {add_erp_from_species}. Must be one of {all_species}")
            return
        
        comp_df = df[df['component_id'].isin(all_components)]
        
        try:
            anno_df = pd.read_csv(anno_file, 
                                 usecols=[f'{add_erp_from_species}_jxnHash', 'ERP'], 
                                 low_memory=False)
            anno_df = anno_df.rename(columns={f'{add_erp_from_species}_jxnHash': 'jxnHash'})
            
            erp_mapping = dict(zip(anno_df['jxnHash'], anno_df['ERP']))
            
            component_erps = {}
            for comp_id in all_components:
                comp_jxns = comp_df[
                    (comp_df['component_id'] == comp_id) & 
                    (comp_df['source'] == add_erp_from_species)
                ]['jxnHash'].values
                
                erp_values = [erp_mapping[jxn] for jxn in comp_jxns if jxn in erp_mapping]
                erp_value = min(erp_values) if erp_values else None
                
                component_erps[comp_id] = erp_value
            
            output_df['ERP'] = output_df['component_id'].map(component_erps)
            
            erp_count = output_df['ERP'].notna().sum()
            print(f"Added ERP values for {erp_count}/{len(output_df)} components")
            
        except FileNotFoundError:
            print(f"ERROR: Could not find annotation file: {anno_file}")
            print("Saving output without ERP column")
        except Exception as e:
            print(f"ERROR adding ERP: {e}")
            print("Saving output without ERP column")
    
    output_df.to_csv(output_file, index=False)
    print(f"\nSaved component list to: {output_file}")
    if 'ERP' in output_df.columns:
        print(f"  Component list includes ERP column from {add_erp_from_species}")
    
    return output_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Find all components for a gene from any species, including misassigned UJCs'
    )
    
    parser.add_argument('--component_file', required=True, 
                       help='Path to component_map_by_node.csv')
    parser.add_argument('--target_gene', required=True, 
                       help='Target geneID (e.g., FBgn0264270)')
    parser.add_argument('--source_species', required=True,
                       choices=['dmel6', 'dsim2', 'dyak2', 'dsan1', 'dser1'],
                       help='Source species for the target gene')
    parser.add_argument('--output', required=True,
                       help='Output file path for component list CSV')
    parser.add_argument('--add_erp', default=None,
                       choices=['dmel6', 'dsim2', 'dyak2', 'dsan1', 'dser1'],
                       help='Optional: Add ERP column from specified species annotation file')
    parser.add_argument('--anno_file', default=None,
                       help='annotation file of species with ERP column (required if --add_erp is used)')
    
    args = parser.parse_args()
    
    if args.add_erp and not args.anno_file:
        parser.error("--anno_file is required when --add_erp is specified")
    
    find_components_with_misassignments(
        component_file=args.component_file,
        target_gene=args.target_gene,
        source_species=args.source_species,
        output_file=args.output,
        add_erp_from_species=args.add_erp,
        anno_file=args.anno_file
    )