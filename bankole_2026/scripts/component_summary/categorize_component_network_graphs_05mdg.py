import pandas as pd
from collections import Counter

EDGE_PATH   = '/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/link_files/edges.csv'
NODE_PATH   = '/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/link_files/component_map_by_node.csv'
ANNOT_DIR   = '/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary'
OUTFILE    = '/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/component_network_graph_categories.csv'

species_list = ['dmel6', 'dsim2', 'dser1', 'dsan1', 'dyak2']

# the four non-serrata species, used to identify the4_no_ser components
THE4 = ['dmel6', 'dsim2', 'dsan1', 'dyak2']

def calculate_flags(component_nodes, edge_set):
    """
    calculates flag_1_node_per_species by looping over species node dictionary
    items
    calculates flag_all_linked by looking over set of tuples where each pair of
    nodes (UJC) is conencted with an edge
    """
    # creates dictionary with node counts (UJC) for each species per component
    species_counts = Counter(component_nodes['source'])

    # create the set of species with at least one node
    present = {species for species, count in species_counts.items() if count > 0}

    # flag_1_node_per_sp: 1 if every species that is present has exactly one node, 0 otherwise
    # this is used to identify one-to-one ortholog relationships
    flag_1_node_per_sp = int(all(species_counts[sp] == 1 for sp in present))

    # build a list of (jxnHash, species) tuples for every node in the component
    nodes = list(component_nodes[['jxnHash', 'source']].itertuples(index=False, name=None))

    # check whether every pair of nodes from different species is connected by an edge
    # this determines whether the component is fully linked across species
    all_linked = True
    for first_node_idx in range(len(nodes)):
        source_jxnHash, source_species = nodes[first_node_idx]
        for second_node_idx in range(first_node_idx + 1, len(nodes)):
            target_jxnHash, target_species = nodes[second_node_idx]
            # skip pairs from the same species because edges only connect across species
            if source_species == target_species:
                continue
            # checks if either source_jxnHash:target_jxnHash OR target_jxnHash:source_jxnHash is in the edge set because links are independent of direction
            if frozenset((source_jxnHash, target_jxnHash)) not in edge_set:
                all_linked = False
                break
        if not all_linked:
            break

    # flag_all_linked: 1 if all cross-species node pairs are linked, 0 otherwise
    flag_all_linked = int(all_linked)
    return species_counts, present, flag_1_node_per_sp, flag_all_linked

def classify_species(species_counts):
    """
    calculate "species" column by checking the species with at least 1 node
    (UJC) present in the component.
    returns the values "all5", "the4_no_ser", "mel_sim_only", "san_yak_only",
    each of the individual species names if only 1 present, and "other" if
    anything else.
    """
    # present is the set of species that have at least one node in this component
    present = {species for species, count in species_counts.items() if count > 0}    # species_counts is a dictionary maping node counts to each species in a component

    # all five species represented
    if all(sp in present for sp in species_list) and len(present) == 5:
        return 'all5'

    # exactly the four non-serrata species, serrata absent
    if all(sp in present for sp in THE4) and 'dser1' not in present and len(present) == 4:
        return 'the4_no_ser'

    # only mel and sim nodes present
    if present == {'dmel6', 'dsim2'}:
        return 'mel_sim_only'

    # only san and yak nodes present
    if present == {'dsan1', 'dyak2'}:
        return 'san_yak_only'

    # exactly one species present
    if len(present) == 1:
        species = next(iter(present))
        if species == 'dmel6':
            return 'mel'
        if species == 'dsim2':
            return 'sim'
        if species == 'dser1':
            return 'ser'
        if species == 'dsan1':
            return 'san'
        if species == 'dyak2':
            return 'yak'

    # anything that does not match a named pattern
    return 'other'

def is_complete_clique(jxns, edge_set):
    # returns True if every pair of jxnHashes in the list shares an edge
    # used to check whether a subset of nodes forms a complete subgraph
    for first_node_idx in range(len(jxns)):
        for second_node_idx in range(first_node_idx + 1, len(jxns)):
            if frozenset((jxns[first_node_idx], jxns[second_node_idx])) not in edge_set:
                return False
    return True

def detect_pendant_ser_patterns(node_df, edge_set, species_counts):
    # this function handles a specific topology category where:
    # - all 5 species are present with exactly 1 node each
    # - the 4 non-serrata nodes form a complete clique (k4)
    # - serrata connects to only a subset of those 4 nodes
    # the return value describes which subset serrata connects to

    # build a lookup from species label to its single jxnHash in this component
    sp_to_jxn = {}
    for jxn, species in node_df[['jxnHash', 'source']].itertuples(index=False, name=None):
        sp_to_jxn[species] = jxn

    # confirm that every species has exactly 1 node before continuing
    if not all(species_counts.get(sp, 0) == 1 for sp in species_list):
        return None

    # pull the single jxnHash for each species
    mel = sp_to_jxn.get('dmel6')
    sim = sp_to_jxn.get('dsim2')
    ser = sp_to_jxn.get('dser1')
    san = sp_to_jxn.get('dsan1')
    yak = sp_to_jxn.get('dyak2')

    # if any species is missing from the lookup, something is wrong, return None
    if None in (mel, sim, ser, san, yak):
        return None

    # check that the four non-serrata nodes are all mutually connected
    k4_nodes = [mel, sim, san, yak]
    if not is_complete_clique(k4_nodes, edge_set):
        return None

    # collect all nodes that serrata shares an edge with
    neighbors = set()
    for edge in edge_set:
        if ser in edge:
            source_jxnHash, target_jxnHash = tuple(edge)
            neighbors.add(target_jxnHash if source_jxnHash == ser else source_jxnHash)

    # match the neighbor set to a named pattern
    target_mel_sim = {mel, sim}
    target_san_yak = {san, yak}

    if neighbors == target_mel_sim:
        return 'k4_ser_to_mel_sim'
    if neighbors == target_san_yak:
        return 'k4_ser_to_san_yak'

    # single-species connections
    if neighbors == {mel}:
        return 'k4_ser_to_mel'
    if neighbors == {sim}:
        return 'k4_ser_to_sim'
    if neighbors == {san}:
        return 'k4_ser_to_san'
    if neighbors == {yak}:
        return 'k4_ser_to_yak'

    # serrata connects to some other combination that is not a named pattern
    return None

def classify_topology(species_label, flag_1_node_per_sp, flag_all_linked, species_counts,
                      present, node_df, edge_set):
    # assign a topology label to a component based on its species composition and edge structure

    # k5: all 5 species present, one node each, all cross-species pairs linked
    # this is a complete graph on 5 nodes
    if species_label == 'all5' and flag_1_node_per_sp == 1 and flag_all_linked == 1:
        return 'k5'

    # k4: the 4 non-serrata species present, one node each, all cross-species pairs linked
    # this is a complete graph on 4 nodes with serrata absent
    if species_label == 'the4_no_ser' and flag_1_node_per_sp == 1 and flag_all_linked == 1:
        return 'k4'

    # k4_ser_to_* patterns: all 5 species, one node each, but not fully linked
    # the k4 subgraph is complete and serrata connects to only a named subset
    if species_label == 'all5' and flag_1_node_per_sp == 1 and flag_all_linked == 0:
        k4_ser_top = detect_pendant_ser_patterns(node_df, edge_set, species_counts)
        if k4_ser_top is not None:
            return k4_ser_top

    # anything that does not match a named topology
    return 'other'

def classify_components(nodes_by_component, edges_by_component, suffix):
    # iterate over every component and compute its species, flag, and topology columns
    # suffix is either None (full annotation pass) or 'data' (data-filtered pass)
    # when suffix is provided the output column names get a trailing _<suffix>
    results = []
    total = len(nodes_by_component)
    for component_idx, (comp_id, comp_nodes) in enumerate(nodes_by_component, 1):
        if component_idx % 10000 == 0:
            print(f"  Processing component {component_idx}/{total}...")

        # get the edge rows that belong to this component, or an empty frame if none
        if comp_id in edges_by_component.groups:
            comp_edges_df = edges_by_component.get_group(comp_id)
        else:
            comp_edges_df = pd.DataFrame()

        # represent edges as a set of frozensets for fast membership testing
        edge_set = set(
            frozenset((row.source_jxnHash, row.target_jxnHash))
            for row in comp_edges_df.itertuples()
        )

        # compute per-component summary values
        species_counts, present, flag_1_node_per_sp, flag_all_linked = calculate_flags(comp_nodes, edge_set)
        species_label = classify_species(species_counts)
        topology = classify_topology(
            species_label, flag_1_node_per_sp, flag_all_linked,
            species_counts, present, comp_nodes, edge_set
        )

        # build the output row, appending suffix to column names if provided
        s = f"_{suffix}" if suffix else ""
        results.append({
            'componentID': comp_id,
            f'species{s}': species_label,
            f'flag_1_node_per_species{s}': flag_1_node_per_sp,
            f'flag_all_linked{s}': flag_all_linked,
            f'topology{s}': topology,
            f'node_count{s}': len(comp_nodes),
            f'edge_count{s}': len(edge_set)
        })
    return pd.DataFrame(results)

def main():
    print("Loading data...")
    edges         = pd.read_csv(EDGE_PATH, low_memory=False)
    component_map = pd.read_csv(NODE_PATH, low_memory=False)

    # standardize column names for consistency
    component_map = component_map.rename(columns={'jxnhash': 'jxnHash',
                                                  'component_id': 'componentID'})

    # build the set of annotated jxnHashes that are in the data
    # jxnHashes not in this set are annotation-only and will be excluded in the data column calculations
    print("Loading annotation tables and building data jxnHash set...")
    data_jxns = set()
    for species in species_list:
        annot_path  = f'{ANNOT_DIR}/fiveSpecies_{species}_full_annotation_w_component.csv'
        jxnhash_col = f'{species}_jxnHash'
        # only load the two columns needed to identify data-supported jxnHashes
        annot = pd.read_csv(annot_path, usecols=[jxnhash_col, 'flag_jxnHash_in_data'])
        sp_jxns = set(annot.loc[annot['flag_jxnHash_in_data'] == 1, jxnhash_col])
        print(f"  {species}: {len(sp_jxns):,} jxnHashes with data")
        data_jxns.update(sp_jxns)
    print(f"  Total across all species: {len(data_jxns):,} jxnHashes with data")

    # map components to edges and filter edges that connecet nodes (UJC) with different componentID
    print("Pre-filtering intra-component edges...")
    jxn_to_component = dict(zip(component_map['jxnHash'], component_map['componentID']))
    edges['source_componentID'] = edges['source_jxnHash'].map(jxn_to_component)
    edges['target_componentID'] = edges['target_jxnHash'].map(jxn_to_component)
    edges = edges[edges['source_componentID'] == edges['target_componentID']].copy()

    # classify every component using all annotation nodes and edges
    print("Running full classification (all nodes)...")
    full_df = classify_components(
        component_map.groupby('componentID'), edges.groupby('source_componentID'), suffix=None
    )

    # restrict nodes and edges to only annotated UJC with data
    print("Running data classification (data nodes only)...")
    component_map_data = component_map[component_map['jxnHash'].isin(data_jxns)].copy()
    edges_data = edges[
        edges['source_jxnHash'].isin(data_jxns) &
        edges['target_jxnHash'].isin(data_jxns)
    ]
    data_df = classify_components(
        component_map_data.groupby('componentID'), edges_data.groupby('source_componentID'), suffix='data'
    )

    # Outer merge
    print("Merging full and data classifications...")
    results_df = full_df.merge(data_df, on='componentID', how='outer', indicator='merge')
    
    # mark every component that had at least one data-supported node
    results_df['flag_has_data'] = (results_df['merge'] != 'left_only').astype(int)
    results_df.drop(columns=['merge'], inplace=True)

    # Fill NaNs with defaults ('' if string and 0 if float or int)
    for col in data_df.columns:
        if col != 'componentID':
            fill_val = '' if results_df[col].dtype == object else 0
            results_df[col] = results_df[col].fillna(fill_val)

    print(f"Saving results to {OUTFILE}")
    results_df.to_csv(OUTFILE, index=False)

if __name__ == "__main__":
    main()

# Full topology summary
# topology
# other                38901
# k5                    8801
# k4                    3560
# k4_ser_to_mel_sim     1160
# k4_ser_to_san_yak     1133
# k4_ser_to_sim         1023
# k4_ser_to_mel          833
# k4_ser_to_yak          617
# k4_ser_to_san          342
# Name: count, dtype: int64
# 
# Data topology summary
# topology_data
#                      28494
# other                21720
# k5                    3802
# k4                    1891
# k4_ser_to_mel_sim      159
# k4_ser_to_san_yak      139
# k4_ser_to_sim           62
# k4_ser_to_mel           49
# k4_ser_to_san           32
# k4_ser_to_yak           22
# Name: count, dtype: int64
# 
# Components with data    : 27,876
# Components without data : 28,494
# Total components         : 56,370