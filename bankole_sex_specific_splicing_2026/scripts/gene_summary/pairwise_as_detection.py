#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 19:39:01 2026

@author: mgaran
"""

def erp_to_set(erp):
    # strip strand prefix (first character), keep index of every '1'
    stripped = erp[2:]
    return set(i for i, c in enumerate(stripped) if c == '1')

def pairwise_as_detection(t1, t2):
    flag_alt_dataOnlyExon = 0
    flag_alt_IR = 0
    flag_alt_donorAcceptor = 0
    flag_alt_5pER = 0
    flag_alt_3pER = 0
    flag_alt_skipER = 0

    # 1) variable comparisons
    if t1['dataOnlyER_ID'] != t2['dataOnlyER_ID']:
        flag_alt_dataOnlyExon = 1

    if t1['IR_ER'] != t2['IR_ER']:
        flag_alt_IR = 1

    if (t1['ERP'] == t2['ERP']
            and t1['dataOnlyER_ID'] == t2['dataOnlyER_ID']
            and t1['IR_ER'] == t2['IR_ER']):
        flag_alt_donorAcceptor = 1

    # 2) set operations
    if t1['ERP'] == t2['ERP']:
        flag_alt_skipER = 0
        flag_alt_5pER = 0
        flag_alt_3pER = 0
    else:
        erp_1_set = erp_to_set(t1['ERP'])
        erp_2_set = erp_to_set(t2['ERP'])

        lower = max(min(erp_1_set), min(erp_2_set))
        upper = min(max(erp_1_set), max(erp_2_set))
        s_int = set(p for p in (erp_1_set | erp_2_set) if lower < p < upper)

        # take the symmetric difference between set 1 and set 2
        differing = (erp_1_set & s_int) ^ (erp_2_set & s_int)
        
        if len(differing) != 0:
            flag_alt_skipER = 1
        if min(erp_1_set) != min(erp_2_set):
            flag_alt_5pER = 1
        if max(erp_1_set) != max(erp_2_set):
            flag_alt_3pER = 1

    return {
        'flag_alt_dataOnlyExon': flag_alt_dataOnlyExon,
        'flag_alt_IR': flag_alt_IR,
        'flag_alt_donorAcceptor': flag_alt_donorAcceptor,
        'flag_alt_5pER': flag_alt_5pER,
        'flag_alt_3pER': flag_alt_3pER,
        'flag_alt_skipER': flag_alt_skipER,
    }