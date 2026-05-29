# -*- coding: utf-8 -*-


import pandas as pd
import os

ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

infile = os.path.join(ind,"allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv")

in_df = pd.read_csv(infile)

out_df = in_df[['ef_id','ef_start']]

out_df = out_df.rename(columns={"ef_id": "featureID", "ef_start": "start"})

out_df.to_csv(os.path.join(ind,"featureID_2_ef_start.csv"), sep=",", index=False)