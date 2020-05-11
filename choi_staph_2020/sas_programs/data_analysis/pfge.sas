
data pfge;
input ID	PatientNumber	PFGE $;
datalines;
4132	0	Different
4130	0	Different
4117	1	Different
4388	1	Identical
4402	1	Identical
5472	4	Identical
5541	4	Identical
4264	5	Identical
4365	5	Identical
4616	7	Different
4912	7	Different
4218	9	Different
4268	9	Different
3979	11	Identical
4124	11	Identical
4471	15	Identical
4590	15	Identical
3867	16	identical
3938	16	identical
4024	16	identical
4725	19	different
5098	19	Identical
5173	19	Identical
4052	20	Identical
4160	20	Identical
5852	21	Identical
5901	21	Identical
3758	22	Identical
3809	22	Identical
3907	23	Identical
4026	23	Identical
3978	24	Identical
4020	24	Identical
3800	25	Identical
3922	25	Identical
3765	27	Different
3825	27	Different
5321	28	Different
5630	28	Different
3984	29	Identical
4167	29	Different
4519	29	Identical
5132	30	Different
5462	30	Different
6018	31	Identical
6040	31	Identical
5587	32	Identical
5639	32	Identical
4002	34	Identical
4053	34	Identical
4739	35	Identical
4697	35	Identical
4885	36	Different
4977	36	Different
5394	37	Identical
5463	37	Identical
3930	38	Identical
3998	38	Identical
3843	39	Different
3942	39	Different
4187	40	Different
4480	40	Identical
4974	40	Identical
4824	41	Identical
4846	41	Identical
3855	42	Different
5434	42	Identical
5744	42	Identical
5019	43	Identical
5261	43	Identical
4359	44	Different
4838	44	Different
3746	45	Different
4168	45	Different
4355	46	Identical
4372	46	Identical
4307	47	Different
4789	47	Different
4503	48	Different
4757	48	Different
4451	49	Different
4533	49	Different
5700	50	Identical
5794	50	Identical
4255	51	Identical
4649	51	Identical
5682	52	Different
5780	52	Different
3869	53	Identical
4025	53	Identical
4159	54	Identical
4557	54	Identical
4834	54	Identical
4572	55	Different
5873	55	Different
4530	56	Identical
4696	56	Identical
4523	57	Different
4777	57	Different
3803	58	Identical
4702	58	Identical
3928	59	Different
4072	59	Different
5926	60	Identical
5940	60	Identical
4278	61	Identical
4426	61	Identical
4643	61	Identical
5635	62	Different
5824	62	Different
5107	63	Identical
5250	63	Identical
5254	64	Identical
5767	64	Identical
3834	65	23-24 ID
3850	65	23-24 ID
3899	65	25-26 ID
4122	65	25-26 ID
4279	66	Identical
5225	66	Identical
4447	67	Identical
4543	67	Identical
4012	68	Identical
4021	68	Identical
4886	69	Identical
4910	69	Identical
5040	69	Identical
5051	70	Identical
5692	70	Identical
5067	71	unknown
5102	71	unknown
4066	72	Different
5468	72	Different
4500	73	Identical
4553	73	Identical
4695	73	Identical
4732	73	Identical
5157	74	Different
5411	74	Different
5743	74	Different
5376	76	Different
5574	76	Different
5236	77	Identical
5354	77	Identical
5820	79	Different
5974	79	Different
3739	80	Identical
3775	80	Identical
5029	81	Identical
5091	81	Identical
3815	82	Identical
3847	82	Identical
3932	82	Identical
3740	83	Identical
3833	83	Identical
3838	83	Identical
4554	83	Identical
3798	84	Different
4201	84	Different
5524	86	Different
5673	86	Different
;
run;

proc freq data=pfge noprint;
by patientNumber;
tables PFGE/out=count_patient;
run;

proc freq data=count_patient;
tables patientnumber/out=count_multiple;
run;

data list_multiple;
set count_multiple;
where COUNT >1;
run;

proc sort data=list_multiple;
by patientnumber;

proc sort data=pfge;
by patientnumber;


data multiple single;
merge pfge (in=in1) list_multiple(in=in2);
by patientnumber;
if in2 then output multiple;
else output single;

run;

proc freq data=single;
tables patientnumber*PFGE/out=count_pfge;
run;

proc freq data=count_pfge;
tables pfge;
run;
