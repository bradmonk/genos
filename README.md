## genos

The Alzheimer's Disease Sequencing Project ([ADSP](https://www.niagads.org/adsp/content/home)) has made available WES data in vcf format. The following is a step-by-step analysis walkthrough of this dataset. Here the ultimate goal is to develop a platform for Alzheimer's Disease diagnosis based on genome sequencing information. The data for this analysis can be downloaded here:

* [ADSP_WES_VCF_LATEST_RELEASE.mat](https://drive.google.com/open?id=1e3tIbQhcDUF1vofAwf4oYVOl8XriC44v).

Both the dataset and the following analysis code are in MATLAB format. Note however there is an [R package](https://cran.r-project.org/web/packages/R.matlab/index.html) for reading .mat files. All MATLAB code and custom functions used in the demo walkthrough below can be downloaded on github: 

* [GENOS GITHUB CODE REPO](https://github.com/subroutines/genos) (aka, here)

A comprehensive genomics analysis and step-by-step code walkthrough (tailored towards machine learning for AD diagnosis) of this dataset can be found on my wiki: 

* http://bradleymonk.com/ADSP

This wiki hosts the most up-to-date information about these resources, and explains the code in this repo in great detail. Here I will just briefly give a few details about the dataset. Oh, and naturally, we are always looking for people to collaborate with so please reach-out if you're interested: bmonk @ ucsd . edu



## ADSP Dataset Compilation

ADSP has combined sequencing data from 25 different consortium studies, and makes for a fun genomics quest. On one hand, combinding data from 25 different sources provides the largest assembled Alzheimer's WES dataset to-date. On the other, this compilation is flush with stratification issues. Here is a summary table of those 25 consortium studies, and their share of participant contributions...


|	COHID	|	CONSORTIUM	|	STUDY	|	FAM	|	CASECTRL	|	ENRICH	|	CASE	|	CTRL	|	TOTAL	|
|	----	|	----	|	----	|	----	|	----	|	----	|	----	|	----	|	----	|
|	1	|	DGC	|	ACT	|	0	|	1277	|	0	|	323	|	945	|	1268	|
|	2	|	ADGC	|	ADC	|	0	|	3263	|	0	|	2438	|	817	|	3255	|
|	3	|	CHARGE	|	ARIC	|	0	|	57	|	0	|	39	|	18	|	57	|
|	4	|	CHARGE	|	ASPS	|	0	|	126	|	0	|	121	|	5	|	126	|
|	5	|	ADGC	|	CHAP	|	0	|	231	|	0	|	27	|	204	|	231	|
|	6	|	CHARGE	|	CHS	|	0	|	840	|	0	|	250	|	583	|	833	|
|	7	|	ADGC	|	CUHS	|	361	|	0	|	331	|	160	|	171	|	331	|
|	8	|	CHARGE	|	ERF	|	30	|	45	|	0	|	45	|	0	|	45	|
|	9	|	CHARGE	|	FHS	|	0	|	581	|	0	|	157	|	424	|	581	|
|	10	|	ADGC	|	GDF	|	0	|	207	|	0	|	111	|	96	|	207	|
|	11	|	ADGC	|	LOAD	|	80	|	113	|	363	|	367	|	109	|	476	|
|	12	|	ADGC	|	MAP	|	0	|	415	|	0	|	138	|	277	|	415	|
|	13	|	ADGC	|	MAYO	|	0	|	349	|	0	|	250	|	99	|	349	|
|	14	|	ADGC	|	MIA	|	61	|	181	|	19	|	186	|	14	|	200	|
|	15	|	ADGC	|	MIR	|	0	|	284	|	47	|	316	|	15	|	331	|
|	16	|	ADGC	|	MPD	|	0	|	20	|	0	|	0	|	20	|	20	|
|	17	|	ADGC	|	NCRD	|	18	|	108	|	52	|	160	|	0	|	160	|
|	18	|	ADGC	|	RAS	|	0	|	0	|	46	|	46	|	0	|	46	|
|	19	|	ADGC	|	ROS	|	0	|	351	|	0	|	154	|	197	|	351	|
|	20	|	CHARGE	|	RS	|	0	|	1089	|	0	|	276	|	813	|	1089	|
|	21	|	ADGC	|	TARC	|	0	|	144	|	0	|	132	|	12	|	144	|
|	22	|	ADGC	|	TOR	|	0	|	0	|	9	|	9	|	0	|	9	|
|	23	|	ADGC	|	VAN	|	6	|	235	|	1	|	210	|	26	|	236	|
|	24	|	ADGC	|	WCAP	|	0	|	147	|	3	|	34	|	116	|	150	|
|	25	|	ADGC 	|	UPN	|	40	|	0	|	3	|	0	|	0	|	0	|
