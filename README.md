# Anticholinergic_ISACS
This repository contains the code to reproduce the results presented in the article: "Systems pharmacology identifies centrally acting anticholinergic activity of drugs and their aging-related signatures".

Access to raw data from databases:
DrugBank (under the academic license agreement): https://go.drugbank.com/releases/latest \
DGIdb: https://dgidb.org/downloads \
TTD: https://ttd.idrblab.cn/full-data-download \
BindingDB: https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp \
ChEMBL: https://www.ebi.ac.uk/chembl/ \
Drug Repurposing Hub (Broad Institute): https://repo-hub.broadinstitute.org/repurposing \
PSICHIC: https://github.com/huankoh/PSICHIC/tree/main/examples \
HPA: https://www.proteinatlas.org/ \
ADMETLAB3.0: https://admetlab3.scbdd.com/ \
ONSIDES: https://onsidesdb.org/ \
SIDER: https://sideeffects.embl.de/ \
DrugMatrix Affymetrix Superseries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57822 \
DrugMatrix CodeLink Superseries: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59927 \

Access to raw data from publications:
Cmax (Smit et al., 2020): https://doi.org/10.1021/acs.chemrestox.0c00294 \
Chew's SAA list (Chew et al., 2008): https://doi.org/10.1111/j.1532-5415.2008.01737.x \
ACSBC (Al Rihani et al., 2021): https://doi.org/10.1007/s40266-021-00895-x \
ORCA (Simal et al., 2025): https://doi.org/10.1093/ageing/afaf313 \
ABC, ADS, ACB, ARS, AAS, and ALS (Lozano-Ortega et al., 2020): https://doi.org/10.1016/j.archger.2019.05.010 \
ABS (Yamada et al., 2023): https://doi.org/10.1111/ggi.14619 \
ATS (Xu et al., 2017): https://doi.org/10.1177/2042098617725267 \
AEC (Bishara et al., 2016): https://doi.org/10.1002/gps.4507 \

Below is the list of R scripts:\
0_HelperFunctions.R - helper functions for network medicine approach.\
1_ConstructingPPI.R - harmonizing PPI networks from Gross et al., 2026 (https://github.com/BnayaGross/Longevity-module/) and Zhou et al., 2020 (https://github.com/ChengF-Lab/2019-nCoV). \
2_DrugBank_DataReconciliation.R - processing raw DrugBank data.\
2_DrugBank_DataReconciliation_ATC_OnSIDES.R - process ATC codes and ONSIDES.\
3_BindingAffinity.R - harmonizing binding affinity data across BindingDB, TTD, and CHEMBL.\
4_NetworkSetup.R - Loading all necessary files for calculating network distance and z-score.\
5_1_ModelOptimization_Loop.R - network model development and optimization.\
5_2_ModelOptimization_Viz.R - comparing different network models.\
5_3_Run7000Drugs_Anticholinergic.R - calculate network distance and z-score for anticholinergic properties.\
6_1_RedefiningAnticholinergics_Setup.R - Loading all necessary files for "RedefiningAnticholinergic" scripts.\
6_2_RedefiningAnticholinergics_StandardAnalysis.R - network analysis on anticholinergic properties.\
6_3_RedefiningAnticholinergics_Classification.R - development and creation of ISACS and drug-centric anticholinergic effect association analysis.\
6_4_RedefiningAnticholinergics_CompareScales.R - comparing ISACS with existing anticholinergic scales and characterizing safety parameters.\
7_1_Compiling_DrugMatrix.R - DrugMatrix data ascertainment from GEO.\
7_2_tAge_analysis.R - tAge calculation from DrugMatrix dataset.\
7_3_tAge_visualization.R - tAge association analysis with ISACS on DrugMatrix.\
8_1_NHANES_setup.R - NHANES data ascertainment and processing.\
8_2_NHANES_analysis_and_visualization.R - phenotypic age gap and cognitive outcome association analysis with ISACS on NHANES cohort.\
