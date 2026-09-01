# Anticholinergic_ISACS
This repository contains the data files as well as code to reproduce the results presented in the article: "Systems pharmacology identifies centrally acting anticholinergic activity of drugs and their aging-related signatures".

Below is the list of R scripts:\
0_HelperFunctions.R - helper functions for network medicine approach.\
1_ConstructingPPI.R - harmonizing PPI networks from Gross et al., 2026 (https://github.com/BnayaGross/Longevity-module/) and Zhou et al., 2020 (https://github.com/ChengF-Lab/2019-nCoV).\
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
