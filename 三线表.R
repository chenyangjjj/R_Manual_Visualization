graphics.off()
rm(list=ls())
library(stringr)
library(dplyr)
library(rprojroot)
library(ggseg)
library(ggpubr)
library(table1)

########################################################################################################## 
## Cutoffs of AD biomarkers based on ROC negative CN Vs. positive AD + MCI
########################################################################################################## 
Cutoff_CSF_ABETA42 = 978
Cutoff_CSF_PTAU = 23
Cutoff_CSF_TAU = 234
Cutoff_CSF_ABETA42.ABETA40 = 0.054  # 0.051 for Gaussian-mixture model; 0.054 for ROC analysis youden index
Cutoff_CSF_PTAU.ABETA40 = 0.0012
Cutoff_aHCV = 6787
Cutoff_Meta_ROI_thickness = 2.6
Cutoff_FDG_Ref_PON = 1.21
Cutoff_FTP_Inf_CER_GM_PVC_Meta_TEM_ROI_SUVR = 1.55
Cutoff_FTP_Inf_CER_GM_Non_PVC_Meta_TEM_ROI_SUVR = 1.25
Cutoff_FTP_Inf_CER_GM_PVC_ENTORHINAL_SUVR = 1.90
Cutoff_FTP_Inf_CER_GM_Non_PVC_ENTORHINAL_SUVR = 1.21
Cutoff_FTP_WM_PVC_Meta_TEM_ROI_SUVR = 2.05
Cutoff_FTP_WM_Non_PVC_Meta_TEM_ROI_SUVR = 1.07
Cutoff_FTP_WM_PVC_ENTORHINAL_SUVR = 2.23
Cutoff_FTP_WM_Non_PVC_ENTORHINAL_SUVR = 1.06
## Still cutoffs form LONI and GMM method for ABETA42/40 ratio
Cutoff_CSF_MASS_ABETA42 = 1079
Cutoff_CSF_MASS_ABETA42.ABETA40 = 0.138
Cutoff_CSF_MASS_ABETA42.ABETA38 = 0.627
Cutoff_SUVR_CER = 1.11
Cutoff_SUVR_Big_Ref = 0.82
Cutoff_SUVR_FBB_CER = 1.08
Cutoff_Centiloid = 24.4  # this is the new Aβ PET cut-off


Data_ADNI_CSF_PET_Abeta_Tau = 
  read.csv("D:/Project/guo_20210206/CSF_PET_tau/Data/Data_ADNI_AV1451_CSF_PET_Abeta_Tau_02_04_21.csv", 
           header=TRUE, sep=",",stringsAsFactors=FALSE)
Data_ADNI_CSF_PET_Abeta_Tau = 
  Data_ADNI_CSF_PET_Abeta_Tau[,2:ncol(Data_ADNI_CSF_PET_Abeta_Tau)]
Data_ADNI_CSF_PET_Abeta_Tau$Diag_AV1451_v1 = 
  factor(Data_ADNI_CSF_PET_Abeta_Tau$Diag_AV1451_v1 ,levels = c("CU","MCI","AD"))

########################################################################################################## 
# set new cut-off
Data_ADNI_CSF_PET_Abeta_Tau$Positivity_CSF_ABETA42_ABETA40_Closest_AV1451_v1_New =
  ifelse(Data_ADNI_CSF_PET_Abeta_Tau$CSF_ABETA42.ABETA40_Closest_AV1451_v1 < 
           Cutoff_CSF_ABETA42.ABETA40,"CSF_ABETA42_P","CSF_ABETA42_N")
Data_ADNI_CSF_PET_Abeta_Tau$Positivity_PET_ABETA_Closest_AV1451_v1_New =
  ifelse(Data_ADNI_CSF_PET_Abeta_Tau$PET_ABETA_CL_Closest_AV1451_v1 > 
           Cutoff_Centiloid,"COMPOSITE_P","COMPOSITE_N")

########################################################################################################## 
# define stages  CSF-/PET-  CSF+/PET-  CSF-/PET+  CSF+/PET+
Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta=NA
Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta[
  (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_PET_ABETA_Closest_AV1451_v1_New=="COMPOSITE_N")&
    (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_CSF_ABETA42_ABETA40_Closest_AV1451_v1_New==
       "CSF_ABETA42_N")] = "CSF-/PET-"
Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta[
  (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_PET_ABETA_Closest_AV1451_v1_New=="COMPOSITE_N")&
    (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_CSF_ABETA42_ABETA40_Closest_AV1451_v1_New==
       "CSF_ABETA42_P")] = "CSF+/PET-"
Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta[
  (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_PET_ABETA_Closest_AV1451_v1_New=="COMPOSITE_P")&
    (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_CSF_ABETA42_ABETA40_Closest_AV1451_v1_New==
       "CSF_ABETA42_N")] = "CSF-/PET+"
Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta[
  (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_PET_ABETA_Closest_AV1451_v1_New=="COMPOSITE_P")&
    (Data_ADNI_CSF_PET_Abeta_Tau$Positivity_CSF_ABETA42_ABETA40_Closest_AV1451_v1_New==
       "CSF_ABETA42_P")] = "CSF+/PET+"

Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta=
  factor(Data_ADNI_CSF_PET_Abeta_Tau$Stage_CSF_PET_Abeta,
         levels = c("CSF-/PET-","CSF+/PET-","CSF-/PET+","CSF+/PET+"))


Data_ADNI_CSF_PET_Abeta_Tau_v1 = subset(Data_ADNI_CSF_PET_Abeta_Tau, AV1451_Scan_ID==1)
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1 = subset(Data_ADNI_CSF_PET_Abeta_Tau_v1, (Diag_AV1451_v1 != "AD"))

table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1$Stage_CSF_PET_Abeta)
nrow(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1)


#table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1$Stage_CSF_PET_Abeta)
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 = 
  subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_v1, 
         !((Stage_CSF_PET_Abeta == "CSF-/PET-")&
             ((Positivity_CSF_PTAU_ABETA40_Closest_AV1451_v1 == "CSF_PTAU_ABETA40_P")|
                (Positivity_AV1451_ENTORHINAL_Non_PVC_Inf_CER_GM_v1 == "FTP_P")|
                (Positivity_AV1451_Meta_TEM_ROI_Non_PVC_Inf_CER_GM_v1 == "FTP_P"))))

########################################################################################################## 
# remove PACC and aHCV NA value
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 = 
  subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1, 
         !(((aHCV_Closest_AV1451_v1 == "NA")|
              (PACC_Digit_LONI_Closest_AV1451_v1 == "NA")|
              (Meta_ROI_thickness_Closest_AV1451_v1 == "NA"))))

########################################################################################################## 


nrow(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A =
  subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1,
         !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 1.6)&(Stage_CSF_PET_Abeta == "CSF+/PET-")))
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1 =
  subset(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1_A,
         !((ENTORHINAL_SUVR_FTP_Inf_CER_GM_Non_PVC_v1 > 2)&(Stage_CSF_PET_Abeta == "CSF+/PET+")))
table(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Stage_CSF_PET_Abeta)
nrow(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)



Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Gender=
  factor(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Gender,
         levels = c("Female","Male"))
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$APOE_status=
  factor(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$APOE_status,
         levels = c("APOE4_carrier","APOE4_Non_carrier"))
Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Diag_AV1451_v1=
  factor(Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1$Diag_AV1451_v1,
         levels = c("CU","MCI"))
table1(~ Diag_AV1451_v1+Age_AV1451_v1 + Gender + Education + 
         APOE_status+aHCV_Closest_AV1451_v1+PACC_Digit_LONI_Closest_AV1451_v1
       |Stage_CSF_PET_Abeta,
       data=Data_ADNI_CSF_PET_Abeta_Tau_CU_MCI_RM_Ref_T_P_v1)

