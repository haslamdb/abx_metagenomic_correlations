# Clinical Data Request: BSI Etiology Analysis

## Study Overview
Analysis of bloodstream infection (BSI) organism etiology in relation to prior antibiotic exposure, with the hypothesis that time since last antibiotic affects which organisms cause BSI (differential recovery kinetics).

**Time Frame:** 2013 - present (13 years)
**Population:** All patients with positive blood cultures

---

## 1. PATIENT DEMOGRAPHICS

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| MRN | Medical record number (de-identified OK) | String | Required |
| DOB or Age | Date of birth or age at encounter | Date/Integer | Required |
| Gender | Male/Female/Other | Categorical | Required |
| Race/Ethnicity | Self-reported race and ethnicity | Categorical | Optional |
| Height | Height in cm or inches | Numeric | Optional |
| Weight | Weight in kg or lbs (closest to admission) | Numeric | Optional |

---

## 2. ENCOUNTER/ADMISSION DATA

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Encounter_ID | Unique identifier for hospitalization | String | Required |
| Admission_Date | Date of hospital admission | Date | Required |
| Admission_Time | Time of admission (if available) | Time | Optional |
| Discharge_Date | Date of hospital discharge | Date | Required |
| Discharge_Disposition | Home, SNF, deceased, hospice, etc. | Categorical | High |
| Admission_Source | ED, transfer, direct admit, OR | Categorical | Medium |
| Admission_Service | Medicine, Surgery, Oncology, etc. | Categorical | High |
| ICU_Admission_Date | Date of first ICU admission (if any) | Date | High |
| ICU_Discharge_Date | Date of ICU discharge | Date | High |
| ICU_LOS | Total ICU days during encounter | Integer | High |
| Hospital_LOS | Total hospital length of stay | Integer | Required |
| Readmission_30d | Readmission within 30 days (Y/N) | Boolean | Medium |

---

## 3. BLOOD CULTURE DATA (you may already have this)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Culture_ID | Unique identifier for culture | String | Required |
| Collection_Date | Date blood culture collected | Date | Required |
| Collection_Time | Time of collection | Time | High |
| Organism | Species identified | String | Required |
| Organism_Code | Standardized organism code if available | String | Optional |
| Polymicrobial | Multiple organisms in same culture (Y/N) | Boolean | High |
| Bottle_Type | Aerobic/Anaerobic | Categorical | Medium |
| Bottles_Positive | Number of bottles positive | Integer | High |
| Time_to_Positivity | Hours to positive signal | Numeric | Medium |
| Contaminant_Flag | Suspected contaminant (Y/N) | Boolean | Medium |
| Source_Documented | Documented source of BSI (line, urine, etc.) | Categorical | High |

---

## 4. SUSCEPTIBILITY DATA (you may already have this)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Culture_ID | Links to blood culture | String | Required |
| Antibiotic_Tested | Name of antibiotic | String | Required |
| MIC | Minimum inhibitory concentration | Numeric | High |
| Interpretation | S/I/R | Categorical | Required |
| ESBL_Screen | ESBL positive (Y/N) | Boolean | High |
| CRE_Flag | Carbapenem-resistant (Y/N) | Boolean | High |
| VRE_Flag | Vancomycin-resistant Enterococcus | Boolean | High |
| MRSA_Flag | Methicillin-resistant S. aureus | Boolean | High |

---

## 5. ANTIBIOTIC EXPOSURE (you have this)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| MRN | Patient identifier | String | Required |
| Drug_Name | Antibiotic name | String | Required |
| Drug_Class | Antibiotic class (if coded) | String | Medium |
| Route | IV, PO, IM, topical, etc. | Categorical | Required |
| Dose | Dose amount | Numeric | Medium |
| Dose_Unit | mg, g, etc. | String | Medium |
| Frequency | Daily, BID, TID, etc. | Categorical | Medium |
| Start_Date | Date antibiotic started | Date | Required |
| Stop_Date | Date antibiotic stopped | Date | Required |
| Duration_Days | Total days on antibiotic | Integer | High |

---

## 6. OTHER MEDICATIONS

### 6a. High Priority Medications

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| **PPI_Exposure** | Proton pump inhibitor use | Boolean + Dates | High |
| **H2_Blocker** | H2 receptor antagonist use | Boolean + Dates | High |
| **Immunosuppressants** | Tacrolimus, cyclosporine, MMF, sirolimus | Boolean + Dates | High |
| **Systemic_Steroids** | Prednisone, dexamethasone, hydrocortisone (>10mg/day) | Boolean + Dates | High |
| **Chemotherapy** | Any antineoplastic agent | Boolean + Dates | High |
| **G-CSF** | Filgrastim, pegfilgrastim | Boolean + Dates | High |
| **TPN** | Total parenteral nutrition | Boolean + Dates | High |
| **Tube_Feeds** | Enteral nutrition | Boolean + Dates | Medium |

### 6b. Other Relevant Medications

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Probiotics | Lactobacillus, Saccharomyces, etc. | Boolean + Dates | Medium |
| Laxatives | Lactulose, PEG, senna | Boolean + Dates | Low |
| Antidiarrheals | Loperamide, etc. | Boolean + Dates | Low |
| Opioids | Morphine, fentanyl, etc. (affect gut motility) | Boolean + Dates | Low |
| Antacids | Calcium carbonate, etc. | Boolean + Dates | Low |

---

## 7. LABORATORY VALUES

### 7a. Hematology (Critical for neutropenia status)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| WBC | White blood cell count | Numeric (K/uL) | Required |
| ANC | Absolute neutrophil count | Numeric (K/uL) | Required |
| ALC | Absolute lymphocyte count | Numeric (K/uL) | High |
| Hemoglobin | Hemoglobin level | Numeric (g/dL) | Medium |
| Platelets | Platelet count | Numeric (K/uL) | Medium |
| Lab_Date | Date of lab draw | Date | Required |

*Request: All values within 7 days before and 3 days after BSI*

### 7b. Chemistry/Metabolic

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Creatinine | Serum creatinine | Numeric (mg/dL) | High |
| BUN | Blood urea nitrogen | Numeric (mg/dL) | Medium |
| eGFR | Estimated GFR | Numeric (mL/min) | High |
| Total_Bilirubin | Total bilirubin | Numeric (mg/dL) | High |
| Direct_Bilirubin | Direct bilirubin | Numeric (mg/dL) | Medium |
| ALT | Alanine aminotransferase | Numeric (U/L) | Medium |
| AST | Aspartate aminotransferase | Numeric (U/L) | Medium |
| Albumin | Serum albumin | Numeric (g/dL) | High |
| Glucose | Blood glucose | Numeric (mg/dL) | Medium |
| Lactate | Serum lactate | Numeric (mmol/L) | High |

### 7c. Inflammatory Markers

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| CRP | C-reactive protein | Numeric (mg/L) | Medium |
| Procalcitonin | Procalcitonin level | Numeric (ng/mL) | Medium |
| ESR | Erythrocyte sedimentation rate | Numeric (mm/hr) | Low |

---

## 8. PROCEDURES AND DEVICES

### 8a. Vascular Access (Critical for line-associated BSI)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Central_Line | Any central venous catheter (Y/N) | Boolean | High |
| CVC_Type | PICC, tunneled, non-tunneled, port | Categorical | High |
| CVC_Insert_Date | Date of insertion | Date | High |
| CVC_Remove_Date | Date of removal | Date | High |
| CVC_Site | Subclavian, IJ, femoral, arm (PICC) | Categorical | Medium |
| Arterial_Line | Arterial catheter (Y/N) | Boolean | Medium |
| HD_Catheter | Hemodialysis catheter (Y/N) | Boolean | High |

### 8b. Other Devices

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Foley_Catheter | Urinary catheter (Y/N) | Boolean | High |
| Foley_Days | Duration of catheterization | Integer | High |
| NG_Tube | Nasogastric tube (Y/N) | Boolean | Medium |
| Mechanical_Vent | On mechanical ventilation (Y/N) | Boolean | High |
| Vent_Days | Duration of ventilation | Integer | High |
| Tracheostomy | Tracheostomy present (Y/N) | Boolean | Medium |

### 8c. Procedures (within 30 days prior to BSI)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| GI_Surgery | Abdominal/GI surgery | Boolean + Date | High |
| GI_Procedure | Endoscopy, colonoscopy, ERCP | Boolean + Date | High |
| Cardiac_Surgery | CABG, valve replacement, etc. | Boolean + Date | Medium |
| Other_Surgery | Any other surgical procedure | Boolean + Date | Medium |
| Procedure_Code | CPT or ICD procedure codes | String | Medium |
| Bone_Marrow_Biopsy | BMB performed | Boolean + Date | Medium |
| Lumbar_Puncture | LP performed | Boolean + Date | Low |

---

## 9. DIAGNOSES AND COMORBIDITIES

### 9a. Admission Diagnosis

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Primary_Dx_Code | Primary diagnosis ICD-10 code | String | Required |
| Primary_Dx_Desc | Primary diagnosis description | String | Required |
| All_Dx_Codes | All diagnosis codes for encounter | String (list) | High |
| DRG | Diagnosis-related group | String | Medium |

### 9b. Comorbidities (Charlson/Elixhauser components)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Diabetes | Diabetes mellitus (Y/N) | Boolean | High |
| Diabetes_Complicated | With end-organ damage | Boolean | High |
| CKD | Chronic kidney disease | Boolean | High |
| CKD_Stage | Stage 1-5 or dialysis | Categorical | High |
| ESRD_Dialysis | On dialysis (Y/N) | Boolean | High |
| Cirrhosis | Liver cirrhosis (Y/N) | Boolean | High |
| CHF | Congestive heart failure | Boolean | Medium |
| COPD | Chronic lung disease | Boolean | Medium |
| Malignancy | Active malignancy | Boolean | High |
| Malignancy_Type | Solid tumor, heme malignancy, etc. | Categorical | High |
| Metastatic | Metastatic cancer | Boolean | High |
| HIV | HIV infection | Boolean | High |
| HIV_AIDS | AIDS-defining illness | Boolean | High |
| Transplant_Solid | Solid organ transplant | Boolean | High |
| Transplant_HSCT | Hematopoietic stem cell transplant | Boolean | High |
| Transplant_Type | Specific organ/type | Categorical | High |
| Transplant_Date | Date of transplant | Date | High |
| Autoimmune | Autoimmune disease requiring immunosuppression | Boolean | Medium |
| IBD | Inflammatory bowel disease | Boolean | High |
| GI_Malignancy | GI tract cancer | Boolean | High |

### 9c. Calculated Scores (if available)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Charlson_Score | Charlson Comorbidity Index | Integer | Medium |
| APACHE_II | APACHE II score (ICU patients) | Integer | Medium |
| SOFA_Score | SOFA score at BSI | Integer | Medium |
| Pitt_Bacteremia | Pitt Bacteremia Score | Integer | High |
| NEWS | National Early Warning Score | Integer | Medium |

---

## 10. CLINICAL STATUS AT TIME OF BSI

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Location_at_BSI | Ward, ICU, ED, OR | Categorical | High |
| Temperature | Temperature at time of culture | Numeric (°C or °F) | High |
| Heart_Rate | Heart rate | Numeric (bpm) | Medium |
| Blood_Pressure_Sys | Systolic BP | Numeric (mmHg) | High |
| Blood_Pressure_Dia | Diastolic BP | Numeric (mmHg) | Medium |
| Respiratory_Rate | Respiratory rate | Numeric (breaths/min) | Medium |
| O2_Saturation | Oxygen saturation | Numeric (%) | Medium |
| FiO2 | Fraction of inspired oxygen | Numeric (%) | Medium |
| Vasopressors | On vasopressor support (Y/N) | Boolean | High |
| Sepsis_Diagnosis | Sepsis/severe sepsis/septic shock | Categorical | High |

---

## 11. OUTCOMES

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Mortality_Inpatient | Died during hospitalization | Boolean | Required |
| Mortality_30d | Died within 30 days of BSI | Boolean | Required |
| Mortality_90d | Died within 90 days of BSI | Boolean | High |
| Death_Date | Date of death (if applicable) | Date | High |
| ICU_After_BSI | Required ICU admission after BSI | Boolean | High |
| Recurrent_BSI | Another BSI within 30 days | Boolean | High |
| LOS_After_BSI | Hospital days after BSI | Integer | High |
| Appropriate_Abx_24h | Appropriate antibiotic within 24h | Boolean | High |
| Time_to_Appropriate | Hours to appropriate antibiotic | Numeric | High |

---

## 12. PRIOR MICROBIOLOGY (within 90 days before BSI)

| Variable | Description | Format | Priority |
|----------|-------------|--------|----------|
| Prior_BSI | Previous BSI in past 90 days | Boolean | High |
| Prior_BSI_Organism | Organism from prior BSI | String | High |
| Prior_Positive_Cx | Any prior positive culture | Boolean | High |
| Prior_Cx_Site | Site of prior culture (urine, sputum, wound) | Categorical | High |
| Prior_Cx_Organism | Organism from prior culture | String | High |
| Colonization_MRSA | Known MRSA colonization | Boolean | High |
| Colonization_VRE | Known VRE colonization | Boolean | High |
| Colonization_CRE | Known CRE colonization | Boolean | High |
| Colonization_ESBL | Known ESBL colonization | Boolean | High |
| C_diff_History | Prior C. difficile infection | Boolean | High |

---

## Data Format Recommendations

1. **One row per BSI episode** for main analysis table
2. **Separate tables** for:
   - Time-varying labs (long format: MRN, Date, Lab_Name, Value)
   - Medications (long format: MRN, Drug, Start_Date, Stop_Date)
   - Procedures (long format: MRN, Procedure, Date)
   - All diagnosis codes (long format: MRN, Encounter, ICD_Code)

3. **Date handling**: All dates in consistent format (YYYY-MM-DD preferred)

4. **Missing data**: Please indicate missing vs. not measured vs. negative

---

## Minimum Viable Dataset

If extraction is challenging, the **absolute minimum** for meaningful analysis:

1. Patient ID (MRN)
2. BSI culture date and organism
3. Hospital admission and discharge dates
4. Antibiotic exposures with dates (already have)
5. Age, gender
6. ICU status at BSI
7. WBC/ANC closest to BSI
8. 30-day mortality
9. Central line present (Y/N)
10. Primary diagnosis code

---

## Contact

For questions about this data request, please contact: [Your name/email]

Analysis will be conducted in accordance with IRB protocol: [Protocol number]
