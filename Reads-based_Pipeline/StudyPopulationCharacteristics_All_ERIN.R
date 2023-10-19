########## R Code for Study Characteristics Analysis - All ERIN Metagenome Samples ##########

# Author: Zoe Hansen
# Last Modified: 2020.12.29

# This script is meant to perform summaries of study participant demographics using the 
# 'ERIN_full_metadata" spreadsheet

# The metadata will be used to identify various factors associated with Case, Controls, etc. for a 
# brief summary in the Results section of Publication #2. 

#########################################################

library(tidyverse)

ERIN_full_meta <- read.csv('D://Manning_ERIN/ERIN_MetaData/ERIN_full_metadata.csv', header = TRUE)

ERIN_full_meta <- ERIN_full_meta %>%
  filter(!grepl('\\<23\\>', ID))%>% # We need to use this text to ensure that R is removing *exact* matches, not just the values which could match in other numbers 
  filter(!grepl('\\<66\\>', ID)) %>%
  filter(!grepl('\\<85\\>', ID))%>%
  filter(!grepl('\\<86\\>', ID)) %>%
  filter(Case.status != 'Missing')

# Drop the unused levels (Follow-up and Missing)
ERIN_full_meta$Case.status <- factor(ERIN_full_meta$Case.status)
ERIN_full_meta$Age.group <- factor(ERIN_full_meta$Age.group)

Status_age <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Age.group)
prop.table(Status_age)

Status_gender <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Gender)
prop.table(Status_gender)

Gender_age <- table(ERIN_full_meta$Gender, ERIN_full_meta$Age.group)

Status_race <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Race)
prop.table(Status_race)
Status_hispanic <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Hispanic.Ethnicity)
Status_arabic <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Arab.Ethnicity)

Status_ResidenceType <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Residence.type)
prop.table(Status_ResidenceType)

Status_CareType <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Care.type)

Status_AbdPain <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Abdominal.Pain)
Status_Diarrhea <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Diarrhea)
Status_Vomit <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Vomiting)
Status_Nausea <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Nausea)
Status_Blood <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Bloody.stool)
Status_Fever <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Fever)
Status_Stool <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Stool.type)
Status_Hospital <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Hospitalized.)


Status_Abx <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Antibiotics)
Status_AbxType <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Antibiotic.type)

Status_travel_all <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Travel.past.month)
Status_travel_US <- table(ERIN_full_meta$Case.status, ERIN_full_meta$travel.in.US)
Status_travel_int <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Travel.out.US)

Status_animal <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Animal.contact.past.week)

Status_water <- table(ERIN_full_meta$Case.status, ERIN_full_meta$Water_home)

###### We can also break down our spreadsheet by Case Status #######

Case_meta <- Campy_full_meta %>%
  filter(Case.status == 'Case')
Control_meta <- Campy_full_meta %>%
  filter(Case.status == 'Control')

# Breakdown of age groups
Case_age <- table(Case_meta$Age.group)
Control_age <- table(Control_meta$Age.group)

# Breakdown of genders
table(Case_meta$Gender)
table(Control_meta$Gender)

# Combination of Age group and Genders
table(Case_meta$Age.group, Case_meta$Gender)
table(Control_meta$Age.group, Control_meta$Gender)

#Counties
table(Case_meta$County)
table(Control_meta$County)

# Care Type, Hospitalization and Length of Hospital Stays
Case_CareType <- table(Case_meta$Care.type)
Case_Hospital <- table(Case_meta$Hospitalized.)
Case_LengthStay <- table(Case_meta$Hospital.days)

# Antibiotic Treatment and Antibiotic Type
Case_Abx <- table(Case_meta$Antibiotics)
Case_AbxType <- table(Case_meta$Antibiotic.type)
