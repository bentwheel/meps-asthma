# This R Script is only intended to be called from the 'README.Rmd' Preprocessing Chunk and not intended to run on a standalone basis.

# This chunk consists of all package installs and other environment setup considerations that should be run only once. The use of if statements below ensures that packages which are already installed are not re-installed - however, if you have a very old version of these packages, you might need to reinstall with a newer version.

# Invoke required libraries - install them in the User package library if not already installed.
if(!require(devtools)) { install.packages("devtools"); library(devtools) }
if(!require(survey)) { install.packages("survey"); library(survey) }
if(!require(foreign)) { install.packages("foreign"); library(foreign) }
if(!require(haven)) { install.packages("haven"); library(haven) }
if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse) }
if(!require(scales)) { install.packages("scales"); library(scales) }
if(!require(srvyr)) { install.packages("srvyr"); library(srvyr) }

# Use devtools to download custom MEPS dataset interface package from Github if it is not already installed.
# Fore more information on the MEPS package and how to use it, refer to the following URL:
# https://github.com/HHS-AHRQ/MEPS/tree/master/R#all-data-years-using-the-meps-package
if(!require(MEPS)) { devtools::install_github("e-mitchell/meps_r_pkg/MEPS"); library(MEPS) }


# Set options to deal with lonely PSUs. A PSU is a Primary Sampling Unit. Primary Sampling Units are divided among several sampling strata. The MEPS survey design package provides appropriate sampling weights for each strata in order for each sampling unit to be re-weighted in proportion to the POI (Population of Interest - in this case, the entire US).

# In some cases, our analysis might necessitate drilling down to very small subpopulations, which will whittle down the membership of some strata to 1 or fewer members. In this case, we might encounter some strata with a single member - a "lonely PSU" - at which point the Survey package will error out when computing the sample variance to determine the standard error for a particular statistic (mean, total sum, etc.). Setting this option will adjust the data in the single-PSU stratum so that it is centered at the entire sample mean instead of the particular stratum mean, which tends to be a more conservative computation of the variance and has the effect of contributing to a wider standard error estimate for any statistic of interest in the analysis.

# In short, the following line will conservatively re-center certain data points so that a standard error for statistics of interest (mean, total sum, etc.) is computable for all strata - even those containing a single PSU, with the tradeoff of a larger (and more conservative) magnitude of standard error.

# An excellent article that goes into more detail about this process (and expresses some concern about the magnitude of overconservatism that R's survey package employs in re-centering the lonely PSU mean) can be read here:
# https://www.practicalsignificance.com/posts/bugs-with-singleton-strata/

options(survey.lonely.psu='adjust')

# The following lines of code make use of the custom MEPS R package developed by the MEPS staff. This package is not on CRAN and must be downloaded via Github using the "devtools" package, which is done in the "Initial_Setup" chunk of code in the 'README.Rmd' markdown file.

# Due to the desire to have statistically meaningful results on display in later sections of this analysis, I have decided to pull MEPS data files spanning the full years 2018 through 2021.

# Download the 2021 Full Year Consolidated Data File
# https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-233
fyc21 = MEPS::read_MEPS(year = 2021, type = "FYC")

# Download the 2020 Full Year Consolidated Data File
# https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-224
fyc20 = MEPS::read_MEPS(year = 2020, type = "FYC")

# Download the 2019 Full Year Consolidated Data File
# https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-216
fyc19 = MEPS::read_MEPS(year = 2019, type = "FYC")

# Download the 2018 Full Year Consolidated Data File
# https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-209
fyc18 = MEPS::read_MEPS(year = 2018, type = "FYC")

# From the MEPS website at https://meps.ahrq.gov/mepsweb/data_stats/download_data_files.jsp:
# "The pooled linkage file contains the standardized variance strata and PSU variables for a pooled analysis of multiple years of MEPS data."
linkage = MEPS::read_MEPS(type = "Pooled linkage")

# Custom function creation for subsetting, renaming columns based on the year, and adding a year identifier row for datasets from different years, 2018 - 2021.
process_fyc_data <- function(fyc_data, year) {
  # Extract the last two digits of the year for regex and replacement operations
  year_suffix <- substr(as.character(year), 3, 4)
  
  # Regex pattern to find year-specific column names
  year_pattern <- as.character(year_suffix)
  
  # Function to replace two-digit year in column names
  replace_year_suffix <- function(name) {
    sub(year_pattern, "YY", name)
  }
  
  # List of columns to keep (adapted to match only the two last digits of the year)
  cols_to_keep <- c("DUID", "PID", "DUPERSID", "PANEL", "SEX", "AGELAST", 
                    paste0("REGION", year_suffix), "RACETHX", 
                    paste0("POVLEV", year_suffix),
                    paste0("RXEXP", year_suffix), paste0("TOTEXP", year_suffix), 
                    paste0("RXTOT", year_suffix), "VARSTR", "VARPSU", 
                    paste0("PERWT", year_suffix, "F"), "ASTHDX", 
                    "ASATAK31", "ASDALY31", "ASTHEP31", "ASTHAGED")
  
  # Select the specified columns and rename them to end with 'YY'
  fyc_data %>%
    select(all_of(cols_to_keep)) %>%
    rename_with(.fn = replace_year_suffix, .cols = contains(year_suffix)) %>%
    mutate(MEPS_DATA_YEAR = year) # Adding a column to identify the data year
}

# Subset and process all FYC files to only those fields containing survey analytics, demographic, cost, and asthma-related information.
fyc21_processed = process_fyc_data(fyc21, 2021)
fyc20_processed = process_fyc_data(fyc20, 2020)
fyc19_processed = process_fyc_data(fyc19, 2019)
fyc18_processed = process_fyc_data(fyc18, 2018)

# Subset the pooled linkage file to just the fields we need
linkage_sub <- linkage %>% 
  select(DUPERSID, PANEL, STRA9621, PSU9621)

# Stack these files, then calculate pooled weights as 1/4 (because 4 years of pooling) of current values, and append pooled variances from subsetted linkage files.
fyc_processed <- bind_rows(fyc18_processed, fyc19_processed, fyc20_processed, fyc21_processed) %>% 
  mutate(POOLWT = PERWTYYF / 4) %>% 
  left_join(linkage_sub, by = c("DUPERSID", "PANEL"))

# Finally, I neeed to provide some natural language descriptions for some of the coded response columns in the combined 2018 - 2021 datasets.
fyc_w_desc <- fyc_processed %>% 
  arrange(MEPS_DATA_YEAR) %>% 
  mutate(MEPS_DATA_YEAR = forcats::fct_inorder(as.character(MEPS_DATA_YEAR))) %>% 
  mutate(SEX_DSC = if_else(SEX == 1, "Male", if_else(SEX == 2, "Female", "Unknown"))) %>% 
  mutate(SEX_DSC = forcats::fct_inorder(SEX_DSC)) %>% 
  mutate(AGE_GRP_2 = if_else(AGELAST >= 65, "65 and over", "Under 64")) %>% 
  mutate(AGE_GRP_3 = case_when(AGELAST < 18 ~ "Under 18",
                               AGELAST < 65 ~ "18 - 64",
                               AGELAST >= 65 ~ "65 and over",
                               T ~ as.character(AGELAST))) %>%
  mutate(AGE_GRP_5 = case_when(AGELAST < 5 ~ "Under 5",
                               AGELAST <= 17 ~ "5 - 17",
                               AGELAST <= 44 ~ "18 - 44",
                               AGELAST <= 64 ~ "45 - 64",
                               AGELAST >= 65 ~ "65 and over",
                               T ~ as.character(AGELAST))) %>%
  mutate(AGE_GRP_9 = case_when(AGELAST < 5 ~ "Under 5",
                               AGELAST <= 17 ~ "5 - 17",
                               AGELAST <= 29 ~ "18 - 29",
                               AGELAST <= 39 ~ "30 - 39",
                               AGELAST <= 49 ~ "40 - 49",
                               AGELAST <= 59 ~ "50 - 59",
                               AGELAST <= 69 ~ "60 - 69",
                               AGELAST <= 79 ~ "70 - 79",
                               AGELAST >= 80 ~ "80 and over",
                               T ~ as.character(AGELAST))) %>%
  arrange(AGELAST) %>% 
  mutate(AGE_GRP_9 = forcats::fct_inorder(AGE_GRP_9),
         AGE_GRP_5 = forcats::fct_inorder(AGE_GRP_5),
         AGE_GRP_3 = forcats::fct_inorder(AGE_GRP_3),
         AGE_GRP_2 = forcats::fct_inorder(AGE_GRP_2)) %>% 
  mutate(REGION_DSC = case_when(REGIONYY == -1 ~ "Inapplicable",
                             REGIONYY == 1 ~ "Northeast",
                             REGIONYY == 2 ~ "Midwest",
                             REGIONYY == 3 ~ "South",
                             REGIONYY == 4 ~ "West",
                             T ~ as.character(REGIONYY))) %>% 
  arrange(REGIONYY) %>% 
  mutate(REGION_DSC = forcats::fct_inorder(REGION_DSC)) %>% 
  mutate(RACETHX_DSC = case_when(RACETHX == 1 ~ "Hispanic",
                                 RACETHX == 2 ~ "Non-Hispanic White Only",
                                 RACETHX == 3 ~ "Non-Hispanic Black Only",
                                 RACETHX == 4 ~ "Non-Hispanic Asian Only",
                                 RACETHX == 5 ~ "Non-Hispanic Other Race or Multiple Race",
                                 T ~ as.character(RACETHX))) %>% 
  arrange(RACETHX) %>% 
  mutate(RACETHX_DSC = forcats::fct_inorder(RACETHX_DSC)) %>% 
  mutate(POVLEV_DSC = case_when(POVLEVYY < 100 ~ "Less than 100%",
                              POVLEVYY >= 100 & POVLEVYY < 138 ~ "100% to less than 138%",
                              POVLEVYY >= 138 & POVLEVYY < 150 ~ "138% to less than 150%",
                              POVLEVYY >= 150 & POVLEVYY < 400 ~ "150% to less than 400%",
                              POVLEVYY >= 400 ~ "400% or more",
                              T ~ as.character(POVLEVYY))) %>% 
  arrange(POVLEVYY) %>% 
  mutate(POVLEV_DSC = forcats::fct_inorder(POVLEV_DSC, ordered=T))  %>% 
  mutate(ASTHDX_DSC = case_when(ASTHDX == 1 ~ "Diagnosed with asthma", 
                              ASTHDX == 2 ~ "Not diagnosed with asthma", 
                              T ~ "Unknown or Inapplicable")) %>% 
  mutate(ASTHDX_DSC = forcats::fct_inorder(ASTHDX_DSC, ordered=T)) 

# Now we create the survey object and use this to build basic survey data tables
pooled_svydsgn <- fyc_w_desc %>% 
  as_survey_design(
    id = PSU9621,
    strata = STRA9621,
    weights = POOLWT,
    nest = T
  )

# Using the srvyr package, we can interact with the survey design object using dyplr-like verbs for subsetting and performing estimates of key statistics (means, quantiles, proportions, etc.) on by-groups.

# Asthma diagnoses by age group
asthma_by_age.data <- pooled_svydsgn %>% 
  filter(AGELAST > 0,
         POOLWT > 0) %>%  # To avoid warning, does not change results
  group_by(AGE_GRP_5, SEX_DSC, ASTHDX_DSC) %>% 
  summarize(asthma_diag_prop = survey_prop(vartype=c("se", "ci"))) %>% 
  as_tibble() %>%  #convert back to tibble so we can use ordinary dplyr verbs and not srvyr verbs
  mutate(rse = asthma_diag_prop_se / asthma_diag_prop) %>% 
  filter(ASTHDX_DSC != "Unknown or Inapplicable")


