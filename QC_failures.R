# eGTEx samples that failed QC. Data provided by Lindsay Rizzardi on Slack
# Peter Hickey
# 2018-09-26

library(here)

# Donors whose WGBS and WGS genotypes don't match:
# https://jhu-genomics.slack.com/archives/GB1UH5DB2/p1537195797000100
genotype_failures <- c(
  "GTEX-13G51-0011-R6b-SM-AFUJX",
  "GTEX-13G51-0011-R10b-SM-AFUJU",
  "GTEX-13N1W-0826-SM-ATE4J",
  "GTEX-13O3Q-0011-R3a-SM-AFUKX",
  "GTEX-T2IS-0011-R10a-SM-AFUOS",
  "GTEX-13N1W-0011-R5b-SM-AFUQW",
  "GTEX-13N1W-0011-R1b-SM-AFUQX",
  "GTEX-13N1W-0011-R7b-SM-AFUQT",
  "GTEX-13O3Q-0011-R6a-SM-AFUKS",
  "GTEX-WZTO-0426-SM-ATE6F",
  "GTEX-13N1W-0726-SM-ATE66",
  "GTEX-13FHO-1026-SM-ATE65",
  "GTEX-13N1W-0011-R4a-SM-AFUQV",
  "GTEX-13NYB-0011-R1b-SM-9VIIM",
  "GTEX-11GSP-0011-R7b-SM-9VIIR",
  "GTEX-XV7Q-2526-SM-9VXYG",
  "GTEX-13FTW-2626-SM-9VXZI",
  "GTEX-13OW8-1426-SM-9VXZL",
  "GTEX-U3ZN-2226-SM-9VXZT",
  "GTEX-1211K-0826-SM-9VXYJ",
  "GTEX-XV7Q-0326-SM-9VXYI",
  "GTEX-1122O-2426-SM-9VXZJ",
  "GTEX-1399R-0726-SM-9VXZ3",
  "GTEX-U3ZN-0326-SM-9VXYQ",
  "GTEX-YEC4-0526-SM-9VXYP",
  "GTEX-U3ZN-2226-SM-9VXZT",
  "GTEX-1211K-2026-SM-9VXZN",
  "GTEX-13FTW-0626-SM-9VXYR",
  "GTEX-13OW8-1426-SM-9VXZL",
  "GTEX-11DXX-0326-SM-9VXZM")

# Outliers based on PCA of mCG:
# https://jhu-genomics.slack.com/archives/GB1UH5DB2/p1537196090000100
mCG_PCA_outliers <- c(
  "SM-AFUJL",
  "SM-AFUKH",
  "SM-AFUJM",
  "SM-AFUJG",
  "SM-AFUN2")

save(
  genotype_failures,
  mCG_PCA_outliers,
  file = here("miscellaneous_objects", "QC_failures.rda"))
