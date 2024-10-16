sxtTools::setwd_project()
rm(list=ls())

###treat data
treatment = 
readr::read_csv("data/metadata/Treatment.csv") %>% 
  as.data.frame()

setwd("data/metabolomics/data_preparation/")

library(tidyverse)

####load data
load("../data_cleaning/POS_NEG/expression_data")
load("../data_cleaning/POS_NEG/sample_info")
load("../data_cleaning/POS_NEG/variable_info")

# ###annotation result
# load("../metabolite_annotation/POS/NCE25_rp_result")
# load("../metabolite_annotation/POS/NCE50_rp_result")
# 
# annotation_table_pos = NCE25_rp_result
# nce50_pos = NCE50_rp_result
# names(nce50_pos)
# names(annotation_table_pos)
# 
# annotation_table_pos =
#   metID::get_identification_table_all(annotation_table_pos,
#                                       nce50_pos,
#                                       candidate.num = 1)
# 
# annotation_table_pos = 
# annotation_table_pos %>% 
#   dplyr::select(name:rt, MS2.spectra.name:Level) %>% 
#   dplyr::filter(!is.na(Level)) %>% 
#   dplyr::filter(Level != 3)
# 
# save(annotation_table_pos, file = "annotation_table_pos")
# 
# 
# load("../metabolite_annotation/NEG/NCE25_rp_result")
# load("../metabolite_annotation/NEG/NCE50_rp_result")
# 
# annotation_table_neg = NCE25_rp_result
# nce50_neg = NCE50_rp_result
# names(nce50_neg)
# names(annotation_table_neg)
# 
# annotation_table_neg =
#   metID::get_identification_table_all(annotation_table_neg,
#                                       nce50_neg,
#                                       candidate.num = 1)
# 
# annotation_table_neg = 
#   annotation_table_neg %>% 
#   dplyr::select(name:rt, MS2.spectra.name:Level) %>% 
#   dplyr::filter(!is.na(Level)) %>% 
#   dplyr::filter(Level != 3)
# 
# save(annotation_table_neg, file = "annotation_table_neg")

load("annotation_table_pos")
load("annotation_table_neg")

which(annotation_table_pos$Compound.name == "Metronidazole")
which(annotation_table_neg$Compound.name == "Metronidazole")

grep("dihydroartemisinin",annotation_table_pos$Compound.name)
grep("piperaquine",annotation_table_pos$Compound.name)
grep("dihydroartemisinin",annotation_table_neg$Compound.name)
grep("piperaquine",annotation_table_neg$Compound.name)

grep("sulfadoxine",annotation_table_pos$Compound.name)
grep("pyrimethamine",annotation_table_pos$Compound.name)
grep("sulfadoxine",annotation_table_neg$Compound.name)
grep("pyrimethamine",annotation_table_neg$Compound.name)

annotation_table_pos[113,]$name

load("metlinDatabase0.0.2")
load("../metabolite_annotation/POS/NCE25_rp_result")
load("../metabolite_annotation/POS/NCE50_rp_result")

library(metID)
names(NCE25_rp_result)
plot = 
metID::ms2plot(object = NCE25_rp_result[[4]], database = metlinDatabase0.0.2,
               which.peak = "M172T193_4_POS")

ggsave(plot, filename = "Metronidazole_ms2_plot.pdf", width = 9, height = 7)

annotation_table_pos = 
annotation_table_pos %>% 
  dplyr::select(name, MS2.spectra.name:Level)

annotation_table_neg = 
  annotation_table_neg %>% 
  dplyr::select(name, MS2.spectra.name:Level)

annotation_table_pos$name

annotation_table = rbind(annotation_table_pos, annotation_table_neg) %>% 
  dplyr::filter(Level != 3)

variable_info = 
variable_info %>% 
  dplyr::left_join(annotation_table, by = c("variable_id" = "name"))

dir.create("peak")
dir.create("metabolite")


####load sample information
clinical_info_pos = readr::read_csv("clinical.worklist.pos.csv")
clinical_info_neg = readr::read_csv("clinical.worklist.neg.csv")

clinical_info_pos = 
clinical_info_pos %>% 
  dplyr::mutate(sample.name = stringr::str_replace(sample.name, "posi\\_", ""))

clinical_info_neg = 
  clinical_info_neg %>% 
  dplyr::mutate(sample.name = stringr::str_replace(sample.name, "neg\\_", ""))

intersect_sample_id =
  intersect(clinical_info_pos$sample.name,
            clinical_info_neg$sample.name)

clinical_info_pos =
  clinical_info_pos %>% 
  dplyr::filter(sample.name %in% intersect_sample_id) %>% 
  dplyr::arrange(sample.name)

clinical_info_neg =
  clinical_info_neg %>% 
  dplyr::filter(sample.name %in% intersect_sample_id) %>% 
  dplyr::arrange(sample.name)

clinical_info_pos$RandomNumber == clinical_info_neg$RandomNumber
clinical_info_pos$New.ID == clinical_info_neg$New.ID
clinical_info_pos$id.x == clinical_info_neg$id.x
clinical_info_pos$SampleDate.x == clinical_info_neg$SampleDate.x
clinical_info_pos$SGA == clinical_info_neg$SGA

clinical_info = clinical_info_pos

clinical_info = 
clinical_info %>% 
  dplyr::select(RandomNumber:sample.name, id.x:sample_GA, group) %>% 
  dplyr::select(-c(id.y, SampleDate.y)) %>% 
  dplyr::rename(id = id.x, SampleDate = SampleDate.x)

sample_info = 
sample_info %>% 
  dplyr::left_join(clinical_info, by = c("sample_id" = "sample.name"))

sample_info = 
  sample_info %>% 
  dplyr::rename(subject_id = id)

sample_info = 
sample_info %>% 
  dplyr::left_join(treatment, by = c("subject_id" = "Maternal_ID")) %>% 
  dplyr::rename(treatment = "Treatment group")


table(sample_info$matmalaria, sample_info$treatment)

save(expression_data, file = "peak/expression_data")
save(variable_info, file = "peak/variable_info")
save(sample_info, file = "peak/sample_info")

###only remain the metabolites
variable_info = 
variable_info %>% 
  dplyr::filter(!is.na(Level))

expression_data = 
  expression_data[variable_info$variable_id,]

# ##remove duplicated
# variable_info$Compound.name
# library(plyr)
# variable_info = 
# variable_info %>% 
#   plyr::dlply(.variables = .(Compound.name)) %>% 
#   purrr::map(function(x){
#     x = 
#       x %>% 
#       dplyr::filter(SS == max(SS)) %>% 
#       head()
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# expression_data = 
#   expression_data[variable_info$variable_id,]

save(expression_data, file = "metabolite/expression_data")
save(variable_info, file = "metabolite/variable_info")
save(sample_info, file = "metabolite/sample_info")
