sxtTools::setwd_project()
rm(list=ls())
source("code/tools.R")
library(tidyverse)

sxtTools::setwd_project()
setwd("data/metabolomics/ga_biomarker")

####load data
load("../data_preparation/metabolite/expression_data")
load("../data_preparation/metabolite/sample_info")
load("../data_preparation/metabolite/variable_info")

####remove sampleGA which are NA
sample_info = 
sample_info %>% 
  dplyr::filter(!is.na(sample_GA)) %>% 
  dplyr::mutate(SGA = case_when(
    SGA == 1 ~ "SGA",
    SGA == 0 ~ "AGA"
  ))

expression_data =
  expression_data[,sample_info$sample_id]

# ###remove some metabolites
# variable_info = 
#   variable_info %>% 
#   dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2", "msDatabase_rplc0.0.2",
#                                 "nistDatabase0.0.2"))

expression_data = expression_data[variable_info$variable_id,]

###find the biomakers which are different in trimester2 and trimester3

fc_p = 
  expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
purrr::map(function(x){
  x = as.numeric(x)
  ###adjusted ga
  control_idx = which(sample_info$timepoint == "2nd trimester")
  case_idx = which(sample_info$timepoint == "3rd trimester")
  x1 = x[control_idx]
  x2 = x[case_idx]
  p.value = 
  wilcox.test(x1, x2)$p.value
  fc = mean(x[case_idx])/mean(x[control_idx])
  c(fc = fc, p = p.value)
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(p.adjust = p.adjust(p, method = "BH"))

rownames(fc_p) == variable_info$variable_id

fc_p =
  data.frame(variable_info, fc_p)

sum(fc_p$p.adjust < 0.05)

#####volcano plot

plot = 
volcano_plot(
  fc = fc_p$fc,
  p_value = fc_p$p.adjust,
  p.cutoff = 0.05,
  fc.cutoff = 1,
  text = TRUE, 
  variable_id = fc_p$Compound.name, 
  point.size = 4
)
plot
ggsave(plot, filename = "volocano_plot.pdf", width = 7, height = 7)


marker =
  fc_p %>%
  dplyr::filter(p.adjust < 0.05)

temp_sample_info =
  sample_info

dim(marker)
###THDOC
grep("TH", marker$Compound.name, value = TRUE)
grep("TH", variable_info$Compound.name, value = TRUE)
####estriol-16-glucuronide
grep("Estriol", marker$Compound.name, value = TRUE)
grep("estriol", marker$Compound.name, value = TRUE)
##progesterone
grep("Progesterone", marker$Compound.name, value = TRUE)
##PE(P-16:0e/0:0)
grep("PE", marker$Compound.name, value = TRUE)
##DHEA-S(Dehydroepiandrosterone sulfate)
grep("DHEA", marker$Compound.name, value = TRUE)

######prediction model
temp_expression_data =
  expression_data[marker$variable_id,temp_sample_info$sample_id] %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

#####5 fold cross validation model
set.seed(seed = "123")
# length(unique(temp_sample_info$subject_id))
# idx =
# sample(1:length(unique(temp_sample_info$subject_id)), length(unique(temp_sample_info$subject_id)))
# 
# sample_idx =
# seq(1,length(idx), by = 28) %>% lapply(function(x){idx[x:(x+28)]}) %>% 
#   purrr::map(function(x){
#     x[!is.na(x)]
#   })
# 
# save(sample_idx, file = "sample_idx")
load("sample_idx")
sample_idx %>% boxplot()

result = vector(mode = "list", length = length(sample_idx))
for(i in 1:length(sample_idx)){
  cat(i, " ")
  validation_idx = sample_idx[[i]] %>%
    sort
  
  discovery_idx = sample_idx[-i] %>%
    unlist() %>%
    sort
  
  validation_data = t(temp_expression_data[,validation_idx])
  validation_sample_info = temp_sample_info[validation_idx,]
  
  discovery_data = t(temp_expression_data[,discovery_idx])
  discovery_sample_info = temp_sample_info[discovery_idx,]
  
  discovery_sample_info$timepoint[discovery_sample_info$timepoint == "2nd trimester"] = 0
  discovery_sample_info$timepoint[discovery_sample_info$timepoint == "3rd trimester"] = 1
  discovery_sample_info$timepoint = as.numeric(discovery_sample_info$timepoint)
  validation_sample_info$timepoint[validation_sample_info$timepoint == "2nd trimester"] = 0
  validation_sample_info$timepoint[validation_sample_info$timepoint == "3rd trimester"] = 1
  validation_sample_info$timepoint = as.numeric(validation_sample_info$timepoint)
  library(randomForest)
  # Fitting Random Forest to the train dataset
  set.seed(120)  # Setting seed
  classifier_RF = randomForest(x = as.matrix(discovery_data),
                               y = discovery_sample_info$timepoint,
                               ntree = 500)
  
  y_pred = predict(object = classifier_RF, newdata = as.matrix(validation_data))
  result[[i]] = 
    data.frame(sample_id = validation_sample_info$sample_id,
               true = validation_sample_info$timepoint,
               pred = y_pred,
               fold = i)
}

result = 
  do.call(rbind, result) %>% 
  as.data.frame()

library(pROC)
pROC_obj <- pROC::roc(result$true,result$pred,
                      smoothed = TRUE,
                      # arguments for ci
                      ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                      # arguments for plot
                      plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                      print.auc=TRUE, show.thres=TRUE)


sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")

plot(sens.ci, type="bars")

library(plotROC)

rocplot <- 
  ggplot(result, aes(m = pred, d = true))+ 
  geom_roc(n.cuts=20,labels=FALSE)
rocplot
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill = "pink") +
  base_theme

library(gghalves)
plot = 
result %>%
  dplyr::mutate(true = as.character(true)) %>%
  dplyr::mutate(true =
                  case_when(
                    true == "1" ~ "3rd trimester",
                    true == "0" ~ "2nd trimester"
                  )) %>% 
  ggplot(aes(true, pred)) +
  gghalves::geom_half_boxplot(aes(color = true),
                              show.legend = FALSE) +
  gghalves::geom_half_violin(aes(color = true), side = "r", 
                             show.legend = FALSE) +
  gghalves::geom_half_point(aes(fill = true),
                            shape = 21,
                            size = 3, alpha = 0.5,
                            side = "r", show.legend = FALSE) +
  # gghalves::geom_half_dotplot(side = "r", aes(fill = true), binwidth = 0.5 /
  #                               30, show.legend = FALSE) +
  scale_color_manual(values = trimester_color) +
  scale_fill_manual(values = trimester_color) +
  labs(x = "", y = "Predicetd value") +
  # scale_x_discrete(expand = expansion(mult = 10)) +
  base_theme
plot
ggsave(plot, filename = "predicetd_ga.pdf", width = 8, height = 7)
  

























