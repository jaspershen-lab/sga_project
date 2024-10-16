sxtTools::setwd_project()
rm(list=ls())
source("code/tools.R")
library(tidyverse)

sxtTools::setwd_project()
setwd("data/metabolomics/sga_biomarker")

####load data
load("../data_preparation/metabolite/expression_data")
load("../data_preparation/metabolite/sample_info")
load("../data_preparation/metabolite/variable_info")

####remove SGA which are NA
sample_info$group
sample_info = 
sample_info %>% 
  dplyr::filter(group != "QC") %>% 
  dplyr::filter(!is.na(timepoint)) %>% 
  dplyr::mutate(SGA = case_when(
    SGA == 1 ~ "SGA",
    SGA == 0 ~ "AGA"
  ))

expression_data =
  expression_data[,sample_info$sample_id]

###remove some metabolites
variable_info = 
variable_info %>% 
  dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2", "msDatabase_rplc0.0.2",
                                "nistDatabase0.0.2"))

expression_data = expression_data[variable_info$variable_id,]

###find the biomakers which are different in AGA and SGA
###trimester 2
t2_fc_p = 
  expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
purrr::map(function(x){
  x = as.numeric(x)
  ###adjusted ga
  control_idx = which(sample_info$timepoint == "2nd trimester" & sample_info$group == "AGA")
  case_idx = which(sample_info$timepoint == "2nd trimester" & sample_info$group == "SGA")
  x1 = x[control_idx]
  x2 = x[case_idx]
  ga1 = sample_info$sample_GA[control_idx]
  ga2 = sample_info$sample_GA[case_idx]
  lm_test1 = lm(formula = c(x1) ~ c(ga1))
  lm_test2 = lm(formula = c(x2) ~ c(ga2))
  x1 = lm_test1$residuals
  x2 = lm_test2$residuals
  p.value = 
  wilcox.test(x1, x2)$p.value
  fc = mean(x[case_idx])/mean(x[control_idx])
  c(fc = fc, p = p.value)
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(p.adjust = p.adjust(p, method = "BH"))


rownames(t2_fc_p) == variable_info$variable_id
t2_fc_p =
  data.frame(variable_info, t2_fc_p)

###trimester 3
t3_fc_p = 
  expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    x = as.numeric(x)
    ###adjusted ga
    control_idx = which(sample_info$timepoint == "3rd trimester" & sample_info$group == "AGA")
    case_idx = which(sample_info$timepoint == "3rd trimester" & sample_info$group == "SGA")
    x1 = x[control_idx]
    x2 = x[case_idx]
    ga1 = sample_info$sample_GA[control_idx]
    ga2 = sample_info$sample_GA[case_idx]
    lm_test1 = lm(formula = c(x1) ~ c(ga1))
    lm_test2 = lm(formula = c(x2) ~ c(ga2))
    x1 = lm_test1$residuals
    x2 = lm_test2$residuals
    p.value = 
      wilcox.test(x1, x2)$p.value
    fc = mean(x[case_idx])/mean(x[control_idx])
    c(fc = fc, p = p.value)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(p.adjust = p.adjust(p, method = "BH"))

rownames(t3_fc_p) == variable_info$variable_id
t3_fc_p =
  data.frame(variable_info, t3_fc_p)

sum(t2_fc_p$p.adjust < 0.05)
sum(t3_fc_p$p.adjust < 0.05)

t2_fc_p %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::select(Compound.name, fc, p.adjust, SS, Total.score) %>% 
  dplyr::arrange(Compound.name)

t3_fc_p %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::select(Compound.name, fc, p.adjust, SS, Total.score) %>% 
  dplyr::arrange(Compound.name)

t2_fc_p = 
t2_fc_p %>% 
  dplyr::filter(Compound.name != "5-Aminopentanoic acid") %>% 
  dplyr::filter(!variable_id %in% c("M373T570_POS", "M407T569_NEG"))

t3_fc_p = 
  t3_fc_p %>% 
  dplyr::filter(Compound.name != "5-Aminopentanoic acid") %>% 
  dplyr::filter(!variable_id %in% c("M373T570_POS", "M407T569_NEG"))

sum(t2_fc_p$p.adjust < 0.05)
sum(t3_fc_p$p.adjust < 0.05)

#####volcano plot

plot1 = 
volcano_plot(
  fc = t2_fc_p$fc,
  p_value = t2_fc_p$p.adjust,
  p.cutoff = 0.05,
  fc.cutoff = 1,
  text = TRUE, 
  variable_id = t2_fc_p$Compound.name, 
  point.size = 4
)
plot1
ggsave(plot1, filename = "volocano_t2_plot.pdf", width = 7, height = 7)

plot2 = 
volcano_plot(
  fc = t3_fc_p$fc,
  p_value = t3_fc_p$p.adjust,
  p.cutoff = 0.05,
  fc.cutoff = 1,
  text = TRUE, 
  variable_id = t3_fc_p$Compound.name
)
plot2
ggsave(plot2, filename = "volocano_t3_plot.pdf", width = 7, height = 7)

###output result
t2_fc_p %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::select(Compound.name, fc, p.adjust) %>% 
  dplyr::arrange(Compound.name) %>% 
  dplyr::left_join(t3_fc_p %>% 
                     dplyr::filter(p.adjust < 0.05) %>% 
                     dplyr::select(Compound.name, fc, p.adjust) %>% 
                     dplyr::arrange(Compound.name), by = "Compound.name")

###show the overlap between markers in two different trimester
temp_data =
  rbind(
    data.frame(
      t2_fc_p %>% dplyr::select(variable_id, Compound.name, fc:p.adjust),
      trimester = "2"
    ),
    data.frame(
      t3_fc_p %>% dplyr::select(variable_id, Compound.name, fc:p.adjust),
      trimester = "3"
    )
  ) %>%
  dplyr::mutate(
    fc = log(fc, 2),
    p.adjust = -log(p.adjust, 10),
    class = case_when(
      p.adjust > -log(0.05, 10) & fc > log(1, 2) ~ "Increase",
      p.adjust > -log(0.05, 10) &
        fc < log(1 / 1, 2) ~ "Decrease",
      TRUE ~ "No"
    )
  )

segment_data =
  temp_data %>% 
  dplyr::select(-c(Compound.name, class, p, p.adjust)) %>% 
  tidyr::pivot_wider(names_from = trimester, values_from = c(fc)) %>% 
  dplyr::rename(y = "2", yend = "3") %>% 
  dplyr::mutate(x = "2", xend = "3") %>% 
  dplyr::left_join(temp_data %>% dplyr::filter(trimester == "2") %>% dplyr::select(variable_id, class),
                   by = "variable_id") %>% 
  dplyr::rename(class2 = class) %>% 
  dplyr::left_join(temp_data %>% dplyr::filter(trimester == "3") %>% dplyr::select(variable_id, class),
                   by = "variable_id") %>% 
  dplyr::rename(class3 = class) %>% 
  dplyr::filter(class2 != "No" | class3 != "No") %>% 
  dplyr::mutate(class = case_when(
    class2 == class3 ~ "Same",
    class2 != class3 ~ "Difference"
  )) %>% 
  dplyr::left_join(variable_info[,c("variable_id", "Compound.name")], by = "variable_id")


library(ggnewscale)
plot = 
ggplot(data = temp_data, aes(trimester, fc)) +
  geom_hline(yintercept = 0) +
  geom_segment(
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend,
      color = class
    ),
    data = segment_data,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Same" = "red",
                                "Different" = "grey")) +
  ggnewscale::new_scale_color() +
  geom_point(
    shape = 16,
    aes(color = class,
        size = p.adjust),
    show.legend = TRUE,
    alpha = 0.7,
  ) +
  scale_color_manual(
    values = c(
      "Increase" = ggsci::pal_aaas()(10)[2],
      "Decrease" = ggsci::pal_aaas()(10)[1],
      "No" = "black"
    )
  ) +
  ggrepel::geom_text_repel(
    aes(
      x = trimester,
      y = fc,
      label = Compound.name,
      color = class
    ),
    data = temp_data %>% dplyr::filter(trimester == 2, class != "No"),
    direction = "y",
    hjust = "right",
    nudge_x = -0.2
  ) +
  ggrepel::geom_text_repel(
    aes(
      x = trimester,
      y = fc,
      label = Compound.name,
      color = class
    ),
    data = temp_data %>% dplyr::filter(trimester == 3, class != "No"),
    direction = "y",
    hjust = "left",
    nudge_x = 0.2
  ) +
  scale_x_discrete(expand = expansion(mult = c(2))) +
  base_theme +
  labs(x = "Trimester", y = "log2 (Fold change)")

plot
ggsave(plot, filename = "biomarker.pdf", width = 7, height = 12)


final_result = 
  t2_fc_p %>% 
  dplyr::left_join(t3_fc_p %>% dplyr::select(variable_id, fc: p.adjust), by = "variable_id") %>% 
  dplyr::rename(fc.t2 = fc.x,
                p.t2 = p.x,
                p.adjust.t2 = p.adjust.x,
                fc.t3 = fc.y,
                p.t3 = p.y,
                p.adjust.t3 = p.adjust.y) %>% 
  dplyr::select(variable_id, fc.t2:p.adjust.t3, everything())

increase_marker = 
final_result %>% 
  dplyr::filter(fc.t2 > 1 & fc.t3 > 1) %>% 
  dplyr::filter(p.adjust.t2 < 0.05 | p.adjust.t3 < 0.05)

decrease_marker = 
final_result %>% 
  dplyr::filter(fc.t2 < 1 & fc.t3 < 1) %>% 
  dplyr::filter(p.adjust.t2 < 0.05 | p.adjust.t3 < 0.05)

marker = 
  rbind(increase_marker, decrease_marker)

dim(marker)
dim(increase_marker)


###show the overlap between markers in two different trimester
temp_data =
  rbind(
    data.frame(
      t2_fc_p %>% dplyr::select(variable_id, Compound.name, fc:p.adjust),
      trimester = "2"
    ),
    data.frame(
      t3_fc_p %>% dplyr::select(variable_id, Compound.name, fc:p.adjust),
      trimester = "3"
    )
  ) %>%
  dplyr::mutate(
    fc = log(fc, 2),
    p.adjust = -log(p.adjust, 10),
    class = case_when(
      p.adjust > -log(0.05, 10) & fc > log(1, 2) ~ "Increase",
      p.adjust > -log(0.05, 10) &
        fc < log(1 / 1, 2) ~ "Decrease",
      TRUE ~ "No"
    )
  ) %>% 
  dplyr::filter(variable_id %in% marker$variable_id)

segment_data =
  temp_data %>% 
  dplyr::select(-c(Compound.name, class, p, p.adjust)) %>% 
  tidyr::pivot_wider(names_from = trimester, values_from = c(fc)) %>% 
  dplyr::rename(y = "2", yend = "3") %>% 
  dplyr::mutate(x = "2", xend = "3") %>% 
  dplyr::left_join(temp_data %>% dplyr::filter(trimester == "2") %>% dplyr::select(variable_id, class),
                   by = "variable_id") %>% 
  dplyr::rename(class2 = class) %>% 
  dplyr::left_join(temp_data %>% dplyr::filter(trimester == "3") %>% dplyr::select(variable_id, class),
                   by = "variable_id") %>% 
  dplyr::rename(class3 = class) %>% 
  dplyr::filter(class2 != "No" | class3 != "No") %>% 
  dplyr::mutate(class = case_when(
    class2 == class3 ~ "Same",
    class2 != class3 ~ "Difference"
  )) %>% 
  dplyr::left_join(variable_info[,c("variable_id", "Compound.name")], by = "variable_id")


library(ggnewscale)
plot = 
  ggplot(data = temp_data, aes(trimester, fc)) +
  geom_hline(yintercept = 0) +
  geom_segment(
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend,
      color = class
    ),
    data = segment_data,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Same" = "red",
                                "Different" = "grey")) +
  ggnewscale::new_scale_color() +
  geom_point(
    shape = 16,
    aes(color = class,
        size = p.adjust),
    show.legend = TRUE,
    alpha = 0.7,
  ) +
  scale_color_manual(
    values = c(
      "Increase" = ggsci::pal_aaas()(10)[2],
      "Decrease" = ggsci::pal_aaas()(10)[1],
      "No" = "black"
    )
  ) +
  ggrepel::geom_text_repel(
    aes(
      x = trimester,
      y = fc,
      label = Compound.name,
      color = class
    ),
    data = temp_data %>% dplyr::filter(trimester == 2, class != "No"),
    direction = "y",
    hjust = "right",
    nudge_x = -0.2
  ) +
  ggrepel::geom_text_repel(
    aes(
      x = trimester,
      y = fc,
      label = Compound.name,
      color = class
    ),
    data = temp_data %>% dplyr::filter(trimester == 3, class != "No"),
    direction = "y",
    hjust = "left",
    nudge_x = 0.2
  ) +
  scale_x_discrete(expand = expansion(mult = c(2))) +
  base_theme +
  labs(x = "Trimester", y = "log2 (Fold change)")

plot
ggsave(plot, filename = "biomarker2.pdf", width = 7, height = 12)


###look like that metronidazole are very high in the SGA group. this drug is used
###to treat malaria
temp_data = 
sample_info %>% 
  dplyr::select(subject_id, matmalaria, group) %>% 
  dplyr::distinct(subject_id, .keep_all = TRUE)

table(temp_data$matmalaria, temp_data$group)

chisq.test(x = c(61, 60), y = c(7, 9))
chisq.test(matrix(c(61,7,60, 9), nrow = 2))

######prediction model
####trimester 2
marker_t2 = 
marker %>% 
  dplyr::filter(p.adjust.t2 < 0.05)

temp_sample_info = 
  sample_info %>% 
  dplyr::filter(timepoint == "2nd trimester")

temp_expression_data =
  expression_data[marker_t2$variable_id,temp_sample_info$sample_id] %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

#####5 fold cross validation model
set.seed(seed = "123")
length(unique(temp_sample_info$subject_id))
idx =
sample(1:length(unique(temp_sample_info$subject_id)), length(unique(temp_sample_info$subject_id)))

sample_idx =
seq(1,length(idx), by = 22) %>% lapply(function(x){idx[x:(x+21)]})

save(sample_idx, file = "sample_idx")
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
  
  discovery_sample_info$group[discovery_sample_info$group == "SGA"] = 1
  discovery_sample_info$group[discovery_sample_info$group == "AGA"] = 0
  discovery_sample_info$group = as.numeric(discovery_sample_info$group)
  validation_sample_info$group[validation_sample_info$group == "SGA"] = 1
  validation_sample_info$group[validation_sample_info$group == "AGA"] = 0
  validation_sample_info$group = as.numeric(validation_sample_info$group)
  library(randomForest)
  # Fitting Random Forest to the train dataset
  set.seed(120)  # Setting seed
  classifier_RF = randomForest(x = as.matrix(discovery_data),
                               y = discovery_sample_info$group,
                               ntree = 500)
  
  y_pred = predict(object = classifier_RF, newdata = as.matrix(validation_data))
  result[[i]] = 
  data.frame(sample_id = validation_sample_info$sample_id,
             true = validation_sample_info$group,
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
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill = "pink") 












######prediction model
####trimester 3
marker_t3 = 
  marker %>% 
  dplyr::filter(p.adjust.t3 < 0.05)

temp_sample_info = 
  sample_info %>% 
  dplyr::filter(timepoint == "3rd trimester")

temp_expression_data =
  expression_data[marker_t3$variable_id,temp_sample_info$sample_id] %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

#####5 fold cross validation model
set.seed(seed = "123")
length(unique(temp_sample_info$subject_id))
idx =
sample(1:length(unique(temp_sample_info$subject_id)), length(unique(temp_sample_info$subject_id)))

sample_idx =
seq(1,length(idx), by = 22) %>% lapply(function(x){idx[x:(x+21)]})

save(sample_idx, file = "sample_idx")
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
  
  discovery_sample_info$group[discovery_sample_info$group == "SGA"] = 1
  discovery_sample_info$group[discovery_sample_info$group == "AGA"] = 0
  discovery_sample_info$group = as.numeric(discovery_sample_info$group)
  validation_sample_info$group[validation_sample_info$group == "SGA"] = 1
  validation_sample_info$group[validation_sample_info$group == "AGA"] = 0
  validation_sample_info$group = as.numeric(validation_sample_info$group)
  library(randomForest)
  # Fitting Random Forest to the train dataset
  set.seed(120)  # Setting seed
  classifier_RF = randomForest(x = as.matrix(discovery_data),
                               y = discovery_sample_info$group,
                               ntree = 500)
  
  y_pred = predict(object = classifier_RF, newdata = as.matrix(validation_data))
  result[[i]] = 
    data.frame(sample_id = validation_sample_info$sample_id,
               true = validation_sample_info$group,
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
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill = "pink") 










