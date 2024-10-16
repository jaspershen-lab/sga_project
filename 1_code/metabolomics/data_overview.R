sxtTools::setwd_project()
rm(list=ls())
source("code/tools.R")

library(tidyverse)

sxtTools::setwd_project()
setwd("data/metabolomics/data_overview")

####load data
load("../data_preparation/metabolite/expression_data")
load("../data_preparation/metabolite/sample_info")
load("../data_preparation/metabolite/variable_info")

dim(variable_info)

###PCA to show the batch effect
##positive
###log
variable_info_pos = 
  variable_info %>% 
  dplyr::filter(stringr::str_detect(variable_id, "POS"))

expression_data_pos =
  expression_data[variable_info_pos$variable_id,]

temp_data <-
  log(expression_data_pos + 1, 2)

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

pca_object <-
  prcomp(x =
           temp_data, center = FALSE, scale. = FALSE)

library(ggfortify)

x <- pca_object$x

###only for subject samples
rownames(x) == sample_info$sample_id

x <- x[,1:2]

x <- data.frame(x,
                sample_info,
                stringsAsFactors = FALSE)

plot <-
  ggplot(x, aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_point(size = 4, 
             aes(color = as.character(batch.pos),
                 shape = class)) +
  ggsci::scale_color_lancet() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = paste("PC1 (", round(summary(pca_object)$importance[2, 1] * 100, 2), "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object)$importance[2, 2] * 100, 2), "%)", sep = "")
  ) 

plot

ggsave(plot, filename = "pca_score_pos.pdf", width = 7, height = 7)

###PCA to show the batch effect
##negative
###log
variable_info_neg = 
  variable_info %>% 
  dplyr::filter(stringr::str_detect(variable_id, "NEG"))

expression_data_neg =
  expression_data[variable_info_neg$variable_id,]

temp_data <-
  log(expression_data_neg + 1, 2)

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

pca_object <-
  prcomp(x =
           temp_data, center = FALSE, scale. = FALSE)

library(ggfortify)

x <- pca_object$x

###only for subject samples
rownames(x) == sample_info$sample_id

x <- x[,1:2]

x <- data.frame(x,
                sample_info,
                stringsAsFactors = FALSE)

plot <-
  ggplot(x, aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_point(size = 4, 
             aes(color = as.character(batch.neg),
                 shape = class)) +
  ggsci::scale_color_lancet() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = paste("PC1 (", round(summary(pca_object)$importance[2, 1] * 100, 2), "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object)$importance[2, 2] * 100, 2), "%)", sep = "")
  ) 

plot

ggsave(plot, filename = "pca_score_neg.pdf", width = 7, height = 7)


####combine them together for PCA show GA and SGA and AGA
###log
temp_data <-
  log(expression_data + 1, 2)

##remove QC
temp_sample_info = 
  sample_info %>% 
  dplyr::filter(class == "Subject") %>% 
  dplyr::filter(!is.na(ultra_GA))

temp_data <-
  log(expression_data[,temp_sample_info$sample_id] + 1, 2)

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x))/sd(x)
  })

pca_object <-
  prcomp(x =
           temp_data, center = FALSE, scale. = FALSE)

library(ggfortify)

x <- pca_object$x

###only for subject samples
rownames(x) == temp_sample_info$sample_id

x <- x[,1:2]

x <- data.frame(x,
                temp_sample_info,
                stringsAsFactors = FALSE)

library(ggfortify)

plot <-
  autoplot(object = pca_object, 
           frame = TRUE, frame.type = 'norm', 
           data = temp_sample_info, 
           colour = "timepoint",
           fill = "timepoint",
           size = 4,
           alpha = 0.7
           ) +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  scale_color_manual(values = trimester_color) +
  scale_fill_manual(values = trimester_color) +
  base_theme 

plot

ggsave(plot, filename = "pca_score_for_sample_ga.pdf", width = 9, height = 7)


plot <-
  autoplot(object = pca_object, 
           frame = FALSE, 
           # frame.type = 'norm', 
           data = temp_sample_info, 
           colour = "sample_GA",
           # fill = "timepoint",
           size = 4,
           alpha = 0.7
  ) +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  # scale_color_manual(values = trimester_color) +
  # scale_fill_manual(values = trimester_color) +
  base_theme 

plot

ggsave(plot, filename = "pca_score_for_sample_ga2.pdf", width = 9, height = 7)




plot <-
  autoplot(object = pca_object, 
           frame = TRUE, 
           frame.type = 'norm', 
           data = temp_sample_info %>% 
             dplyr::mutate(SGA = case_when(
               SGA == 1 ~ "SGA",
               SGA == 0 ~ "AGA"
             )), 
           colour = "SGA",
           fill = "SGA",
           shape = "timepoint",
           size = 4,
           alpha = 0.7
  ) +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  scale_color_manual(values = sga_color) +
  scale_fill_manual(values = sga_color) +
  base_theme 

plot

ggsave(plot, filename = "pca_score_for_sga.pdf", width = 9, height = 7)










