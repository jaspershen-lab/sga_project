sxtTools::setwd_project()
setwd("data/metabolomics/data_cleaning/POS_NEG/")
rm(list=ls())

library(tidyverse)
peak_table_pos = data.table::fread("../POS/expression_data_pos.csv")
peak_table_neg = data.table::fread("../NEG/expression_data_neg.csv")

sample_info_pos = data.table::fread("../POS/sample_info_pos.csv")
sample_info_neg = data.table::fread("../NEG/sample_info_neg.csv")

dim(peak_table_pos)

setdiff(colnames(peak_table_pos), colnames(peak_table_neg))
setdiff(colnames(peak_table_neg), colnames(peak_table_pos))

sample_info_pos =
  sample_info_pos %>% 
  as.data.frame() %>% 
  dplyr::mutate(sample.name = stringr::str_replace(sample.name, "posi\\_", ""))

sample_info_neg =
  sample_info_neg %>% 
  as.data.frame() %>% 
  dplyr::mutate(sample.name = stringr::str_replace(sample.name, "neg\\_", ""))

sum(duplicated(sample_info_pos$sample.name))
sum(duplicated(sample_info_neg$sample.name))

peak_table_pos = as.data.frame(peak_table_pos)
peak_table_neg = as.data.frame(peak_table_neg)

colnames(peak_table_pos) =
  colnames(peak_table_pos) %>% 
  stringr::str_replace("posi\\_", "")

colnames(peak_table_neg) =
  colnames(peak_table_neg) %>% 
  stringr::str_replace("neg\\_", "")

sample_info_pos$sample.name == colnames(peak_table_pos)[-c(1:3)]
sample_info_neg$sample.name == colnames(peak_table_neg)[-c(1:3)]

intersect_sample_id = intersect(sample_info_pos$sample.name,
                                sample_info_neg$sample.name)

sample_info_pos =
  sample_info_pos %>%
  dplyr::filter(sample.name %in% intersect_sample_id) %>% 
  dplyr::arrange(sample.name)

sample_info_neg =
  sample_info_neg %>%
  dplyr::filter(sample.name %in% intersect_sample_id) %>% 
  dplyr::arrange(sample.name)

peak_table_pos = 
  peak_table_pos[,c("name", "mz", "rt", sample_info_pos$sample.name)]

peak_table_neg = 
  peak_table_neg[,c("name", "mz", "rt", sample_info_neg$sample.name)]

dim(peak_table_pos)
dim(peak_table_neg)

colnames(peak_table_pos)==colnames(peak_table_neg)

peak_table=
  rbind(peak_table_pos,
        peak_table_neg)

sample_info_pos=
  sample_info_pos %>% 
  dplyr::rename(injection.order.pos = injection.order,
                batch.pos = batch) %>% 
  dplyr::select(-group)

sample_info_neg=
  sample_info_neg %>% 
  dplyr::rename(injection.order.neg = injection.order,
                batch.neg = batch) %>% 
  dplyr::select(-c(group, class))

sample_info = 
  sample_info_pos %>% 
  dplyr::left_join(sample_info_neg, by = c("sample.name"))

sample_info$sample.name == colnames(peak_table_pos)[-c(1:3)]
sample_info$sample.name == colnames(peak_table_neg)[-c(1:3)]

sample_info =
  sample_info %>% 
  dplyr::rename(sample_id = sample.name)

variable_info =
  peak_table[,c(1:3)]

expression_data =
  peak_table[,-c(1:3)]

dim(expression_data)
dim(variable_info)
dim(sample_info)

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) = variable_info$name

variable_info = 
variable_info %>% 
  dplyr::rename(variable_id = name)

##reassign injection order
sample_info = 
sample_info %>% 
  dplyr::arrange(injection.order.pos) %>% 
  dplyr::mutate(injection.order = 1:nrow(sample_info)) %>% 
  dplyr::arrange(sample_id)

save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")
save(sample_info, file = "sample_info")

###log
temp_data <-
  log(expression_data + 1, 2)

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
             aes(color = injection.order,
                 shape = class)) +
  scale_color_gradient(low = alpha("red", 0.1), high = alpha("red", 1)) +
  ggrepel::geom_text_repel(aes(PC1, PC2, label = ifelse(class == "QC",injection.order, NA))) +
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

ggsave(plot, filename = "pca_score.pdf", width = 7, height = 7)








