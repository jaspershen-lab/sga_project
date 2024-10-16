sxtTools::setwd_project()
setwd("data/metabolomics/data_cleaning/NEG/")
rm(list=ls())

library(tidyverse)
peak_table_neg = data.table::fread("Peak_table_neg.csv")
colnames(peak_table_neg)

sample_info = readxl::read_xlsx("Metflow sample list.xlsx")

sample_info_neg = 
sample_info %>% 
  dplyr::rename(sample.name = "sample name", 
                injection.order = "injection order") %>% 
  dplyr::filter(stringr::str_detect(sample.name, "neg"))
sample_info_neg$sample.name
sample_info_neg$injection.order
peak_table_neg = as.data.frame(peak_table_neg)
sample_info_neg$batch[sample_info_neg$injection.order  >= 287] = 2
table(sample_info_neg$batch)

dim(sample_info_neg)
dim(peak_table_neg)

setdiff(colnames(peak_table_neg), sample_info_neg$sample.name)
setdiff(sample_info_neg$sample.name, colnames(peak_table_neg))

sample_info_neg = 
  sample_info_neg %>% 
  dplyr::filter(!sample.name %in% c(setdiff(sample_info_neg$sample.name, colnames(peak_table_neg))))
  
peak_table_neg = 
  peak_table_neg[,c("name", "mz", "rt", sample_info_neg$sample.name)]

sample_info_neg =
  as.data.frame(sample_info_neg)

###remove QCP
sample_info_neg = 
sample_info_neg %>% 
  dplyr::filter(class != "QCP")

peak_table_neg =
  peak_table_neg[,c("name", "mz", "rt", sample_info_neg$sample.name)]

write.csv(peak_table_neg, file = "peak_table_neg2.csv", row.names = FALSE)
write.csv(sample_info_neg, file = "sample_info_neg2.csv", row.names = FALSE)

library(metflow2)
object <-
  create_metflow_object(
    ms1.data = c("peak_table_neg2.csv"),
    sample.information = "sample_info_neg2.csv",
    path = "."
  )

object2 <- filter_peaks(
  object = object,
  min.fraction = 0.8,
  type = "all",
  according.to = "class",
  which.group = c("QC", "QCP")
)

object2 <- filter_peaks(
  object = object2,
  min.fraction = 0.5,
  type = "all",
  according.to = "class",
  which.group = c("Subject")
)

object2 <- filter_samples(object = object2,
                          min.fraction.peak = 0.5)

get_mv_plot_samples(object = object2, interactive = TRUE)

object2 <- metflow2::impute_mv(object = object2,
                               method = "minimum")

object3 <- normalize_data(object = object2, method = "mean")
object4 <- integrate_data(object = object3, method = "subject.mean")

####PCA plot
peak_table = object4@ms1.data[[1]]
sample_info = object4@sample.info

###log
temp_data <-
  log(peak_table[,-c(1:3)] + 1, 2)

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

rownames(x) == sample_info$sample.name

x <- x[,1:2]

x <- data.frame(x, 
               sample_info,
                stringsAsFactors = FALSE)

plot <-
  ggplot(x, aes(PC1, PC2)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_point(size = 4, 
             aes(color = as.character(batch),
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

###output peak table and sample info
expression_data_neg = object4@ms1.data[[1]]
sample_info_neg = object4@sample.info
write.csv(expression_data_neg, "expression_data_neg.csv", row.names = FALSE)
write.csv(sample_info_neg, "sample_info_neg.csv", row.names = FALSE)





