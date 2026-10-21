#load required packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("purrr")
install.packages("broom")
install.packages("ggplot2")
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)

#paired t test for one protein
paired_data <- assay_df %>% 
  select(Experiment, Sample, ABI3BP) %>%
  pivot_wider( names_from = Sample,values_from = ABI3BP ) %>% 
  drop_na(Control, Treatment) 

paired_t_test <- t.test( paired_data$Treatment,paired_data$Control, 
                         paired = TRUE, alternative = "two.sided" ) 
paired_t_test

#how to extract t.test results
result <- paired_t_test$p.value
result

confidence_interval <- paired_t_test$conf.int
confidence_interval

assay_df_final <- assay_df[,-c("Extension control", "Incubation control", "Amplification control")]
assay_df_final <- subset(assay_df, 
                         select = -c(`Extension control`,`Incubation control`, `Amplification control`))

longer_df <- assay_df_final %>%
  pivot_longer(
    cols = -c(Experiment, Sample),
    names_to = "Protein",
    values_to = "Expression"
  )

wider <- longer_df %>%
  pivot_wider(
    id_cols = c(Experiment, Protein),
    names_from = Sample,
    values_from = Expression
  )

##conducting t test with Benjamini-Hochberg correction for multiple testing

ttest_results_bh <- assay_df_final %>%
  pivot_longer(
    cols = -c(Experiment, Sample),
    names_to = "Protein",
    values_to = "Expression"
  ) %>%
  pivot_wider(
    id_cols = c(Experiment, Protein),
    names_from = Sample,
    values_from = Expression
  ) %>%
  drop_na(Control, Treatment) %>%
  group_by(Protein) %>%
  nest() %>%
  mutate(
    test = map(
      data,
      ~ tryCatch(
        t.test(
          .x$Treatment,
          .x$Control,
          paired = TRUE
        ),
        error = function(e) NULL
      )
    ),
    result = map(test, ~ {
      if (is.null(.x)) {
        tibble(
          estimate = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_
        )
      } else {
        tidy(.x)
      }
    })
  ) %>%
  unnest(result) %>%
  ungroup() %>%
  mutate(
    p_adjusted = p.adjust(p.value, method = "BH")
  ) %>%
  select(
    Protein,
    estimate,
    statistic,
    p.value,
    p_adjusted,
    everything()
  )

ttest_results_bh
library(forcats)
#scatter plot
plot(ttest_results$p.value, ttest_results$p_adjusted)

#creating NPX means for treatment and cotnrol
control_means <- assay_df_final |>
  filter(Sample == "Control") |>
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

treatment_means <- assay_df_final |>
  filter(Sample == "Treatment") |>
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

#calculating difference in mean NPX
log2FC <- treatment_means - control_means

#converting to long format
mean_results <- data.frame(
  Protein = names(log2FC),
  log2FC = as.numeric(log2FC)
)

#joining to ttest_results to make volcano plot
volcano_df <- ttest_results |>
  left_join(mean_results, by = "Protein")

#adding column to show downregulated and upregulated proteins

volcano_df <- volcano_df |>
  mutate(
    regulation = case_when(
      log2FC > 0 & p_adjusted < 0.05 ~ "UPREGULATED",
      log2FC < 0 & p_adjusted < 0.05 ~ "DOWNREGULATED",
      .default = "NS"
    ),
    regulation = fct_relevel(
      regulation,
      "UPREGULATED",
      "DOWNREGULATED",
      "NS"
    )
  )

#creating vector of proteins i want labelled
prot <- c("LPO", "PTPRN2", "BANK1")  

#making a volcano plot
volcano_df |> 
  ggplot(aes(log2FC, -log10(p_adjusted), color = regulation))+
  geom_point(size = 2, 
             alpha = 0.7, 
             shape = 16)+


  geom_vline(
    xintercept = 0, 
    linetype = "dashed"
    )+
  

   geom_hline(
     yintercept = -log10(0.05), 
     linetype = "dashed"
     )+
  
##colour code according to upregulated, downregulated, NS
  scale_color_manual(
    values = c(
      "DOWNREGULATED" = "#fc8d59",
      "UPREGULATED" = "#91bfdb",
      "NS" = "grey70"
    ),
    name = ""
  )+
  
  geom_text_repel(
    data = volcano_df |> filter(Protein %in% prot),
    aes(label = Protein),
    box.padding = 0.5, 
    point.padding = 0.3, 
    max.overlaps = Inf,
    color =  "black",
    )+
 
   theme_classic()+
  
  labs(
    x = "Difference in mean NPX (Treatment − Control)",
    y = "-log10(P-adjusted value)"
  )+
  
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "top"
    )+


  guides(
    color = guide_legend(
      override.aes = list(size = 5)
      )
    )
