library(tidyverse)
library(cmdstanr)
library(ggplot2)
library(latex2exp)

specimens = c("VB1-v2", "VB2-v1", "VB3-v1", "VB4-v1")

read_laplace_csv = function(file) {
  lines = readLines(file)
  starts_with_pound = lapply(lines, function(line) {
    return(substr(line, 1, 1) == "#")
  }) %>%
    unlist
  
  #Optimization terminated with error:
  #Line search failed to achieve a sufficient decrease, no more progress can be made
  
  skip = which(starts_with_pound) %>% tail(1)
  n_max = which(starts_with_pound) %>% tail(1) - which(starts_with_pound) %>% tail(2) %>% head(1) - 2
  
  list(df_opt = read_csv(file, comment = "#", n_max = n_max, col_types = cols(.default = col_double())),
       df_laplace = read_csv(file, comment = "#", skip = skip, col_types = cols(.default = col_double())))
}

df2 = read_laplace_csv(paste0("VB-LA_outputs_VB2-v1/output_VB2-v1_45N_init-cij-tbc1_1.csv"))$df_laplace

df = lapply(specimens, function(specimen) {
  lapply(seq(10, 70, by = 5), function(N) {
    df_laplace = read_laplace_csv(paste0("VB-LA_outputs_", specimen, "/output_", specimen, "_", N, "N_init-cij-tbc1_1.csv"))$df_laplace
    
    tryCatch({
      return(df_laplace %>%
          filter(log_p != 0.0) %>%
            mutate(modes = N,
                   specimen = specimen))
    }, error = function(e) {
      print(paste0("specimen: ", specimen, ", N = ", N))
      return(NULL)
    })
  }) %>% bind_rows
}) %>% bind_rows

reference = tibble(c11_coat = 0.3,
                   c12_coat = 0.05,
                   c13_coat = 0.08,
                   c33_coat = 0.4,
                   c44_coat = 0.10,
                   c66_coat = 0.125) %>%
  mutate(c66_coat = (c11_coat - c12_coat) / 2.0) %>%
  rename("C[11]" = "c11_coat",
         "C[12]" = "c12_coat",
         "C[13]" = "c13_coat",
         "C[33]" = "c33_coat",
         "C[44]" = "c44_coat",
         "C[66]" = "c66_coat") %>%
  gather(param, value) %>%
  mutate(value = 100 * value)

df_khat = df %>%
  group_by(modes, specimen) %>%
  mutate(log_ratio = log_p - log_g) %>%
  select(modes, specimen, log_ratio) %>%
  summarize(khat = loo::psis(log_ratio)$diagnostics$pareto_k)

df %>%
  select(specimen, modes, c11_coat, c12_coat, c13_coat, c33_coat, c44_coat, c66_coat) %>%
  rename("C[11]" = "c11_coat",
         "C[12]" = "c12_coat",
         "C[13]" = "c13_coat",
         "C[33]" = "c33_coat",
         "C[44]" = "c44_coat",
         "C[66]" = "c66_coat") %>%
  gather(param, value, -modes, -specimen) %>%
  group_by(modes, specimen, param) %>%
  mutate(value = value * 100) %>%
  summarize(ql = quantile(value, 0.1),
            qh = quantile(value, 0.9),
            median = median(value)) %>%
  ggplot(aes(modes, median)) +
  theme_set(theme_bw(base_size = 10)) +
  geom_point(aes(color = param, group = param), position = position_dodge(width = 2), size = 1.0) +
  #geom_ribbon(aes(ymin = ql, ymax = qh, group = param, fill = param), alpha = 0.15) +
  geom_errorbar(aes(ymin = ql, ymax = qh, group = param, color = param), width = 0, position = position_dodge(width = 2)) +
  geom_hline(data = reference, aes(yintercept = value, color = param), size = 0.5, linetype = "dashed") +
  ylab("parameter value") +
  theme(text = element_text(size = 20),
        panel.background = element_rect(linetype = "solid"),
        panel.grid.major = element_line(size = 0.35, linetype = 'solid', colour = "gray80"),
        axis.ticks.length=unit(-0.15, "cm"), # negative value puts tick marks to the inside
        axis.text.x.bottom = element_text(margin = margin(t = 0.3, b = 0.0, unit = "cm")), #adjusts the margins after adjusting tick size
        axis.text.y = element_text(margin = margin(r = 0.3, l = 0.0, unit = "cm")), #adjusts the margins after adjusting tick size
        plot.margin = unit(c(0.5,0.4,0.2,0.0),"cm"), #specify margins starting from the top, then right, bottom and left
        axis.text.x.top = element_blank(), # remove tick labels from top axes
        axis.text.y.right = element_blank(), # remove tick labels from right axes
        axis.title.x.top = element_blank(), # remove title from top axes
        axis.title.y.right = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.placement.x = "outside") +
  facet_wrap(. ~ specimen) +
  coord_cartesian(ylim = c(0, 50))# remove title from right axes)

df_khat %>%
  mutate(quality = ifelse(khat < 0.5, "good", ifelse(khat < 0.7, "okay", "bad")),
         inrange = ifelse(khat < 1.09, 16, 4),
         khat = pmin(khat, 1.09)) %>%
  ggplot() +
  theme_set(theme_bw(base_size = 10)) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_point(aes(modes, khat, color = quality, shape = inrange), size = 1.5, stroke = 1.0) +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0), limits = c(0.0, 1.1)) +
  ylab(TeX('$\\hat{k}')) +
  scale_shape_identity() +
  scale_colour_manual(values = c("red", "green2", "blue")) +
  theme(text = element_text(size = 20),
        legend.position = "none",
        panel.background = element_rect(linetype = "solid"),
        panel.grid.major = element_line(size = 0.35, linetype = 'solid', colour = "gray80"),
        axis.ticks.length=unit(-0.15, "cm"), # negative value puts tick marks to the inside
        axis.text.x.bottom = element_text(margin = margin(t = 0.3, b = 0.0, unit = "cm")), #adjusts the margins after adjusting tick size
        axis.text.y = element_text(margin = margin(r = 0.3, l = 0.0, unit = "cm")), #adjusts the margins after adjusting tick size
        plot.margin = unit(c(0.5,0.4,0.2,0.0),"cm"), #specify margins starting from the top, then right, bottom and left
        axis.text.x.top = element_blank(), # remove tick labels from top axes
        axis.text.y.right = element_blank(), # remove tick labels from right axes
        axis.title.x.top = element_blank(), # remove title from top axes
        axis.title.y.right = element_blank(),
        axis.title.x = element_blank(),
        strip.background.x = element_blank(),
        strip.placement.x = "outside") +
  facet_wrap(~ specimen, scales = "free_x")
