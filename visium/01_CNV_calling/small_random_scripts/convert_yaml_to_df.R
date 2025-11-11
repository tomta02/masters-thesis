# convert yaml file to R dataframe, make heatmap to identify which markers 
# were found in which samples 

library(yaml)
library(tidyverse)

# load yaml file containing codex image thresholds
yml = read_yaml("/g/saka/Tatjana/analysis/codex/yamlfiles/thresholds.yaml")

th = yml$thresholds %>%
  map(unlist) %>%                   
  map_dfr(~as.data.frame(as.list(.x)), .id = "samples") %>% 
  mutate(across(-samples, as.numeric))   

rownames(th) = th$samples

# read in marker names as character vector
mk = readLines("/g/saka/Tatjana/analysis/codex/codex_markers.txt")
# marker names in "th" have whitespace and non-alphanumeric chars converted to "."
#NOW: do the same with marker list, so names match
mk = gsub("[^A-Za-z0-9]", ".", mk)

th = th %>%
  select(samples, any_of(mk))

# exclude samples LN0438_MAAFHY1_R1 and LN0193_1ITRJL_R2 for now
th = th[!rownames(th) %in% c("LN0438_MAAFHY1_R1", "LN0193_1ITRJL_R2"), ]


thl = th %>%
  pivot_longer(-samples, names_to = "variable", values_to = "value")

plt = ggplot(thl, aes(x = variable, y = samples, fill = !is.na(value))) +
  geom_tile(color = "grey80") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#3b8132")) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/g/saka/Tatjana/data/04_CODEX/20251009_thresholded_samples.png", plt, width = 10, height = 6, dpi = 300)
