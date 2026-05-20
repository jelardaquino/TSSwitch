library(ggplot2)
library(dplyr)

df <- data.frame(
    category = factor(c("Isoforms", "Genes", "Switches"),
                      levels = c("Isoforms", "Genes", "Switches")),
    value = c(96, 86, 109)
)

ggplot(df, aes(x = category, y = value, fill = category)) +
    geom_col() +
    geom_text(aes(label = value), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + # color-friendly palette
    labs(x = NULL, y = "number of switching features") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu/switching_features.png", width = 4, height = 3, dpi = 300)
