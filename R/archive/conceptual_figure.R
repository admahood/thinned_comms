# make a conceptual figure
library(tidyverse)
library(geomtextpath)


data <- data.frame(x = 1:10, 
                   native = c(33,33,34,36, 37, 42, 47, 48, 51, 51),
                   invasive = c(0, 0, 1, 2, 1.5, 2, 3.5, 4, 6, 7)) |>
  mutate(total = native + invasive) |>
  pivot_longer(cols = c(native, invasive, total)) |>
  mutate(value = value/2)

ggplot(data, aes(x, value, color = name)) +
  geom_line(lwd=2) +
  geom_labelvline(xintercept = 2, lty=2, col ='red',# lwd = 2, 
                 label = 'Treatment', hjust = .3) +
  geom_labelhline(yintercept = 30, label = 'Maximum Cover Post-Treatment', lty = 3) +
  geom_labelhline(yintercept = 16.5, label = 'Maximum Cover Pre-Treatment', hjust = 0.7, lty = 3) +
  theme_classic() +
  geom_labelsegment(x=3, xend = 3, y = 18, yend = 29.5, arrow =arrow(ends = 'both'),
                    label = 'Avail. Res.', color = 'grey30') +
  xlab("Time") +
  ylab("Cover") +
  theme(legend.position = c(1, .3),
        legend.justification = c(1,0), 
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black'))

ggsave('out/conceptual_figure.png', width = 7, height  = 5, bg = 'white')
