ggplot(data = data, aes(x = as.character(x),
                        y = as.numeric(as.character(y)))) +
  
  # boxplot
  geom_boxplot(
    outlier.shape = NA,
    # position = position_dodge(width=0.5),
    aes(fill = as.character(fill)),
    lwd = 1.1,
    fatten = 1) +
  scale_fill_manual(values = c("colors"),
                    name = "name") +
  
  # scatter points
  geom_jitter(
    # position = position_jitter(width = 0.3),
    size = 2.0,
    aes(colour = as.character(fill))) +
  scale_color_manual(values = c("colors"),
                     name = "name") +
  
  # labels
  xlab("lab") +
  ylab("lab") +
  ggtitle("") +
  coord_cartesian(ylim = c(0,1)) +
  # labs(fill = "lab") +
  
  # theme
  theme(
    plot.title = element_text(hjust=0.5, face = "bold", size = 20),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.background = element_blank(), axis.line = element_line(color = "black", size = 1),
    panel.grid.major = element_line(colour = "grey"))
