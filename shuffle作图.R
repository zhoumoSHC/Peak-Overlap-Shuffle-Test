library(ggplot2)
library(readr)
library(dplyr)
library(grid)
# 将随机打乱(shuffle)后的peak重叠次数分布绘制为密度图
# 作为shuffle.sh分析的后续，需要输入的内容包括shuffle的结果，以及原始两个bed文件的overlap数量
# 并在最后自定义文件名
peak <- read_table('overlap_result.txt', col_names = "overlaps")
overlap_count <- 30 # 注意修改这个数值为重叠的peak数

# 自适应 x 轴范围（左右各扩 10%，左界 >= 0）
x_min <- min(peak$overlaps, overlap_count)
x_max <- max(peak$overlaps, overlap_count)
range_expand <- (x_max - x_min) * 0.1
x_min <- max(0, x_min - range_expand)
x_max <- x_max + range_expand
brks <- pretty(c(x_min, x_max), n = 10)

# 计算密度用于确定 y 轴上限
dens  <- density(peak$overlaps, adjust = 2, from = 0)
y_top <- max(dens$y)

# 箭头长度：从 0 到 y_top 的固定比例（比如 20%）
arrow_height_ratio <- 0.2
y_start <- y_top * arrow_height_ratio
y_end   <- 0  # 坐标轴底端

gg <- ggplot(peak, aes(x = overlaps)) +
  geom_density(color = "black", size = 1, adjust = 2, from = 0) +
  labs(
    title = "Distribution of Random Overlaps Compared to Observed Overlaps",
    x = "Number of Overlaps",
    y = "Density"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = brks,
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.005, 0.06))
  ) +
  geom_segment(
    x = overlap_count, xend = overlap_count,
    y = y_start, yend = y_end,
    lineend = "round", linejoin = "round",
    linewidth = 1, colour = "#EC7014",
    arrow = arrow(length = unit(0.04, "snpc"))  # 固定比例箭头头
  )
gg
ggsave("name.pdf", gg,
       device = "pdf", width = 8.78, height = 5.68, units = "in")
