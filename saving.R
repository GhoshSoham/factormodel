library(ggplot2)

re = matrix(0, 10, 21)
for(i in 1:10) {
  re[i, ] = output[[i]]
}
re1 = round(re[, 1:7], 6)
re2 = round(re[, 8:14], 6)
re3 = round(re[, 15:21], 6)
re1 = as.data.frame(re1)
re2 = as.data.frame(re2)
re3 = as.data.frame(re3)
colnames(re1) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")
colnames(re2) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")
colnames(re3) = c("mean", "var", "median", "0.05th", "0.25th", "0.75th", "0.95th")


write.csv(re1, file = "n250d20qt5q6.csv")
write.csv(re2, file = "n250d20qt5q7.csv")
write.csv(re3, file = "n250d20qt5q8.csv")

df = NULL

df1 = data.frame(mean = re1$mean, n = 5000, d = 1000, k=1, r = 1)
df2 = data.frame(mean = re2$mean, n = 5000, d = 1000, k=1, r = 5)
df3 = data.frame(mean = re3$mean, n = 5000, d = 1000, k=1, r = 10)
df = rbind(df1, df2, df3)

# Create the boxplot with adjusted ylim
boxplot(mean ~ r, data = df, 
        main = "Boxplot of Mean for Different n", 
        xlab = "n", 
        ylab = "Mean",
        col = "lightblue",
        border = "blue")  # Set ylim to include horizontal lines


# Determine the range of your data and include the extra values for the lines
y_range <- range(df$mean, 1, 5, 10)  # Include 1, 5, 10 in the range

# Create the boxplot with adjusted ylim
boxplot(mean ~ r, data = df, 
        main = "Boxplot of Mean for Different r", 
        xlab = "r", 
        ylab = "Mean",
        col = "lightblue",
        border = "blue",
        ylim = y_range)  # Set ylim to include horizontal lines

# Add horizontal lines at y = 1, 5, and 10
abline(h = c(1, 5, 10), col = "red", lty = 2)  # Red dashed lines

bp <- ggplot(df, mapping = aes(x=n, y=mean)) + 
  geom_boxplot(aes(fill=n))

boxplot = bp + facet_grid(. ~ n) + ggtitle("k = 2") + special_theme

df$n = as.factor(df$n)
special_theme = theme(plot.title = element_text(size = 15,hjust = .5, color = "black",
                                                vjust = .5),
)

#df$q = as.numeric(df$q)
bp <- ggplot(df, mapping = aes(x=n, y=mean, group=q)) + 
  geom_boxplot(aes(fill=q))
boxplot = bp + facet_grid(. ~ n) + ggtitle("k = 2") + special_theme

# Specify the file path and name with the unknown integer
file_path <- "/Users/sohamghosh/Desktop/factormodel/Results_homo_0.001/d30qt10.png"

ggsave(file_path, plot = boxplot, width = 6, height = 4, dpi = 300)



re1 = read.csv(file = "/Users/sohamghosh/Desktop/factormodel/n200d1000k3r2.csv", row.names = NULL)
oddrows = seq(3, 200, by = 4)
re1 = re1[oddrows, 2]
mean(re1)
sd(re1)


library(ggplot2)

# Example simulated results
set.seed(123)
n_vals <- c(100, 200)
d_vals <- c(500, 1000)
k_vals <- c(1, 3)

# Generate data frame
results_list <- list()
idx <- 1
for (n in n_vals) {
  for (d in d_vals) {
    for (k in k_vals) {
      df = read.csv(file = paste0("/Users/sohamghosh/Desktop/factormodel/n", n, "d", d, "k", k, "r2.csv"), row.names = NULL)
      dfm = df[seq(1, 200, by = 4), 2]
      fm = df[seq(3, 200, by = 4), 2]
      results_list[[idx]] <- data.frame(
        value = c(dfm, fm),  # 10 results per method
        Method = rep(c("DFM", "FM"), each = 50),
        n = n,
        d = d,
        k = k
      )
      idx <- idx + 1
    }
  }
}
data <- do.call(rbind, results_list)

# Create labels and treat k as a facet variable
data$combination <- paste0("n = ", data$n, ", d = ", data$d)
data$k <- factor(data$k, levels = c(1, 3), labels = c("k = 1", "k = 3"))

# Plot
ggplot(data, aes(x = Method, y = value, fill = Method)) +
  geom_boxplot() +
  facet_grid(k ~ combination, scales = "free") +
  labs(
    title = "Comparison of DFM vs FM across (n, d) settings by k",
    x = "Method", y = "Subspace Estimation Error"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(hjust = 1),
    strip.text = element_text(size = 10)
  )

re1 = read.csv(file = "/Users/sohamghosh/Desktop/factormodel/Results_homo_0.001/n1000d30qt10q11.csv")
re2 = read.csv(file = "/Users/sohamghosh/Desktop/factormodel/Results_homo_0.001/n1000d30qt10q12.csv")
re3 = read.csv(file = "/Users/sohamghosh/Desktop/factormodel/Results_homo_0.001/n1000d30qt10q13.csv")

re1 = read.csv(file = "/Users/sohamghosh/Desktop/factormodel/Results_homo_0.001/n5000d30qt10q11.csv")
re2 = read.csv(file = "/Users/sohamghosh/Desktop/factormodel/Results_homo_0.001/n5000d30qt10q12.csv")
re3 = read.csv(file = "/Users/sohamghosh/Desktop/factormodel/Results_homo_0.001/n5000d30qt10q13.csv")


