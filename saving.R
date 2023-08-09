library(ggplot2)

re = matrix(0, 50, 21)
for(i in 1:50) {
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


# write.csv(re1, file = "n250d30qt5q6.csv")
# write.csv(re2, file = "n250d30qt5q7.csv")
# write.csv(re3, file = "n250d30qt5q8.csv")

df = NULL

df1 = data.frame(mean = re1$mean, n = 250, d = 30, qt = 5, q = 6)
df2 = data.frame(mean = re2$mean, n = 250, d = 30, qt = 5, q = 7)
df3 = data.frame(mean = re3$mean, n = 250, d = 30, qt = 5, q = 8)
df = rbind(df, df1, df2, df3)

df$q = as.factor(df$q)
df$n = as.factor(df$n)
special_theme = theme(plot.title = element_text(size = 15,hjust = .5, color = "black",
                                                vjust = .5),
)
#df$q = as.numeric(df$q)
bp <- ggplot(df, mapping = aes(x=q, y=mean, group=q)) + 
  geom_boxplot(aes(fill=q))
boxplot = bp + facet_grid(. ~ n) + ggtitle("d = 20, qt = 10") + special_theme

# Specify the file path and name with the unknown integer
file_path <- "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Plots_homo/d20qt10.png"

ggsave(file_path, plot = boxplot, width = 6, height = 4, dpi = 300)



re1 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n250d20qt10q11.csv")
re2 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n250d20qt10q12.csv")
re3 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n250d20qt10q13.csv")

re1 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n1000d20qt10q11.csv")
re2 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n1000d20qt10q12.csv")
re3 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n1000d20qt10q13.csv")

re1 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n5000d20qt10q11.csv")
re2 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n5000d20qt10q12.csv")
re3 = read.csv(file = "/Users/sohamghosh/Desktop/orderinvariantfactormodelcode/Results_homo/n5000d20qt10q13.csv")


