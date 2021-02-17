# 读入mpa格式的tsv物种表格文件
tax.all.count <- read.csv("tax_count_species.tsv", 
                          sep = "\t", header = T)

# 从全体数据中提取出“种”的信息
# species.count <- tax.all.count[grep("s__", tax.all.count[ , "tax"], invert = F), ]
# 或者不用单独提出“种”信息，直接用全部物种数据
species.count <- tax.all.count
# 把第一列去除成为行名
rownames(species.count) <- species.count[ , 1]
species.count <- species.count[ , 2:11]

# 分组信息
metadata <- data.frame(Group = c('non', 'non', 'non', 'hsct', 'hsct', 
                                 'hsct', 'non', 'non', 'non', 'hsct'), 
                       Series = c('01', '02', '03', '04', '05', 
                                  '06', '07', '08', '09', '10'),
                       Samples = c('BALF', 'BALF', 'BALF', 'BALF', 'BALF',
                                   'BALF', 'BALF', 'BALF', 'BALF', 'BALF'),
                       row.names = colnames(species.count))

# 导入相关的运算库
library(vegan)
library(amplicon)
library(pheatmap)
library(RColorBrewer)

# 计算α多样性指标
shannon.wiener <- diversity(t(species.count), index = "shannon")
simpson <- diversity(t(species.count), index = "simpson")
inverse.simpson <- diversity(t(species.count), index = "inv")
S <- specnumber(t(species.count))
pielou <- shannon.wiener/log(S)
# 建立α多样性数据框
alpha.diversity <- rbind(shannon.wiener, simpson, inverse.simpson, S, pielou)
# 选择出hsct组和non组
alpha.diversity.hsct <- alpha.diversity[ , c(4, 5, 6, 10)]
alpha.diversity.non <- alpha.diversity[ , c(1, 2, 3, 7, 8, 9)]
# 使用Mann-Whitney U非参数检验比较各个α多样性指标
fivenum(alpha.diversity.hsct[5, ])
fivenum(alpha.diversity.non[5, ])
wilcox.test(alpha.diversity.non[5, ], alpha.diversity.hsct[5, ])
# α多样性指标绘制箱线图
boxplot(shannon.wiener ~ metadata$Group, data.frame(shannon.wiener, metadata$Group))
boxplot(pielou ~ metadata$Group, data.frame(pielou, metadata$Group))

# 计算β多样性指标
bray <- vegdist(t(species.count), method = "bray")
bray.matrix <- as.matrix(bray)
# β多样性绘制PCoA图
p <- beta_pcoa(bray.matrix, metadata, "Group")
p
beta_pcoa_stat(bray.matrix, metadata, "Group", result = "beta_pcoa_stat.txt")
pheatmap(bray.matrix, color = colorRampPalette(brewer.pal(7,"RdYlBu"))(100))
