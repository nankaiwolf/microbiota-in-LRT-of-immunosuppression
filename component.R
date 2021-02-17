# 载入运行库
library(amplicon)

# 读入细菌“门”级别的reads数据
phylum.count <- read.table("tax.phylum.count.tsv",
                           header = T, sep = "\t", row.names = 1)
# 计算比例
phylum.abund <- data.frame(abund01 = phylum.count$count_01 / sum(phylum.count$count_01) * 100,
                           abund02 = phylum.count$count_02 / sum(phylum.count$count_02) * 100,
                           abund03 = phylum.count$count_03 / sum(phylum.count$count_03) * 100,
                           abund04 = phylum.count$count_04 / sum(phylum.count$count_04) * 100,
                           abund05 = phylum.count$count_05 / sum(phylum.count$count_05) * 100,
                           abund06 = phylum.count$count_06 / sum(phylum.count$count_06) * 100,
                           abund07 = phylum.count$count_07 / sum(phylum.count$count_07) * 100,
                           abund08 = phylum.count$count_08 / sum(phylum.count$count_08) * 100,
                           abund09 = phylum.count$count_09 / sum(phylum.count$count_09) * 100,
                           abund10 = phylum.count$count_10 / sum(phylum.count$count_10) * 100,
                           All = rowSums(phylum.count) / sum(phylum.count) * 100,
                           row.names = rownames(phylum.count))
# 设置metadata的分组
metadata <- data.frame(Group = factor(c('non', 'non', 'non', 'hsct', 'hsct', 
                                        'hsct', 'non', 'non', 'non', 'hsct')), 
                       Series = c('01', '02', '03', '04', '05', 
                                  '06', '07', '08', '09', '10'),
                       Samples = c('BALF', 'BALF', 'BALF', 'BALF', 'BALF',
                                   'BALF', 'BALF', 'BALF', 'BALF', 'BALF'),
                       row.names = colnames(phylum.abund)[1:10])
# 对数据进行绘图
p <- tax_stackplot(phylum.abund, metadata, groupID = "Group", topN = 11)
p
p <- tax_stackplot(phylum.abund, metadata, groupID = "Group", style = "sample")
# 提取柔膜菌门的构成数据进行比较
