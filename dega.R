library(DESeq2)

# Разница между GeneCounts в STAR и в featureCounts
VA6R3_STAR <- read.delim("genecounts/VA6R3_2ReadsPerGene.out.tab", header = FALSE, skip = 4);
VS4R1_STAR <- read.delim("genecounts/VS4R1_2ReadsPerGene.out.tab", header = FALSE, skip = 4);
VR2R2_STAR <- read.delim("genecounts/VR2R2_2ReadsPerGene.out.tab", header = FALSE, skip = 4);


VA6R3 <- read.delim("genecounts/VA6R3counts.txt", header = FALSE, skip = 2);
VS4R1 <- read.delim("genecounts/VS4R1counts.txt", header = FALSE, skip = 2);
VR2R2 <- read.delim("genecounts/VR2R2counts.txt", header = FALSE, skip = 2);


difVA6R3 <- data.frame(name=VA6R3$V1, name_star=VA6R3_STAR$V1,
                              dif_names=ifelse(VA6R3_STAR$V1 == VA6R3$V1, FALSE, TRUE),
                       StarVA6R3=VA6R3_STAR$V2, VA6R3=VA6R3$V7, counts_dif=VA6R3_STAR$V2-VA6R3$V7);

difVS4R1 <- data.frame(name=VS4R1$V1, name_star=VS4R1_STAR$V1,
                              dif_names=ifelse(VS4R1_STAR$V1 == VS4R1$V1, FALSE, TRUE),
                       StarVS4R1=VS4R1_STAR$V2, VS4R1=VS4R1$V7, counts_dif=VS4R1_STAR$V2-VS4R1$V7);

difVR2R2 <- data.frame(name=VR2R2$V1, name_star=VR2R2_STAR$V1,
                              dif_names=ifelse(VR2R2_STAR$V1 == VR2R2$V1, FALSE, TRUE),
                       StarVR2R2=VR2R2_STAR$V2, VR2R2=VR2R2$V7, counts_dif=VR2R2_STAR$V2-VR2R2$V7);


# Создание матриц GeneCounts и ColData для различных экспериментов

# angularis
sample_names <- c('VA0R2', 'VA2R1', 'VA2R2', 'VA2R3', 'VA4R1', 'VA4R2', 'VA6R1', 'VA6R2', 'VA6R3');

gc_angularis <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_angularis[name] <- sample_i$V2;
}

cd_angularis <- read.csv('conditions/samplesAngularis.txt', sep=";", row.names=1, header=TRUE)

cd_angularis$DAI <- as.factor(cd_angularis$DAI)

dds_angularis <- DESeqDataSetFromMatrix(countData = gc_angularis, colData = cd_angularis, design = ~ DAI)

# radiata
sample_names <- c('VR0R3', 'VR2R1', 'VR2R2', 'VR2R3', 'VR4R1', 'VR4R2', 'VR4R3', 'VR6R1', 'VR6R2', 'VR6R3');

gc_radiata <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_radiata[name] <- sample_i$V2;
}

cd_radiata <- read.csv('conditions/samplesRadiata.txt', sep=";", row.names=1, header=TRUE)

cd_radiata$DAI <- as.factor(cd_radiata$DAI)

dds_radiata <- DESeqDataSetFromMatrix(countData = gc_radiata, colData = cd_radiata, design = ~ DAI)

# reflexopilosa
sample_names <- c('VRP0R2', 'VRP2R2', 'VRP2R3', 'VRP4R2', 'VRP4R3', 'VRP6R2', 'VRP6R3');

gc_reflexopilosa <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_reflexopilosa[name] <- sample_i$V2;
}

cd_reflexopilosa <- read.csv('conditions/samplesReflexopilosa.txt', sep=";", row.names=1, header=TRUE)

cd_reflexopilosa$DAI <- as.factor(cd_reflexopilosa$DAI)

dds_reflexopilosa <- DESeqDataSetFromMatrix(countData = gc_reflexopilosa, colData = cd_reflexopilosa, design = ~ DAI)

# stipulacea
sample_names <- c('VS0R1', 'VS0R2', 'VS0R3', 'VS2R1', 'VS2R2', 'VS2R3', 'VS4R1', 'VS4R2', 'VS4R3', 'VS6R1', 'VS6R2', 'VS6R3');

gc_stipulacea <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_stipulacea[name] <- sample_i$V2;
}

cd_stipulacea <- read.csv('conditions/samplesStipulacea.txt', sep=";", row.names=1, header=TRUE)

cd_stipulacea$DAI <- as.factor(cd_stipulacea$DAI)

dds_stipulacea <- DESeqDataSetFromMatrix(countData = gc_stipulacea, colData = cd_stipulacea, design = ~ DAI)

# radiata vs angularis, reflexopilosa, stipulacea DAI 0
sample_names <- c('VR0R3', 'VA0R2', 'VRP0R2', 'VS0R1', 'VS0R2', 'VS0R3');

gc_radiata_DAI_0 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_radiata_DAI_0[name] <- sample_i$V2;
}

cd_radiata_DAI_0 <- read.csv('conditions/radiata_DAI_0.txt', sep=";", row.names=1, header=TRUE)

cd_radiata_DAI_0$Species <- as.factor(cd_radiata_DAI_0$Species)
cd_radiata_DAI_0$Species <- relevel(cd_radiata_DAI_0$Species, 'radiata')

dds_radiata_DAI_0 <- DESeqDataSetFromMatrix(countData = gc_radiata_DAI_0, colData = cd_radiata_DAI_0, design = ~ Species)

# radiata vs angularis, reflexopilosa, stipulacea DAI 2
sample_names <- c('VR2R1', 'VR2R2', 'VR2R3', 'VA2R1', 'VA2R2', 'VA2R3', 'VRP2R2', 'VRP2R3', 'VS2R1', 'VS2R2', 'VS2R3');

gc_radiata_DAI_2 <- data.frame(gene_id=VA6R3$V1, row.names=1);

dds_stipulacea <- DESeqDataSetFromMatrix(countData = gc_stipulacea, colData = cd_stipulacea, design = ~ DAI)
for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_radiata_DAI_2[name] <- sample_i$V2;
}

cd_radiata_DAI_2 <- read.csv('conditions/radiata_DAI_2.txt', sep=";", row.names=1, header=TRUE)

cd_radiata_DAI_2$Species <- as.factor(cd_radiata_DAI_2$Species)
cd_radiata_DAI_2$Species <- relevel(cd_radiata_DAI_2$Species, 'radiata')

dds_radiata_DAI_2 <- DESeqDataSetFromMatrix(countData = gc_radiata_DAI_2, colData = cd_radiata_DAI_2, design = ~ Species)

# radiata vs angularis, reflexopilosa, stipulacea DAI 4
sample_names <- c('VR4R1', 'VR4R2', 'VR4R3', 'VA4R1', 'VA4R2', 'VRP4R2', 'VRP4R3', 'VS4R1', 'VS4R2', 'VS4R3');

gc_radiata_DAI_4 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_radiata_DAI_4[name] <- sample_i$V2;
}

cd_radiata_DAI_4 <- read.csv('conditions/radiata_DAI_4.txt', sep=";", row.names=1, header=TRUE)

cd_radiata_DAI_4$Species <- as.factor(cd_radiata_DAI_4$Species)
cd_radiata_DAI_4$Species <- relevel(cd_radiata_DAI_4$Species, 'radiata')

dds_radiata_DAI_4 <- DESeqDataSetFromMatrix(countData = gc_radiata_DAI_4, colData = cd_radiata_DAI_4, design = ~ Species)

# radiata vs angularis, reflexopilosa, stipulacea DAI 6
sample_names <- c('VR6R1', 'VR6R2', 'VR6R3', 'VA6R1', 'VA6R2', 'VA6R3', 'VRP6R2', 'VRP6R3', 'VS6R1', 'VS6R2', 'VS6R3');

gc_radiata_DAI_6 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_radiata_DAI_6[name] <- sample_i$V2;
}

cd_radiata_DAI_6 <- read.csv('conditions/radiata_DAI_6.txt', sep=";", row.names=1, header=TRUE)

cd_radiata_DAI_6$Species <- as.factor(cd_radiata_DAI_6$Species)
cd_radiata_DAI_6$Species <- relevel(cd_radiata_DAI_6$Species, 'radiata')

dds_radiata_DAI_6 <- DESeqDataSetFromMatrix(countData = gc_radiata_DAI_6, colData = cd_radiata_DAI_6, design = ~ Species)

# angularis vs reflexopilosa, stipulacea DAI 0
sample_names <- c('VA0R2', 'VRP0R2', 'VS0R1', 'VS0R2', 'VS0R3');

gc_angularis_DAI_0 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_angularis_DAI_0[name] <- sample_i$V2;
}

cd_angularis_DAI_0 <- read.csv('conditions/angularis_DAI_0.txt', sep=";", row.names=1, header=TRUE)

cd_angularis_DAI_0$Species <- as.factor(cd_angularis_DAI_0$Species)

dds_angularis_DAI_0 <- DESeqDataSetFromMatrix(countData = gc_angularis_DAI_0, colData = cd_angularis_DAI_0, design = ~ Species)

# angularis vs reflexopilosa, stipulacea DAI 2
sample_names <- c('VA2R1', 'VA2R2', 'VA2R3', 'VRP2R2', 'VRP2R3', 'VS2R1', 'VS2R2', 'VS2R3');

gc_angularis_DAI_2 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_angularis_DAI_2[name] <- sample_i$V2;
}

cd_angularis_DAI_2 <- read.csv('conditions/angularis_DAI_2.txt', sep=";", row.names=1, header=TRUE)

cd_angularis_DAI_2$Species <- as.factor(cd_angularis_DAI_2$Species)

dds_angularis_DAI_2 <- DESeqDataSetFromMatrix(countData = gc_angularis_DAI_2, colData = cd_angularis_DAI_2, design = ~ Species)

# angularis vs reflexopilosa, stipulacea DAI 4
sample_names <- c('VA4R1', 'VA4R2', 'VRP4R2', 'VRP4R3', 'VS4R1', 'VS4R2', 'VS4R3');

gc_angularis_DAI_4 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_angularis_DAI_4[name] <- sample_i$V2;
}

cd_angularis_DAI_4 <- read.csv('conditions/angularis_DAI_4.txt', sep=";", row.names=1, header=TRUE)

cd_angularis_DAI_4$Species <- as.factor(cd_angularis_DAI_4$Species)

dds_angularis_DAI_4 <- DESeqDataSetFromMatrix(countData = gc_angularis_DAI_4, colData = cd_angularis_DAI_4, design = ~ Species)

# angularis vs reflexopilosa, stipulacea DAI 6
sample_names <- c('VA6R1', 'VA6R2', 'VA6R3', 'VRP6R2', 'VRP6R3', 'VS6R1', 'VS6R2', 'VS6R3');

gc_angularis_DAI_6 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_angularis_DAI_6[name] <- sample_i$V2;
}

cd_angularis_DAI_6 <- read.csv('conditions/angularis_DAI_6.txt', sep=";", row.names=1, header=TRUE)

cd_angularis_DAI_6$Species <- as.factor(cd_angularis_DAI_6$Species)

dds_angularis_DAI_6 <- DESeqDataSetFromMatrix(countData = gc_angularis_DAI_6, colData = cd_angularis_DAI_6, design = ~ Species)

# stipulacea vs reflexopilosa DAI 0
sample_names <- c('VS0R1', 'VS0R2', 'VS0R3', 'VRP0R2');

gc_stipulacea_DAI_0 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_stipulacea_DAI_0[name] <- sample_i$V2;
}

cd_stipulacea_DAI_0 <- read.csv('conditions/stipulacea_DAI_0.txt', sep=";", row.names=1, header=TRUE)

cd_stipulacea_DAI_0$Species <- as.factor(cd_stipulacea_DAI_0$Species)
cd_stipulacea_DAI_0$Species <- relevel(cd_stipulacea_DAI_0$Species, 'stipulacea')

dds_stipulacea_DAI_0 <- DESeqDataSetFromMatrix(countData = gc_stipulacea_DAI_0, colData = cd_stipulacea_DAI_0, design = ~ Species)

# stipulacea vs reflexopilosa DAI 2
sample_names <- c('VS2R1', 'VS2R2', 'VS2R3', 'VRP2R2', 'VRP2R3');

gc_stipulacea_DAI_2 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_stipulacea_DAI_2[name] <- sample_i$V2;
}

cd_stipulacea_DAI_2 <- read.csv('conditions/stipulacea_DAI_2.txt', sep=";", row.names=1, header=TRUE)

cd_stipulacea_DAI_2$Species <- as.factor(cd_stipulacea_DAI_2$Species)
cd_stipulacea_DAI_2$Species <- relevel(cd_stipulacea_DAI_2$Species, 'stipulacea')

dds_stipulacea_DAI_2 <- DESeqDataSetFromMatrix(countData = gc_stipulacea_DAI_2, colData = cd_stipulacea_DAI_2, design = ~ Species)

# stipulacea vs reflexopilosa DAI 4
sample_names <- c('VS4R1', 'VS4R2', 'VS4R3', 'VRP4R2', 'VRP4R3');

gc_stipulacea_DAI_4 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_stipulacea_DAI_4[name] <- sample_i$V2;
}

cd_stipulacea_DAI_4 <- read.csv('conditions/stipulacea_DAI_4.txt', sep=";", row.names=1, header=TRUE)

cd_stipulacea_DAI_4$Species <- as.factor(cd_stipulacea_DAI_4$Species)
cd_stipulacea_DAI_4$Species <- relevel(cd_stipulacea_DAI_4$Species, 'stipulacea')

dds_stipulacea_DAI_4 <- DESeqDataSetFromMatrix(countData = gc_stipulacea_DAI_4, colData = cd_stipulacea_DAI_4, design = ~ Species)

# stipulacea vs reflexopilosa DAI 6
sample_names <- c('VS6R1', 'VS6R2', 'VS6R3', 'VRP6R2', 'VRP6R3');

gc_stipulacea_DAI_6 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_stipulacea_DAI_6[name] <- sample_i$V2;
}

cd_stipulacea_DAI_6 <- read.csv('conditions/stipulacea_DAI_6.txt', sep=";", row.names=1, header=TRUE)

cd_stipulacea_DAI_6$Species <- as.factor(cd_stipulacea_DAI_6$Species)
cd_stipulacea_DAI_6$Species <- relevel(cd_stipulacea_DAI_6$Species, 'stipulacea')

dds_stipulacea_DAI_6 <- DESeqDataSetFromMatrix(countData = gc_stipulacea_DAI_6, colData = cd_stipulacea_DAI_6, design = ~ Species)


# DESeq

dds_angularis <- DESeq(dds_angularis)
dds_radiata <- DESeq(dds_radiata)
dds_reflexopilosa <- DESeq(dds_reflexopilosa)
dds_stipulacea <- DESeq(dds_stipulacea)

dds_radiata_DAI_0 <- DESeq(dds_radiata_DAI_0)
dds_radiata_DAI_2 <- DESeq(dds_radiata_DAI_2)
dds_radiata_DAI_4 <- DESeq(dds_radiata_DAI_4)
dds_radiata_DAI_6 <- DESeq(dds_radiata_DAI_6)

dds_angularis_DAI_0 <- DESeq(dds_angularis_DAI_0)
dds_angularis_DAI_2 <- DESeq(dds_angularis_DAI_2)
dds_angularis_DAI_4 <- DESeq(dds_angularis_DAI_4)
dds_angularis_DAI_6 <- DESeq(dds_angularis_DAI_6)

dds_stipulacea_DAI_0 <- DESeq(dds_stipulacea_DAI_0)
dds_stipulacea_DAI_2 <- DESeq(dds_stipulacea_DAI_2)
dds_stipulacea_DAI_4 <- DESeq(dds_stipulacea_DAI_4)
dds_stipulacea_DAI_6 <- DESeq(dds_stipulacea_DAI_6)

print(resultsNames(dds_angularis))
print(resultsNames(dds_radiata))
print(resultsNames(dds_reflexopilosa))
print(resultsNames(dds_stipulacea))

print(resultsNames(dds_radiata_DAI_0))
print(resultsNames(dds_radiata_DAI_2))
print(resultsNames(dds_radiata_DAI_4))
print(resultsNames(dds_radiata_DAI_6))

print(resultsNames(dds_angularis_DAI_0))
print(resultsNames(dds_angularis_DAI_2))
print(resultsNames(dds_angularis_DAI_4))
print(resultsNames(dds_angularis_DAI_6))

print(resultsNames(dds_stipulacea_DAI_0))
print(resultsNames(dds_stipulacea_DAI_2))
print(resultsNames(dds_stipulacea_DAI_4))
print(resultsNames(dds_stipulacea_DAI_6))

res_angularis <- results(dds_angularis)
res_radiata <- results(dds_radiata)
res_reflexopilosa <- results(dds_reflexopilosa)
res_stipulacea <- results(dds_stipulacea)

res_radiata_DAI_0 <- results(dds_radiata_DAI_0)
res_radiata_DAI_2 <- results(dds_radiata_DAI_2)
res_radiata_DAI_4 <- results(dds_radiata_DAI_4)
res_radiata_DAI_6 <- results(dds_radiata_DAI_6)

res_angularis_DAI_0 <- results(dds_angularis_DAI_0)
res_angularis_DAI_2 <- results(dds_angularis_DAI_2)
res_angularis_DAI_4 <- results(dds_angularis_DAI_4)
res_angularis_DAI_6 <- results(dds_angularis_DAI_6)

res_stipulacea_DAI_0 <- results(dds_stipulacea_DAI_0)
res_stipulacea_DAI_2 <- results(dds_stipulacea_DAI_2)
res_stipulacea_DAI_4 <- results(dds_stipulacea_DAI_4)
res_stipulacea_DAI_6 <- results(dds_stipulacea_DAI_6)

# Визулизация результатов

# Изменение экспрессии генов с течением времени
sample_names <- c('VA0R2', 'VA2R1', 'VA2R2', 'VA2R3', 'VA4R1', 'VA4R2', 'VA6R1', 'VA6R2', 'VA6R3',
                  'VR0R3', 'VR2R1', 'VR2R2', 'VR2R3', 'VR4R1', 'VR4R2', 'VR4R3', 'VR6R1', 'VR6R2', 'VR6R3',
                  'VRP0R2', 'VRP2R2', 'VRP2R3', 'VRP4R2', 'VRP4R3', 'VRP6R2', 'VRP6R3',
                  'VS0R1', 'VS0R2', 'VS0R3', 'VS2R1', 'VS2R2', 'VS2R3', 'VS4R1', 'VS4R2', 'VS4R3', 'VS6R1', 'VS6R2', 'VS6R3');

gc_time <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_time[name] <- sample_i$V2;
}

cd_time <- read.csv('conditions/samples.txt', sep=";", row.names=1, header=TRUE)

cd_time$DAI <- as.factor(cd_time$DAI)
cd_time$Species <- as.factor(cd_time$Species)

dds_time <- DESeqDataSetFromMatrix(countData = gc_time, colData = cd_time, design = ~ Species + DAI + Species:DAI )

dds_time <- DESeq(dds_time, test = "LRT", reduced = ~ Species + DAI)

res_time <- results(dds_time)

fiss <- plotCounts(dds_time, which.min(res_time$padj), 
                   intgroup = c("DAI","Species"), returnData = TRUE)

fiss$DAI <- as.numeric(as.character(fiss$DAI))

pdf(file="time_plot.pdf")
ggplot(fiss,
       aes(x = DAI, y = count, color = Species, group = Species)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()
dev.off()

# Volcano Plot
library(EnhancedVolcano)

pdf(file="volcano_reflexopilosa_DAI_2_vs_0.pdf")

EnhancedVolcano(res_reflexopilosa,
                lab = rownames(res_reflexopilosa),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Reflexopilosa DAI 2 versus DAI 0')

dev.off()

pdf(file="volcano_stipulacea_DAI_2_vs_0.pdf")

EnhancedVolcano(res_reflexopilosa,
                lab = rownames(res_stipulacea),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Stipulacea DAI 2 versus DAI 0')

dev.off()

pdf(file="volcano_angularis_vs_radiata_DAI_2.pdf")

EnhancedVolcano(res_radiata_DAI_2,
                lab = rownames(res_radiata_DAI_2),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Stipulacea DAI 2 versus DAI 0')

dev.off()

res_radiata_DAI_0_df <- as.data.frame(res_radiata_DAI_0)
res_radiata_DAI_2_df <- as.data.frame(res_radiata_DAI_2)
res_radiata_DAI_4_df <- as.data.frame(res_radiata_DAI_4)
res_radiata_DAI_6_df <- as.data.frame(res_radiata_DAI_6)

res_angularis_DAI_0_df <- as.data.frame(res_angularis_DAI_0)
res_angularis_DAI_2_df <- as.data.frame(res_angularis_DAI_2)
res_angularis_DAI_4_df <- as.data.frame(res_angularis_DAI_4)
res_angularis_DAI_6_df <- as.data.frame(res_angularis_DAI_6)

res_stipulacea_DAI_0_df <- as.data.frame(res_stipulacea_DAI_0)
res_stipulacea_DAI_2_df <- as.data.frame(res_stipulacea_DAI_2)
res_stipulacea_DAI_4_df <- as.data.frame(res_stipulacea_DAI_4)
res_stipulacea_DAI_6_df <- as.data.frame(res_stipulacea_DAI_6)

