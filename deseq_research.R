library(DESeq2)
library("plyr")
library("dplyr")


VA6R3 <- read.delim("genecounts/VA6R3counts.txt", header = FALSE, skip = 2)

# НЕТ диф. экспр. angularis vs reflexopilosa

# angularis vs reflexopilosa DAI 2
sample_names <- c('VA2R1', 'VA2R2', 'VA2R3', 'VRP2R2', 'VRP2R3');

gc_ang_rp_DAI_2 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_ang_rp_DAI_2[name] <- sample_i$V2;
}

cd_ang_rp_DAI_2 <- read.csv('conditions/research/angularis_DAI_2.txt', sep=";", row.names=1, header=TRUE)

cd_ang_rp_DAI_2$Species <- as.factor(cd_ang_rp_DAI_2$Species)

dds_ang_rp_DAI_2 <- DESeqDataSetFromMatrix(countData = gc_ang_rp_DAI_2, colData = cd_ang_rp_DAI_2, design = ~ Species)

dds_ang_rp_DAI_2 <- DESeq(dds_ang_rp_DAI_2)

res_ang_rp_DAI_2 <- results(dds_ang_rp_DAI_2, name="Species_reflexopilosa_vs_angularis", alpha = 0.05, lfcThreshold = 1)

res_ang_rp_DAI_2_filt  <- res_ang_rp_DAI_2[which(res_ang_rp_DAI_2$padj > 0.05 & abs(res_ang_rp_DAI_2$log2FoldChange) < 1), ]

df_res_ang_rp_DAI_2_filt <- as.data.frame(res_ang_rp_DAI_2_filt)
df_res_ang_rp_DAI_2_filt$Gene <- rownames(df_res_ang_rp_DAI_2_filt)

# angularis vs reflexopilosa, stipulacea DAI 4
sample_names <- c('VA4R1', 'VA4R2', 'VRP4R2', 'VRP4R3');

gc_ang_rp_DAI_4 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_ang_rp_DAI_4[name] <- sample_i$V2;
}

cd_ang_rp_DAI_4 <- read.csv('conditions/research/angularis_DAI_4.txt', sep=";", row.names=1, header=TRUE)

cd_ang_rp_DAI_4$Species <- as.factor(cd_ang_rp_DAI_4$Species)

dds_ang_rp_DAI_4 <- DESeqDataSetFromMatrix(countData = gc_ang_rp_DAI_4, colData = cd_ang_rp_DAI_4, design = ~ Species)

dds_ang_rp_DAI_4 <- DESeq(dds_ang_rp_DAI_4)

res_ang_rp_DAI_4 <- results(dds_ang_rp_DAI_4, name="Species_reflexopilosa_vs_angularis", alpha = 0.05, lfcThreshold = 1)

res_ang_rp_DAI_4_filt  <- res_ang_rp_DAI_4[which(res_ang_rp_DAI_4$padj > 0.05 & abs(res_ang_rp_DAI_4$log2FoldChange) < 1), ]

df_res_ang_rp_DAI_4_filt <- as.data.frame(res_ang_rp_DAI_4_filt)
df_res_ang_rp_DAI_4_filt$Gene <- rownames(df_res_ang_rp_DAI_4_filt)

# angularis vs reflexopilosa, stipulacea DAI 6
sample_names <- c('VA6R1', 'VA6R2', 'VA6R3', 'VRP6R2', 'VRP6R3');

gc_ang_rp_DAI_6 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_ang_rp_DAI_6[name] <- sample_i$V2;
}

cd_ang_rp_DAI_6 <- read.csv('conditions/research/angularis_DAI_6.txt', sep=";", row.names=1, header=TRUE)

cd_ang_rp_DAI_6$Species <- as.factor(cd_ang_rp_DAI_6$Species)

dds_ang_rp_DAI_6 <- DESeqDataSetFromMatrix(countData = gc_ang_rp_DAI_6, colData = cd_ang_rp_DAI_6, design = ~ Species)

dds_ang_rp_DAI_6 <- DESeq(dds_ang_rp_DAI_6)

res_ang_rp_DAI_6 <- results(dds_ang_rp_DAI_6, name="Species_reflexopilosa_vs_angularis", alpha = 0.05, lfcThreshold = 1)

res_ang_rp_DAI_6_filt  <- res_ang_rp_DAI_6[which(res_ang_rp_DAI_6$padj > 0.05 & abs(res_ang_rp_DAI_6$log2FoldChange) < 1), ]

df_res_ang_rp_DAI_6_filt <- as.data.frame(res_ang_rp_DAI_6_filt)
df_res_ang_rp_DAI_6_filt$Gene <- rownames(df_res_ang_rp_DAI_6_filt)

# Объединение данных
ang_rp_DAI_2_4_6 <- inner_join(df_res_ang_rp_DAI_2_filt, df_res_ang_rp_DAI_4_filt, by=c("Gene" = "Gene"))
ang_rp_DAI_2_4_6 <- inner_join(ang_rp_DAI_2_4_6, df_res_ang_rp_DAI_6_filt, by=c("Gene" = "Gene"))


# ЕСТЬ диф. экспр. radiata vs angularis

# radiata vs angularis DAI 2
sample_names <- c('VR2R1', 'VR2R2', 'VR2R3', 'VA2R1', 'VA2R2', 'VA2R3');

gc_rad_ang_DAI_2 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_rad_ang_DAI_2[name] <- sample_i$V2;
}

cd_rad_ang_DAI_2 <- read.csv('conditions/research/radiata_DAI_2.txt', sep=";", row.names=1, header=TRUE)

cd_rad_ang_DAI_2$Species <- as.factor(cd_rad_ang_DAI_2$Species)
cd_rad_ang_DAI_2$Species <- relevel(cd_rad_ang_DAI_2$Species, 'radiata')

dds_rad_ang_DAI_2 <- DESeqDataSetFromMatrix(countData = gc_rad_ang_DAI_2, colData = cd_rad_ang_DAI_2, design = ~ Species)

dds_rad_ang_DAI_2 <- DESeq(dds_rad_ang_DAI_2)

res_rad_ang_DAI_2 <- results(dds_rad_ang_DAI_2, name="Species_angularis_vs_radiata", alpha = 0.05, lfcThreshold = 1)

res_rad_ang_DAI_2_filt  <- res_rad_ang_DAI_2[which(res_rad_ang_DAI_2$padj < 0.05 & abs(res_rad_ang_DAI_2$log2FoldChange) > 1), ]

df_res_rad_ang_DAI_2_filt <- as.data.frame(res_rad_ang_DAI_2_filt)
df_res_rad_ang_DAI_2_filt$Gene <- rownames(df_res_rad_ang_DAI_2_filt)

# radiata vs angularis DAI 4
sample_names <- c('VR4R1', 'VR4R2', 'VR4R3', 'VA4R1', 'VA4R2');

gc_rad_ang_DAI_4 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_rad_ang_DAI_4[name] <- sample_i$V2;
}

cd_rad_ang_DAI_4 <- read.csv('conditions/research/radiata_DAI_4.txt', sep=";", row.names=1, header=TRUE)

cd_rad_ang_DAI_4$Species <- as.factor(cd_rad_ang_DAI_4$Species)
cd_rad_ang_DAI_4$Species <- relevel(cd_rad_ang_DAI_4$Species, 'radiata')

dds_rad_ang_DAI_4 <- DESeqDataSetFromMatrix(countData = gc_rad_ang_DAI_4, colData = cd_rad_ang_DAI_4, design = ~ Species)

dds_rad_ang_DAI_4 <- DESeq(dds_rad_ang_DAI_4)

res_rad_ang_DAI_4 <- results(dds_rad_ang_DAI_4, name="Species_angularis_vs_radiata", alpha = 0.05, lfcThreshold = 1)

res_rad_ang_DAI_4_filt  <- res_rad_ang_DAI_4[which(res_rad_ang_DAI_4$padj < 0.05 & abs(res_rad_ang_DAI_4$log2FoldChange) > 1), ]

df_res_rad_ang_DAI_4_filt <- as.data.frame(res_rad_ang_DAI_4_filt)
df_res_rad_ang_DAI_4_filt$Gene <- rownames(df_res_rad_ang_DAI_4_filt)

# radiata vs angularis DAI 6
sample_names <- c('VR6R1', 'VR6R2', 'VR6R3', 'VA6R1', 'VA6R2', 'VA6R3');

gc_rad_ang_DAI_6 <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_rad_ang_DAI_6[name] <- sample_i$V2;
}

cd_rad_ang_DAI_6 <- read.csv('conditions/research/radiata_DAI_6.txt', sep=";", row.names=1, header=TRUE)

cd_rad_ang_DAI_6$Species <- as.factor(cd_rad_ang_DAI_6$Species)
cd_rad_ang_DAI_6$Species <- relevel(cd_rad_ang_DAI_6$Species, 'radiata')

dds_rad_ang_DAI_6 <- DESeqDataSetFromMatrix(countData = gc_rad_ang_DAI_6, colData = cd_rad_ang_DAI_6, design = ~ Species)

dds_rad_ang_DAI_6 <- DESeq(dds_rad_ang_DAI_6)

res_rad_ang_DAI_6 <- results(dds_rad_ang_DAI_6, name="Species_angularis_vs_radiata", alpha = 0.05, lfcThreshold = 1)

res_rad_ang_DAI_6_filt  <- res_rad_ang_DAI_6[which(res_rad_ang_DAI_6$padj < 0.05 & abs(res_rad_ang_DAI_6$log2FoldChange) > 1), ]

df_res_rad_ang_DAI_6_filt <- as.data.frame(res_rad_ang_DAI_6_filt)
df_res_rad_ang_DAI_6_filt$Gene <- rownames(df_res_rad_ang_DAI_6_filt)

# Объединение данных
rad_ang_DAI_2_4_6 <- rbind(df_res_rad_ang_DAI_2_filt, df_res_rad_ang_DAI_4_filt)
rad_ang_DAI_2_4_6 <- rbind(rad_ang_DAI_2_4_6, df_res_rad_ang_DAI_6_filt)
rad_ang_DAI_2_4_6 <- rad_ang_DAI_2_4_6[!duplicated(rad_ang_DAI_2_4_6[,c('Gene')]),]


stage_2_join <- inner_join(ang_rp_DAI_2_4_6, rad_ang_DAI_2_4_6, by=c("Gene" = "Gene"))

# stipulacea по дням
sample_names <- c('VS0R1', 'VS0R2', 'VS0R3', 'VS2R1', 'VS2R2', 'VS2R3', 'VS4R1', 'VS4R2', 'VS4R3', 'VS6R1', 'VS6R2', 'VS6R3');

gc_stipulacea <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_stipulacea[name] <- sample_i$V2;
}

cd_stipulacea <- read.csv('conditions/research/samplesStipulacea.txt', sep=";", row.names=1, header=TRUE)

cd_stipulacea$DAI <- as.factor(cd_stipulacea$DAI)

dds_stipulacea <- DESeqDataSetFromMatrix(countData = gc_stipulacea, colData = cd_stipulacea, design = ~ DAI)

dds_stipulacea <- DESeq(dds_stipulacea)

res_stipulacea.2 <- results(dds_stipulacea, name="DAI_2_vs_0", alpha = 0.05, lfcThreshold = 1)
res_stipulacea.4 <- results(dds_stipulacea, name="DAI_4_vs_0", alpha = 0.05, lfcThreshold = 1)
res_stipulacea.6 <- results(dds_stipulacea, name="DAI_6_vs_0", alpha = 0.05, lfcThreshold = 1)

res_stipulacea_filt.2  <- res_stipulacea.2[which(res_stipulacea.2$padj < 0.05 & abs(res_stipulacea.2 $log2FoldChange) > 1), ]
res_stipulacea_filt.4  <- res_stipulacea.4[which(res_stipulacea.4$padj < 0.05 & abs(res_stipulacea.4 $log2FoldChange) > 1), ]
res_stipulacea_filt.6  <- res_stipulacea.6[which(res_stipulacea.6$padj < 0.05 & abs(res_stipulacea.6 $log2FoldChange) > 1), ]

df_stipulacea_filt.2 <- as.data.frame(res_stipulacea_filt.2)
df_stipulacea_filt.4 <- as.data.frame(res_stipulacea_filt.4)
df_stipulacea_filt.6 <- as.data.frame(res_stipulacea_filt.6)

df_stipulacea_filt.2$Gene <- rownames(df_stipulacea_filt.2)
df_stipulacea_filt.4$Gene <- rownames(df_stipulacea_filt.4)
df_stipulacea_filt.6$Gene <- rownames(df_stipulacea_filt.6)

# Объединение данных
stipulacea_2_4_6 <- rbind(df_stipulacea_filt.2, df_stipulacea_filt.4)
stipulacea_2_4_6 <- rbind(stipulacea_2_4_6, df_stipulacea_filt.6)
stipulacea_2_4_6 <- stipulacea_2_4_6[!duplicated(stipulacea_2_4_6[,c('Gene')]),]

stage_3_join <- inner_join(stage_2_join, stipulacea_2_4_6, by=c("Gene" = "Gene"))
write.csv(stage_3_join, file = "stage_3_join.csv")

# Create GO universe

sample_names <- c('VA0R2', 'VA2R1', 'VA2R2', 'VA2R3', 'VA4R1', 'VA4R2', 'VA6R1', 'VA6R2', 'VA6R3',
                  'VR0R3', 'VR2R1', 'VR2R2', 'VR2R3', 'VR4R1', 'VR4R2', 'VR4R3', 'VR6R1', 'VR6R2', 'VR6R3',
                  'VRP0R2', 'VRP2R2', 'VRP2R3', 'VRP4R2', 'VRP4R3', 'VRP6R2', 'VRP6R3',
                  'VS0R1', 'VS0R2', 'VS0R3', 'VS2R1', 'VS2R2', 'VS2R3', 'VS4R1', 'VS4R2', 'VS4R3', 'VS6R1', 'VS6R2', 'VS6R3');

gc_univ <- data.frame(gene_id=VA6R3$V1, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_univ[name] <- sample_i$V2;
}

cd_univ <- read.csv('conditions/samples.txt', sep=";", row.names=1, header=TRUE)

cd_univ$DAI <- as.factor(cd_univ$DAI)
cd_univ$Species <- as.factor(cd_univ$Species)

dds_univ <- DESeqDataSetFromMatrix(countData = gc_univ, colData = cd_univ, design = ~ Species + DAI + Species:DAI)

dds_univ <- DESeq(dds_univ)

res_univ <- results(dds_univ)

df_univ <- as.data.frame(res_univ)

df_univ$Gene <- rownames(df_univ)
