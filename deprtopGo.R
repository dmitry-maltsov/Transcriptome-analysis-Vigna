# topGO

# Добавляем GO аннотацию к генам
tmp_genes <- read.delim("genecounts/VA6R3_2ReadsPerGene.out.tab", header = FALSE, skip = 4);

go_vra <- read.csv(file = 'go.vra.csv', header = TRUE, skip = 8)

go_vra$gene_id <- paste0('gene:',go_vra$gene_id)

go_vra <- aggregate(x=go_vra, by=list(go_vra$gene_id), FUN=paste, collapse = ";")

sample_names <- c('VA0R2', 'VA2R1', 'VA2R2', 'VA2R3', 'VA4R1', 'VA4R2', 'VA6R1', 'VA6R2', 'VA6R3',
                  'VR0R3', 'VR2R1', 'VR2R2', 'VR2R3', 'VR4R1', 'VR4R2', 'VR4R3', 'VR6R1', 'VR6R2', 'VR6R3',
                  'VRP0R2', 'VRP2R2', 'VRP2R3', 'VRP4R2', 'VRP4R3', 'VRP6R2', 'VRP6R3',
                  'VS0R1', 'VS0R2', 'VS0R3', 'VS2R1', 'VS2R2', 'VS2R3', 'VS4R1', 'VS4R2', 'VS4R3', 'VS6R1', 'VS6R2', 'VS6R3');

gc_go <- data.frame(gene_id=tmp_genes$V1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_go[name] <- sample_i$V2;
}

gc_go_trimmed <- gc_go[gc_go$gene_id %in% go_vra$Group.1,]
go_vra_trimmed <- go_vra[go_vra$Group.1 %in% gc_go$gene_id,]


gc_go_trimmed <- data.frame(gc_go_trimmed, row.names=1);

for(name in sample_names){
  sample_i <- read.delim(paste0("genecounts_star/", name, "/ReadsPerGene.out.tab"), header = FALSE, skip = 4);
  gc_time[name] <- sample_i$V2;
}

go_ann <- data.frame(genes=rownames(gc_go_trimmed),
                     go=go_vra_trimmed$go)

cd_go <- read.csv('conditions/samples.txt', sep=";", row.names=1, header=TRUE)

cd_go$DAI <- as.factor(cd_go$DAI)
cd_go$Species <- as.factor(cd_go$Species)

dds_go <- DESeqDataSetFromMatrix(countData = gc_go_trimmed, colData = cd_go, design = ~ Species + DAI + Species:DAI )

dds_go <- DESeq(dds_go, test = "LRT", reduced = ~ Species + DAI)

res_go <- results(dds_go)


go_ann[] <- lapply(go_ann, as.character)

library("sanger-pathogens/deago")

addAnnotationsToDataSet(dds_go, good_ann)


library(topGO)