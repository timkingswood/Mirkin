# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(readxl)
library(ggplot2)
library(dplyr)

# Read the DESeq2 result file
file_path <- "C:/Users/young/Desktop/RNAseq_Data/CD8+_E_SNA_vs_CD8+_N_HSNA_Differential_Expression.xlsx"
deseq2_results <- read_excel(file_path)

# Filter significant genes based on FDR adjusted p-value
significant_genes <- deseq2_results[deseq2_results$`FDR Adj p Value` < 0.05, ]

# Remove rows with NA values
significant_genes <- na.omit(significant_genes)

# Separate upregulated and downregulated genes
upregulated_genes <- significant_genes[significant_genes$`Log2 Fold Change` > 0, ]
downregulated_genes <- significant_genes[significant_genes$`Log2 Fold Change` < 0, ]

# Extract gene symbols for upregulated genes
up_gene_symbols <- upregulated_genes$`Gene Symbol`

# Extract gene symbols for downregulated genes
down_gene_symbols <- downregulated_genes$`Gene Symbol`

# Perform GO enrichment analysis for upregulated genes using Gene Symbol
go_enrichment_up <- enrichGO(gene = up_gene_symbols,
                             OrgDb = org.Mm.eg.db,
                             keyType = "SYMBOL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)

# Perform GO enrichment analysis for downregulated genes using Gene Symbol
go_enrichment_down <- enrichGO(gene = down_gene_symbols,
                               OrgDb = org.Mm.eg.db,
                               keyType = "SYMBOL",
                               ont = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2)

# Ensure Count and GeneRatio are numeric
go_enrichment_up@result$Count <- as.numeric(go_enrichment_up@result$Count)
go_enrichment_down@result$Count <- as.numeric(go_enrichment_down@result$Count)

# Calculate enrichment score as the ratio of observed gene count to the total number of genes in the category
go_enrichment_up@result$EnrichmentScore <- go_enrichment_up@result$Count / as.numeric(sub("/.*", "", go_enrichment_up@result$GeneRatio))
go_enrichment_down@result$EnrichmentScore <- -(go_enrichment_down@result$Count / as.numeric(sub("/.*", "", go_enrichment_down@result$GeneRatio)))

# Select the top 10 most significant GO terms for upregulated and downregulated genes
top_go_up <- go_enrichment_up@result %>% top_n(6, wt = -p.adjust)
top_go_down <- go_enrichment_down@result %>% top_n(6, wt = -p.adjust)

# Combine the results into one data frame
go_results <- rbind(
  data.frame(top_go_up, Regulation = "Upregulated"),
  data.frame(top_go_down, Regulation = "Downregulated")
)

library(ggplot2)
library(extrafont)  # Arial 폰트 사용을 위해 필요할 수도 있음
library(tools)  # toTitleCase 함수 사용을 위한 패키지

ggplot(go_results, aes(x = reorder(toTitleCase(Description), EnrichmentScore), 
                       y = EnrichmentScore * Count / max(Count))) +
  geom_col(aes(fill = -log10(p.adjust)), color = "black", alpha = 0.8) +  # p-value의 로그 변환
  geom_point(aes(y = EnrichmentScore, size = Count), 
             shape = 21, color = "black", fill = "white", stroke = 1.2) +
  coord_flip() + 
  labs(title = "GO Enrichment Analysis: Admix vs NSNA",
       x = "GO Term",
       y = "Adjusted Enrichment Score (Proportional to Gene Count)",
       size = "Gene Count",
       fill = "-log10(p-value)") +  # 범례를 -log10(p-value)로 수정
  theme_minimal(base_size = 14, base_family = "Arial") +  # Arial 폰트 적용
  scale_fill_gradient(low = "blue", high = "red", 
                      limits = c(0, max(-log10(go_results$p.adjust)))) +  # p-value 로그 값의 범위로 색상 설정
  scale_size(range = c(2, 6)) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold", family = "Arial"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold", family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "Arial"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold", family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial")
  )



