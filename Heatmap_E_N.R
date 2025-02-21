# Load necessary libraries
library(readxl)
library(pheatmap)

# Set the new file path
file_path <- "C:/Users/young/Desktop/RNAseq_Data/CD8+_E_SNA_vs_CD8+_N_HSNA_Differential_Expression.xlsx"

# Read the data from the Excel file
deseq2_data <- read_excel(file_path)

# Filter the dataframe to include only significant genes
significant_genes <- deseq2_data[deseq2_data$Significant == "Yes", ]

# Extract the expression data and set row names to gene names
expression_data <- as.data.frame(significant_genes[, c("CD8T_N_1_S11", "CD8T_N_2_S16", "CD8T_N_3_S27", 
                                                       "CD8T_E_1_S10", "CD8T_E_2_S15", "CD8T_E_3_S26")])
rownames(expression_data) <- significant_genes$`Gene Symbol`

# Rename columns for labeling
colnames(expression_data) <- c("N_SNA_1", "N_SNA_2", "N_SNA_3", "E_SNA_1", "E_SNA_2", "E_SNA_3")

# Normalize the data using Z-score
expression_data_normalized <- t(scale(t(expression_data)))

# Define the list of genes to annotate
annotated_genes <- c("Sirpb1b", "Sirpb1a", "Ifitm3", "Tlr13", "Lag3", 
                     "Alox5ap", "Cd244a", "Serpina3g", "Pirb", "Cd7", 
                     "Ptpn6", "Trim39", "Atg16l2", "Rarg", "Socs1", 
                     "Ap1ar", "Ikbke", "Il18r1", "Clcf1", "Txk")

# Create a vector for row labels, showing only selected genes
row_labels <- ifelse(rownames(expression_data) %in% annotated_genes, rownames(expression_data), "")

# Save the heatmap with high DPI (600)
png("C:/Users/young/Desktop/heatmap_Admix_Nsna_optimized.png", width = 10, height = 8, units = "in", res = 600)

# Draw the heatmap
pheatmap(expression_data_normalized, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Disable column clustering
         show_rownames = TRUE,  # Show row names
         labels_row = row_labels,  # Display only selected gene names
         fontsize_row = 5,      # Keep original font size
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cellwidth = 40,         # Keep original cell width
         cellheight = 0.8,         # Increase the height of cells for more space
         treeheight_row = FALSE,  
         annotation_legend = TRUE,  # Show the legend for the annotations
         legend = FALSE)  # Disable the color legend

dev.off()  # Close the graphic device to save the image
