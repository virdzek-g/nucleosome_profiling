# create dataframe of nucleosome profiles, normalize for variation in ctDNA values across samples and find sites that are statisticly different between groups of samples
# input is output from Griffin (Doebley, AL., Ko, M., Liao, H. et al. A framework for clinical cancer subtyping from nucleosome profiling of cell-free DNA. Nat Commun 13, 7475 (2022). https://doi.org/10.1038/s41467-022-35076-w)


# read a single file to get a list of all profiled sites
file <- read.delim("GC_corrected.coverage.tsv" )

# loop across all profiled samples
list_files <- dir() 
data = data.frame(gene = file$site_name)
data_names <- data.frame()

for (i in list_files) {
    # read the data
    file_data <- read.delim(i) #
    file_data <- file_data[,c('central_coverage')]
    name <- sub('GC_corrected.coverage.tsv','',i)
    name <- sub('Griffin_results_','',name)   
    data_names <- rbind(data_names,name)
    data <- cbind(data,file_data)
    print(i)
     }
     
colnames(data)[1] <- 'sites'   
colnames(data)[2:ncol(data)] <- data_names[,1]

df <- data
rownames(df) <- df$sites
df <- df[,-1]

colnames(df) <- sub('_30X','',colnames(df))

# remove RBC samples
rbc <- df[,c(1:4)]
df <- df[,-c(1:4)]


# load ctDNA values
annot <- read.delim(./annotation_data.txt', row.names=1)

all.equal(rownames(annot), colnames(df))

##Quantily normalization
library(preprocessCore)
library(pheatmap)


# Normalize ctDNA values to range between 0 and 1
# Quantile normalization using ctDNA values as a covariate
df_normalized <- normalize.quantiles(as.matrix(df))

# Create a heatmap of the normalized data

colnames(df_normalized) <- colnames(df)
rownames(df_normalized) <- rownames(df)
pheatmap(df_normalized, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE)

### wilcox test between groups

library(preprocessCore)
library(pheatmap)
library(ggplot2)


group <- read.delim('./group.txt')

group <- group[match(colnames(df_normalized), group$sample),]

# check if the samples are correctly matched
if (!all(group$sample == colnames(df_normalized))) {
  stop("Sample names do not match between df and group.")
}


group$group <- as.factor(group$group)


wilcox.test_results <- apply(df_normalized, 1, function(tf) {
  wilcox.test(tf ~ annot$time_point)$p.value #t.test()
})


adjusted_p_values <- p.adjust(wilcox.test_results, method = "BH")


results <- data.frame(
  TranscriptionFactor = rownames(df_normalized),
  p.value = wilcox.test_results,
  adjusted.p.value = adjusted_p_values
)

significant_tfs <- results[results$adjusted.p.value < 0.05, ]
top_tfs <- head(significant_tfs[order(significant_tfs$adjusted.p.value), ], 20)
pheatmap(df_normalized[rownames(top_tfs), ],color = colorRampPalette(c("blue", "white", "red"))(100), scale = "row", cluster_rows = TRUE, cluster_cols = TRUE)
