
#文章复现练习 《Using clusterProfiler to characterize multiomics data》
#来源：S Xu#, E Hu#, Y Cai#, Z Xie#, X Luo#, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu*. Using clusterProfiler to characterise Multi-Omics Data. Nature Protocols. 2024, https://doi.org/10.1038/s41596-024-01020-z. (IF: 13.1)
#完成时间：2026.02.05

#packages
pkgs <- c("aplot", "CelliD", "clusterProfiler", "DESeq2", "enrichplot", "ggfun", "ggplot2", "ggrepel", "ggsc", "MicrobiotaProcess", "Seurat") 
install.packages("BiocManager") 
BiocManager::install(pkgs)

#Download group information
url <- "https://yulab-smu.top/clusterProfiler_protocol/examples" 
IBD_files <- c("mg.meta.csv", "mg.expr.csv", "metabolism_meta.csv", "metabolism_expr.csv") 
for (f in IBD_files) { download.file(file.path(url, "IBD_2_subtypes", f), destfile = f) }

#raw count and the experimental group information
count_files <- c("counts.txt", "group_info.txt")  
url <- "https://yulab-smu.top/clusterProfiler_protocol/examples/Phyllostachys_heterocycla"  
for (f in count_files) { download.file(file.path(url, f), destfile = f) }

#Download PE annotation data
annot_files <- c("regulation_from_motif_CE_Phe.txt", "Phe_TF_list.txt", "Phe_GO_annotation.txt" )  
for (f in annot_files) { download.file(file.path(url, "annot_data", f), destfile = f) }

#Single-cell transcriptomic data and its annotation for cell type identification
#Download the single-cell dataset
url <- "https://yulab-smu.top/clusterProfiler_protocol/examples/single_cell/" 
dir.create("hg19") 
pbmc_files <- c("barcodes.tsv", "genes.tsv", "matrix.mtx") 
for (f in pbmc_files) {
  download.file( 
    file.path(url, "filtered_gene_bc_matrices/hg19/", f), 
    method = "auto", 
    destfile = file.path("hg19", f) 
  ) 
}

#Download gene-cell type annotation data
download.file( file.path(url, "cell_marker_db/c8.all.v2023.1.Hs.symbols.gmt"), 
               method = "auto", 
               destfile = "c8.all.v2023.1.Hs.symbols.gmt" 
               )

###Procedure 1: metabolomics and metagenomics functional enrichment analysis
#Set up the environment and data objects
library(MicrobiotaProcess) 
BiocManager::install("GO.db")
#BiocManager::install("clusterProfiler", dependencies = TRUE) #补充：所有安装依赖包
library(clusterProfiler) 
library(ggplot2) 
library(enrichplot)
meta_mg <- read.csv("mg.meta.csv") 
metagenome <- read.csv("mg.expr.csv", row.names = 1, check.name = FALSE) #用文件的第 1 列作为数据框的行名，关闭 R 对列名的自动修正，保留原始列名


#Metagenomic data differential analysis 
 #Define a function to perform differential analysis:
DA <- function(expr, meta, abundance = "Abundance", group_colname = "Diagnosis",  #function自定义函数的固定关键字
               force = TRUE, 
               relative = FALSE, 
               subset_group, 
               diff_group, filter.p = "pvalue", ...) {
  sign_group_colname <- paste0("Sign_", group_colname)
  mpse <- MPSE(expr) 
  mpse <- mpse |> left_join(meta, by = "Sample") 
  mpse |>   #|>链式传导
    dplyr::filter(!!as.symbol(group_colname) %in% subset_group) |>  #%in% 左侧是否存在于右侧
    mp_diff_analysis( 
      .abundance = !!as.symbol(abundance), 
      .group = !!as.symbol(group_colname), 
      force = force, 
      relative = relative, 
      filter.p = filter.p, 
      ...  #省略参数
    ) |>
    mp_extract_feature() |> 
    dplyr::filter(!!as.symbol(sign_group_colname) == diff_group) |> 
    dplyr::pull(OTU) |> suppressMessages() 
}

#Identify differentially abundant microbial genes:
groups <- c(CD = "CD", UC = "UC")                             #定义了需要分析的疾病分组（CD：克罗恩病；UC：溃疡性结肠炎）
de_gene <- lapply(groups, function(x) {                       #lapply遍历 groups 中的每个分组，批量执行差异分析
  DA(expr = metagenome, 
     meta = meta_mg, 
     subset_group = c(x, "Control"), diff_group = x) 
  })

#Functional analysis of differential microbial genes
library(clusterProfiler)
gene_enrich_result <- compareCluster(geneClusters = de_gene,   #多组基因集富集分析 输入差异基因列表
                                     fun = "enrichKEGG",       #指定使用 KEGG 通路富集算法
                                     organism = "ko")          #指定数据库 KO
options(timeout = 300)

#Visualize the functional enrichment result of differential genes
dotplot(gene_enrich_result, facet = "intersect", showCategory = 10,       #dotpolt绘图函数 facet+split 分离、前10 （allow the separation of common and unique biological pathways that are enriched in different IBD subtypes）
        split = "intersect", label_format = 60) +                         #“label_format”参数允许用户指定路径标识符的显示长度，如果长度超过指定的限制，则自动换行。
  ggtitle("Functional enrichment of intestinal genes") +                  #title
  theme(plot.title = element_text(hjust = 1))                             #修改非图标样式、图标标题、文字、水平右对齐


#################################################################
##Metabolomic data differential analysis
#Import metabolomic data
meta_mb <- read.csv("metabolism_meta.csv") 
metabolism <- read.csv("metabolism_expr.csv", 
                       row.names = 1, 
                       check.name = FALSE)

#Identify differentially abundance metabolites
groups <- c(CD = "CD", UC = "UC") 
de_cpd <- lapply(groups, function(x) { 
  DA(expr = metabolism, 
     meta = meta_mb, 
     subset_group = c(x, "Control"), 
     diff_group = x) 
  })

##Functional analysis of differential metabolites
#Perform enrichment analysis of differential metabolites
cpd_enrich_result <- compareCluster(geneClusters = de_cpd, 
                                    fun = "enrichKEGG", 
                                    organism = "cpd")

#Visualize the functional enrichment results of differential metabolites
dotplot(cpd_enrich_result, facet = "intersect", showCategory = 10, 
        split = "intersect", label_format = 60) + 
  ggtitle("Functional enrichment of chemical compounds") + 
  theme(plot.title = element_text(hjust = 1))



######################
#transcription factor analysis pertaining to cold tolerance in PE
#PREPARE expression matrix
library(DESeq2) 
library(clusterProfiler) 
library(ggplot2) 
library(enrichplot) 
library(aplot) 
library(ggfun)
counts <- read.delim("counts.txt")              #read.delim() 读取.txt 表格文件 、read.csv() 读取.csv文件
group_info <- read.delim("group_info.txt")        #Vars赋值
group_info$group <- factor(group_info$group, levels = c("0h", "2h", "24h", "168h"))   #数据列定义 factor为分类数据 level分类依据

#Calculate logarithmic fold changes (log2FC) using DESeq2

#Create a DESeqDataSet object and perform default differential expression analysis
count_dds <- DESeqDataSetFromMatrix(countData = counts,            #打包DESeq2 专用的分析对象DESeqDataSetFromMatrix
                                    colData = group_info, 
                                    design = ~ group)              #因变量（空）~（受什么影响）自变量1+自变量2...
count_dds <- DESeq(count_dds)                                      #运行 DESeq2 的差异分析流程

#Extract log2FC, and sort genes based on the log2FC values
#从 DESeq2 的差异分析结果中，批量提取不同时间点（168h、24h、2h）与 0h 相比的基因表达变化倍数（log2FC），并按变化幅度从大到小排序。
#Define groups
time_points <- setNames(object = c("168h", "24h", "2h"),                   #object比较的对象
                        nm = c("168h_vs_0h", "24h_vs_0h", "2h_vs_0h"))     #比较名称

#Extract log2FC values
all_result <- lapply(time_points, function(time_point) {                  #lapply(time_points, function(time_point) {...}) 遍历+执行function  该命令总成all_result列表
  result <- results(count_dds, tidy = TRUE,                               #results从对象（count_dds）中提取差异结果 tidy整洁表格
                    contrast = c("group", time_point, "0h"))              #明确提取规则group列中time_point和0h对比
  setNames(object = result$log2FoldChange, nm = result$row) |>            #设置rownames=result$row为log2FoldChange
    sort(decreasing = TRUE)                                               #从大到小排序
  })

##Transcription factor enrichment analysis to identify perturbed transcription factors at different timepoints
#Prepare transcription factor annotation data
tf_db <- read.delim( "regulation_from_motif_CE_Phe.txt",                                                               #TFdatabase
                     header = FALSE, colClasses = c("character", "NULL", "character", rep("NULL", 4))                  #首行header非列名，仅读取1、3列，其余NULL rep重复函数
                     ) |> setNames(c("TF", "targetGene"))                                                              #命名

#Use compareCluster to perform batch GSEA for identifying perturbed transcription factors from each timepoint
perturbed_TF_result <- compareCluster(all_result, fun = "GSEA",                                                        #compareCluster命令批量富集分析
                                      pvalueCutoff = .05,                                                              #fun指定GESA算法；pvalue P值0.05；TERM2GENE自定义富集分析用tf_db注释库；seed设置随机种子，可重复性
                                      TERM2GENE = tf_db, 
                                      seed = 1234)

#Visualize the enrichment result using dotplot
perturbed_TF_plot <- dotplot(perturbed_TF_result,                                  #可视化点图dotplot;aes定义点的形状I使用（22）预设值：填充正方形
                             showCategory = 25) + 
  aes(shape = I(22)) + 
  coord_flip() +                                                                   #默认横向柱转为纵向柱
  theme_minimal() +                                                                #简化主题 no x轴刻度及标签、lab标题NULL
  theme_noxaxis() + 
  xlab(NULL) + 
  set_enrichplot_color(  c("#6C8FAD", "#84ADA7", "#C7B398"),                       #自定义颜色函数，.fun指定渐变方案
                         .fun = ggplot2::scale_fill_gradientn 
                         )

##Predict the biological functions possibly regulated by the perturbed transcription factors

#Construct a named list of transcription factor (TF) target genes and extract a subset of the most significant TFs that are related to cold responses
tf_id <- unique(get_plot_data(perturbed_TF_plot, "ID")[, 1])               #提取点图数据中的ID，unique单一，避免重复
tf_genes <- split(tf_db$targetGene, tf_db$TF)[tf_id]                       #split分组tf_db的targetGene、TF列，仅保留上步[tf_id]

#Import the GO annotations of PE genes
go_db <- read.delim(file = "Phe_GO_annotation.txt")

#Characterize TF functions based on GO enrichment analysis of TF target genes
TF_GO_result <- compareCluster(tf_genes, fun = "enricher",                 #fun指定enricher超几何检验ORA进行富集
                               TERM2GENE = go_db[, c(2, 1)],               #定义变量TERM2GENE、TERM2NAME，从go_db提取2、1列//2、3列
                               TERM2NAME = go_db[, c(2, 3)])

#Visualize TF function enrichment result
TF_GO_plot <- dotplot(TF_GO_result, by = "count",                         #count点的大小by数量决定；标签格式40自动换行
                      showCategory = 3, 
                      label_format = 40) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1,                  #旋转 X 轴标签 30 度，并调整对齐方式和大小
                                   angle = 30, size = 10)) + 
  xlab(NULL) + ggtitle(NULL)                                              #x轴和主标题去除

##Visualization of transcription factor family annotation information

#import TF family annotation information
tf_family <- read.delim("Phe_TF_list.txt", row.names = 1)                 #文件第一列作为行名

#Prepare TF family data for visualization
family_data <- subset(tf_family, Gene_ID %in% tf_id)                      #%in%左是否在右里
family_data <- family_data[order(family_data$Family),]                    #按照 Family列对数据进行排序
family_data$Gene_ID <- factor(family_data$Gene_ID,                        #factor因子类型数据转化；level排列顺序
                              levels = family_data$Gene_ID)

#Visualize the TF family information
tf_family_plot <- ggplot(data = family_data,                                  #ggplot初始化绘图，fill = Family：，同一家族颜色一致。
                         aes(x = Gene_ID, y = 1, fill = Family)) + 
  geom_tile() +                                                               #数据点绘矩形色块                                                       
  scale_fill_discrete(type=c("#B3D3AA", "#4B6B5C", "#E88A71",                 #scale_fill_discrete自定义颜色
                              "#DEAB76", "#CD574D", "#85C1BF", "#BF38AE", 
                              "#176D84", "#7D83B7", "#4040C6", "#994B41")) + 
                              ggfun::theme_nothing()                          #隐藏所有的坐标轴、网格线和背景，只留下纯粹的瓷砖色块

##All-in-one integration to reveal transcription factor perturbation and subsequent biological effects

#Create a composite plot that integrates all the pieces of information

insert_top(tf_family_plot, perturbed_TF_plot, height = 5) |>                   #insert_top/insert_bottom 多图拼接函数
  insert_bottom(TF_GO_plot, height = 50)



###single-cell transcriptomic cell type annotation

#R packages
library(Seurat)                  #单细胞 RNA 测序数据分析：数据质控、标准化、降维、聚类到差异基因分析
library(CelliD)                  #用于无偏的细胞身份识别和跨数据集匹配
library(clusterProfiler)         #基因功能注释和富集分析的核心工具，支持 GO、KEGG、DO 等多种数据库，能对基因列表进行通路富集分析、基因集富集分析（GSEA），并提供火山图、气泡图等可视化结果
library(ggplot2)                 #核心绘图包
library(ggrepel)                 #图表中文本标签重叠
library(ggsc)
remove.packages("Seurat")
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))


#Import example PBMC data
pbmc_counts <- Read10X(data.dir = "hg19")
pbmc <- CreateSeuratObject(counts = pbmc_counts, project = "pbmc3k", min.cells = 3)  #Genes expressed in less than three cells are removed. #对象、项目、最小细胞数量

##Data preprocessing workflow

#Quality control by removing unwanted cells
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")                 #calculates the percentage of mitochondrial genes PercentageFeatureSet(pbmc, pattern = "^MT-")：计算每个细胞中线粒体基因的表达比例，^MT- 是匹配线粒体基因的正则表达式。
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 )  #subset()：根据设定的阈值过滤低质量细胞 &=or

#Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")                  #对数转化归一化
pbmc <- ScaleData(pbmc)                                                             #归一结果缩放处理

##Dimensionality reduction
#Identify highly variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)      #variance stabilizing transformation (vst) method & identifies the top 2,000 highly variable features

#Perform dimensional reduction by typing
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))                    #PCA可视化
pbmc <- RunUMAP(pbmc, dims = 1:10)                                                  #前十PCA主成分进行UMAP降维

##Cluster cells and identify markers of cell clusters

#Construct a K-nearest neighbor graph
pbmc <- FindNeighbors(pbmc, dims = 1:10)                                            #k-近邻法

#Cluster cells into different groups
pbmc <- FindClusters(pbmc, resolution = 0.5)                                        #成簇基因集

#Find cluster biomarkers
pbmc <- RunMCA(pbmc)                                                           #runs multiple correspondence analysis (MCA) to represent both cells and genes in the same space.
cluster_markers <- GetGroupGeneSet(pbmc, n.features = 20)                  #extracts the top 20 genes that are closest to each of the cell clusters

##Cell type annotation

#Import cell marker gene set
cell_marker_db <- read.gmt("c8.all.v2023.1.Hs.symbols.gmt")

#Perform cell type enrichment analysis by typing
cell_type_enrich_result <- compareCluster(cluster_markers, fun = "enricher", TERM2GENE = cell_marker_db )

#Predict cell type
predict_cell_type <- function(enrich_result) {
  enrich_result <- as.data.frame(enrich_result) 
  result <- split( 
    enrich_result, enrich_result$Cluster 
    ) |> vapply(function(x) { 
      x$ID[which.min(x$p.adjust)]
      }, FUN.VALUE = character(1)) 
  cell_type <- gsub("_", " ", result) |> 
    yulab.utils::str_wrap(18) 
  names(cell_type) <- names(result) 
  return(cell_type) 
  }  
cell_type_predict <- predict_cell_type(cell_type_enrich_result)

#Assign cell type identity to clusters
pbmc <- RenameIdents(pbmc, cell_type_predict)

#Visualize predicted cell type identity
cols <- c("#B3D3AA","#E88A71", "#DEAB76", "#CD574D",
           "#BF38AE", "#176D84", "#7D83B7", "#4040C6", "#994B41") 
           sc_dim(pbmc) + sc_dim_geom_label(geom = ggrepel::geom_text_repel, 
                                            color = "black", bg.color = "white") + 
             scale_color_discrete(type=cols) + 
             theme(legend.position = "none")

##总结
#基础R包参数学习：read.delim() 读取.txt 表格文件 、read.csv() 读取.csv文件；%in% 左侧是否存在于右侧；可视化点图dotplot;aes定义点的形状I使用（22）预设值：填充正方形...
#运行结果与文献结果存在偏差:可能原因：版本更迭、数据变更、数据库变化等，需注重Robust and Reproducible
           #代码精简准确，可重复性高。




























