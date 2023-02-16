# PORCUPINE
**P**rincipal Components Analysis to **O**btain **R**egulatory **C**ontributions **U**sing **P**athway-based Interpretation of **N**etwork **E**stimates is an R package to identify biological pathway which drive inter-tumour heterogeneity in a population of gene regulatory networks. It is a Principal Components Analysis (PCA)-based approach that can be used to determine whether a specific set of variables—for example a set of genes in a specific pathway—have coordinated variability in their regulation.

## Method
PORCUPINE uses as input individual patient networks, for example networks modeled using PANDA and LIONESS, as well as a .gmt file that includes biological pathways and the genes belonging to them. For each pathway, it extracts all edges connected to the genes belonging to that pathway and scales each edge across individuals. It then performs a PCA analysis on these edge weights, as well as on a null background that is based on random pathways. To identify significant pathways, PORCUPINE applies a one-tailed t-test and calculates the effect size (ES). 
Here we provide example how we analysed heterogeneity among single gene regulatory sample networks in Leiomyosarcomas (LMS). Networks were obtained with PANDA and LIONESS algorithms. Our example dataset contains data for 20 LMS patient specific gene regulatory networks. 

## Usage
```{r}
library(remotes)
install_github('kuijjerlab/PORCUPINE')
library("PORCUPINE")
```
For this analysis, we need to load the following packages:
```{r}
require("data.table")
require("fgsea")
require("dplyr")
require("plyr")
require("purrr")
require("stats")
require("parallel")
require("lsr")
```
First, we load the network data. Network file is quite big and it takes time to read it in. Our example contains 20 networks, represented by 11,151,077 edges.

```{r}
net <- fread(file.path(data_dir, "20LMS_nets.txt"))
net[1:5, 1:5]
#    0E244FE2-7C17-4642-A51F-2CCA796D9C70 75435ED8-93E8-45FB-8480-98D8EB2EF8CB
# 1:                                 0.76                                 0.10
# 2:                                 0.94                                 1.43
# 3:                                 1.09                                 2.78
# 4:                                 1.13                                 2.60
# 5:                                -0.71                                -1.42
#    B6D11678-15A9-4F43-A0A2-225067DCAF1C B7F5A41E-9559-4329-81F5-1B88A74730B7
# 1:                                -1.27                                 0.01
# 2:                                 0.30                                 0.91
# 3:                                 1.01                                 2.13
# 4:                                 1.66                                 1.71
# 5:                                 0.02                                 0.27
#    04823F53-A12D-4852-8F34-77B9DCBB7DF0
# 1:                                -7.18
# 2:                                -5.69
# 3:                                -6.09
# 4:                                -6.04
# 5:                                 5.31

dim(net)
# [1] 11151077       20

net <- t(net)

```
Next, we need to load in edges information, corresponding to our network data. Edges table should contain columns "tar" and "reg". In our case, edge table contains 623 Tfs and 17,899 target genes.

```{r}
edges <- fread(file.path(data_dir, "edges.txt"))
head(edges)
#      reg  tar prior
# 1:  AIRE A1BG     0
# 2:  ALX1 A1BG     0
# 3:  ALX3 A1BG     0
# 4:  ALX4 A1BG     0
# 5:    AR A1BG     0
# 6: ARID2 A1BG     0

length(unique(edges$reg))
# [1] 623
length(unique(edges$tar))
# [1] 17899
```
Then, we need to load in pathway file (.gmt file). Gmt files can be downloaded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
```{r}
pathways <- load_gmt(file.path(data_dir, "c2.cp.reactome.v7.1.symbols.gmt"))
length(pathways)
# [1] 1532
```
We need to filter pathways in order to include only pathways with genes that are present in our network file. 
```{r}
pathways <- filter_pathways(pathways, edges)
length(pathways)
```
Then, we filter pathways based on their size, and all pathways with less than 10 and more than 150 genes are filtered out. 
```{r}
pathways_to_use <- filter_pathways_size(pathways, minSize = 5, maxSize = 150 )
```
We select the top 10 pathways for the analysis (just to reduce computational time).
```{r}
pathways_to_use <- pathways_to_use[1:10]
```
Next, we perform PCA analysis for the 10 selected pathawys and extact the information on the variance explained by the first principal component in each pathway.
```{r}
pca_res_pathways <- pca_pathway(pathways_to_use, net, edges, ncores = 10, scale_data = TRUE)
head(pca_res_pathways)

#     pathway
# 1  REACTOME_GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS
# 2  REACTOME_PYRIMIDINE_CATABOLISM
# 3  REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS
# 4  REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE
# 5  REACTOME_APOPTOTIC_CLEAVAGE_OF_CELLULAR_PROTEINS
# 6  REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX
      pc1 n_edges pathway_size
# 1 22.911    8099           15
# 2 40.107    7476           12
# 3 31.126   12460           20
# 4 29.538   32396           55
# 5 24.147   21182           38
# 6 23.830   20559           33
```
Then we perform a PCA analysis based on random gene sets. In this case we create 500 random gene sets for each pathway and run PCA for each gene set.

```{r}
pca_res_random <- pca_random(net, edges, pca_res_pathways, pathways, n_perm = 500, ncores = 10, scale_data = TRUE)
head(pca_res_random)

#   pathway    pc1 n_edges pathway_size
# 1       1 40.891    7476           15
# 2       2 39.796    7476           15
# 3       3 34.702    8099           15
# 4       4 42.187    7476           15
# 5       5 20.051    8099           15
# 6       6 37.628    8722           15

```
Then to identify significant pathways we run PORCUPINE, which compares the observed PCA score for a pathway to a set of PCA scores of random gene sets of the same size as pathway. Calculates p-value and effect size.
```{r}

res_porcupine <- porcupine(pca_res_pathways, pca_res_random)
res_porcupine$p.adjust <- p.adjust(res_porcupine$pval, method = "fdr")
res_porcupine
#   pathway
# 1  REACTOME_BILE_ACID_AND_BILE_SALT_METABOLISM
# 2  REACTOME_BASE_EXCISION_REPAIR
# 3  REACTOME_PROCESSING_OF_INTRONLESS_PRE_MRNAS
# 4  REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX
# 5  REACTOME_APOPTOTIC_CLEAVAGE_OF_CELLULAR_PROTEINS
# 6  REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE
# 7  REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS
# 8  REACTOME_PYRIMIDINE_CATABOLISM
# 9  REACTOME_GAP_JUNCTION_DEGRADATION
# 10 REACTOME_GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS

# pathway_size          pval        es      p.adjust
# 1            43  9.999989e-01 0.2143374  1.000000e+00
# 2            91  4.780368e-05 0.1759077  1.593456e-04
# 3            19  1.000000e+00 1.1084317  1.000000e+00
# 4            33  1.000000e+00 2.0303499  1.000000e+00
# 5            38  1.000000e+00 1.9390723  1.000000e+00
# 6            55  1.000000e+00 0.3111659  1.000000e+00
# 7            20  9.994867e-01 0.1476974  1.000000e+00
# 8            12  4.812550e-97 1.1803796  2.406275e-96
# 9            12 2.044926e-259 3.1134399 2.044926e-258
# 10           15  1.000000e+00 1.9375777  1.000000e+00

```
