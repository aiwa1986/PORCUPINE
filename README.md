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
require("ggplot2")
require("gridExtra")
require(cowplot)
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
pathways_to_use <- filter_pathways_size(pathways, minSize = 5, maxSize = 150)
```
We select the top 10 pathways for the analysis (just to reduce computational time).
```{r}
pathways_to_use <- pathways_to_use[1:10]
```
Next, we perform PCA analysis for the 10 selected pathawys and extact the information on the variance explained by the first principal component in each pathway.
```{r}
pca_res_pathways <- pca_pathway(pathways_to_use, net, edges, ncores = 10, scale_data = TRUE)

head(pca_res_pathways)

                                                                  pathway
# 1                                           REACTOME_2_LTR_CIRCLE_FORMATION
# 2  REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS
# 3                                              REACTOME_ABACAVIR_METABOLISM
# 4                                 REACTOME_ABACAVIR_TRANSMEMBRANE_TRANSPORT
# 5                                REACTOME_ABACAVIR_TRANSPORT_AND_METABOLISM
# 6                           REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT
# 7                                        REACTOME_ABC_TRANSPORTER_DISORDERS
# 8                            REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS
# 9    REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT
# 10                     REACTOME_ACETYLCHOLINE_BINDING_AND_DOWNSTREAM_EVENTS
      pc1 n_edges pathway_size
# 1  26.146    4361            7
# 2  34.897   16198           26
# 3  36.172    3115            5
# 4  38.643    3115            5
# 5  36.113    6230           10
# 6  29.476   60431           97
# 7  31.932   44233           71
# 8  32.061   11214           18
# 9  31.695   13083           21
# 10 26.531    7476           12 
```
Then we perform a PCA analysis based on random gene sets. In this case we create 500 random gene sets for each pathway and run PCA for each gene set.

```{r}
pca_res_random <- pca_random(net, edges, pca_res_pathways, pathways, n_perm = 500, ncores = 10, scale_data = TRUE)
head(pca_res_random)

# pathway    pc1 n_edges pathway_size
# 1       1 26.748    4361            7
# 2       2 37.978    4361            7
# 3       3 39.253    4361            7
# 4       4 41.951    4361            7
# 5       5 41.715    4361            7
# 6       6 36.611    4361            7


```
Then to identify significant pathways we run PORCUPINE, which compares the observed PCA score for a pathway to a set of PCA scores of random gene sets of the same size as pathway. Calculates p-value and effect size.
```{r}

res_porcupine <- porcupine(pca_res_pathways, pca_res_random)
res_porcupine$p.adjust <- p.adjust(res_porcupine$pval, method = "fdr")
res_porcupine

#                                                                     pathway
# 1                      REACTOME_ACETYLCHOLINE_BINDING_AND_DOWNSTREAM_EVENTS
# 2    REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT
# 3                            REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS
# 4                                        REACTOME_ABC_TRANSPORTER_DISORDERS
# 5                           REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT
# 6                                REACTOME_ABACAVIR_TRANSPORT_AND_METABOLISM
# 7                                              REACTOME_ABACAVIR_METABOLISM
# 8                                 REACTOME_ABACAVIR_TRANSMEMBRANE_TRANSPORT
# 9  REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS
# 10                                          REACTOME_2_LTR_CIRCLE_FORMATION
#    pathway_size         pval         es     p.adjust
# 1            12 1.000000e+00 1.13481023 1.000000e+00
# 2            21 5.567509e-06 0.19851693 1.113502e-05
# 3            18 2.062801e-02 0.09151077 3.438001e-02
# 4            71 4.837247e-35 0.59456188 2.418624e-34
# 5            97 1.000000e+00 0.34364969 1.000000e+00
# 6            10 1.659546e-24 0.47864279 5.531820e-24
# 7             5 9.806794e-01 0.09272831 1.000000e+00
# 8             5 5.365370e-07 0.22087831 1.341342e-06
# 9            26 4.330923e-63 0.86664440 4.330923e-62
# 10            7 1.000000e+00 1.36020060 1.000000e+00

```
Get individual scores on PC1 and PC2 for a set of selected pathways
```{r}
ind_res <- get_pathway_ind_scores(pathways_to_use, net, edges,  scale_data = TRUE)
head(ind_res)

# $REACTOME_2_LTR_CIRCLE_FORMATION
#                                           Dim.1      Dim.2
# 0E244FE2-7C17-4642-A51F-2CCA796D9C70  -0.7609225  -1.474329
# 75435ED8-93E8-45FB-8480-98D8EB2EF8CB -26.2117484   5.977143
# B6D11678-15A9-4F43-A0A2-225067DCAF1C  -8.3714736  -4.078419
# B7F5A41E-9559-4329-81F5-1B88A74730B7  -8.7981762 -11.451997
# 04823F53-A12D-4852-8F34-77B9DCBB7DF0 -11.3803383 -21.492905
# 49684C2B-D31C-4B45-A400-3497C3CCEC01  -4.4498841 -38.838733
# FFDD7A12-DDEF-4974-8D60-64B7EEAAC994  -8.5446809   4.646410
# 830DFA6F-A85A-4317-82B2-791FAB998A01 -24.5263467  54.711749
# 58578614-E4A3-4655-BBAB-F65851625E0A   7.7084589   3.788487
# 4139E0C9-F712-4A25-8B59-587533B93B3E 122.3065345  20.052205
# 0C375B2F-67BE-4708-BEB2-544DEC812DCA  -2.3274860  -1.823275
# AB6324A1-19AB-400C-8001-54765D190E27 -11.9367843  -5.456453

clin_feature <- data.table("subtype" = c(rep("LMS", 10), rep("STLMS", 10)),
                  "count" =  rnorm(n = 20, mean = 0, sd = 1))
```
Plot PCA results for each pathway, color individuals according to selected feature,
first according to subtype, then according to count
```{r}
titles <- gsub("REACTOME_", "", names(ind_res))
titles <- substr(titles, start = 1, stop = 35)
titles_list <- as.list(titles)

plot_list <- list()
for (i in 1:length(ind_res)) {
      data <- as.data.frame(ind_res[[i]])
      plot_list[[i]] <- 
      plot_clinical_association(data, clin_feature$subtype)
} 
plot_list <- lapply(seq_along(plot_list), function(i) {
  ggdraw(plot_list[[i]]) +
    draw_label(titles_list[[i]], fontface = "bold", x = 0.5, y = 0.95)
})

grid.arrange(grobs = plot_list, ncol = 2)
<embed src="/images/subtype.pdf" type="application/pdf" width="50%" height="300px" />
```
Similar plot but for clinical feature "count"

```{r}
plot_list <- list()
for (i in 1:length(ind_res)) {
      data <- as.data.frame(ind_res[[i]])
      plot_list[[i]] <- 
      plot_clinical_association(data, clin_feature$count)
} 
plot_list <- lapply(seq_along(plot_list), function(i) {
  ggdraw(plot_list[[i]]) +
    draw_label(titles_list[[i]], fontface = "bold", x = 0.5, y = 0.95)
})
grid.arrange(grobs = plot_list, ncol = 2)
<embed src="/images/count.pdf" type="application/pdf" width="50%" height="300px" />
```
Perform Kruskal-Wallis or Pearson correlation tests (categorical and numerical features, respectively) between individual scores on PC1 (or PC2) and clinical feature.
Output is p-value for each pathway

For PC1: clinical_correlation(data, 1, clinical_feature)
For PC2: clinical_correlation(data, 2, clinical_feature)


```{r}
pval_dat_all <- NULL
for (i in 1:length(ind_res)) {
      data <- as.data.frame(ind_res[[i]])
      pval <- clinical_correlation(data, 1, clin_feature$count)
      pval_dat <- data.table("pathway" = names(ind_res[i]), "pval" = pval)
      pval_dat_all <- rbind(pval_dat, pval_dat_all)
      pval_dat_all
      } 
pval_dat_all
 #                                                                   pathway
 # 1:                     REACTOME_ACETYLCHOLINE_BINDING_AND_DOWNSTREAM_EVENTS
 # 2:   REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT
 # 3:                           REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS
 # 4:                                       REACTOME_ABC_TRANSPORTER_DISORDERS
 # 5:                          REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT
 # 6:                               REACTOME_ABACAVIR_TRANSPORT_AND_METABOLISM
 # 7:                                REACTOME_ABACAVIR_TRANSMEMBRANE_TRANSPORT
 # 8:                                             REACTOME_ABACAVIR_METABOLISM
 # 9: REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS
# 10:                                          REACTOME_2_LTR_CIRCLE_FORMATION
         pval
 # 1: 0.3226068
 # 2: 0.2487968
 # 3: 0.1296399
 # 4: 0.2317462
 # 5: 0.2233637
 # 6: 0.1246274
 # 7: 0.1606479
 # 8: 0.1058429
 # 9: 0.2389477
# 10: 0.1094883

pval_dat_all <- NULL
for (i in 1:length(ind_res)) {
      data <- as.data.frame(ind_res[[i]])
      pval <- clinical_correlation(data, 2, clin_feature$count)
      pval_dat <- data.table("pathway" = names(ind_res[i]), "pval" = pval)
      pval_dat_all <- rbind(pval_dat, pval_dat_all)
      pval_dat_all
      } 
pval_dat_all
```