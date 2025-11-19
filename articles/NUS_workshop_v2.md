# Unlocking single cell spatial omics analyses with scdney

**Presenting authors**

Harry Robertson$^{1,2,3}$, Lijia Yu$^{1,2}$, Beilei Bian$^{1,2}$, Andrew
Zhang$^{1,4}$, Jean Yang$^{1,2,3}$.

**Contributors**

Yue Cao$^{1,2}$, Lijia Yu$^{1,2}$, Andy Tran$^{1,2}$, Daniel
Kim$^{1,3}$, Dario Strbenac$^{1,2}$, Nicholas Robertson$^{1}$, Helen
Fu$^{3,4}$, Jean Yang$^{1,2,4}$.

$^{1}$ Sydney Precision Data Science Centre, University of Sydney,
Australia  
$^{2}$ School of Mathematics and Statistics, University of Sydney,
Australia  
$^{3}$ Charles Perkins Centre, University of Sydney, Australia  
$^{4}$ School of Computer Science, University of Sydney, Australia

  
Contact: jean.yang@sydney.edu.au

## Overview

The emergence of high-resolution spatial omics technologies has
revolutionized our ability to map cellular ecosystems in situ. This
workshop explores approaches in multi-sample spatial data analysis. We
cover end-to-end workflow considerations—from experimental design and QC
to spatial feature interpretation—with case studies in disease
prediction.

- [Our dataset](#tabset-1-1)
- [Assumed knowledge](#tabset-1-2)
- [Learning objectives](#tabset-1-3)

&nbsp;

- A subset of the widely-known METABRIC breast cancer cohort has
  recently been analyzed using imaging mass cytometry: [Imaging Mass
  Cytometry and Multiplatform Genomics Define the Phenogenomic Landscape
  of Breast Cancer](https://www.nature.com/articles/s43018-020-0026-6).
  Patient clinical data was sourced from the Supplementary Table 5 of
  [Dynamics of Breast-cancer Relapse Reveal Late-recurring ER-positive
  Genomic Subgroups](https://www.nature.com/articles/s41586-019-1007-8).

- Experience with R.
- Familiarity with the [SingleCellExperiment
  class](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html).
- Basic knowledge in [single cell data
  analysis](https://bioconductor.org/books/release/OSCA/index.html). You
  can access our [previous
  workshops](https://github.com/SydneyBioX/scdney#scdney-workshops-series)
  for a quick tutorial.
- Ability to install all required R packages, please check `sessionInfo`
  at the end of this document to ensure you are using the correct
  versions.
- Familiarity with our previous workshop vignette on [Introduction to
  Single Cell RNA-seq
  Analysis](https://sydneybiox.github.io/BIS2019_SC/).

- Describe and visualise spatial omics datasets.
- Understand approaches for quality control of spatial proteomics data.
- Calculate spatial statistics at the cell-type level using
  `scFeatures`.
- Perform multi-view disease outcome prediction with the package
  `ClassifyR`.
- Develop an understanding of:
  - evaluation of classification and survival models .
  - evaluate cohort heterogeneity given a survival model.
- Explore various strategies for disease outcome prognosis using spatial
  omics data.

**Our question of interest**

We want to know if the risk of recurrence in the METABRIC breast cancer
cohort can be accurately estimated to inform how aggressively they need
to be treated.

## Part 0: Loading an example image in R

To begin, we will start with loading an example image in R. It is not
neccessary for you to be familiar with this part, as we will provide you
with pre-processed data for the rest of the workshop.

system.file(“extdata”, “MB0002_1_345_fullstack.tiff”, package =
“scdneySpatial_ABACBS2025”)

- [DNA staining image](#tabset-2-1)
- [Mask](#tabset-2-2)
- [Overlay](#tabset-2-3)

&nbsp;

- Code
  ``` r
  suppressPackageStartupMessages({
    library(cytomapper)
    library(EBImage)
    library(S4Vectors)
  })


  img_path <- system.file("extdata/raw_image", "MB0002_1_345_fullstack.tiff",
                          package = "scdneySpatialABACBS2025")

  mask_path <- system.file("extdata/raw_image", "MB0002_1_345_cellmask.tiff",
                          package = "scdneySpatialABACBS2025")

  sample_id <- "MB0002_1_345"

  #Load image (multi-channel)
  images <- loadImages(img_path, single_channel = FALSE)
  names(images) <- sample_id
  #Load mask
  mask_img <- readImage(mask_path,as.is =T)
  masks <- CytoImageList(list(mask_img))
  names(masks) <- sample_id

  # Add metadata rows so img_id can be matched by COLUMN NAME
  mcols(images) <- DataFrame(sample_id = names(images))
  mcols(masks)  <- DataFrame(sample_id = names(masks))

  # Build channel names from CSVs and assign them
  metal_path <- system.file("extdata/raw_image", "channel_to_metal_order.csv",
                          package = "scdneySpatialABACBS2025")

  metal_order <- read.csv(metal_path,
                          header = FALSE, stringsAsFactors = FALSE)
  names(metal_order) <- "metal"
  metal_order$Channel <- seq_len(nrow(metal_order))

  marker_path <-  system.file("extdata/raw_image", "metal_to_marker.csv",
                          package = "scdneySpatialABACBS2025")
  metal_map <- read.csv(marker_path,
                        header = TRUE, stringsAsFactors = FALSE)

  merged <- merge(metal_order, metal_map, by = "metal", all.x = TRUE)
  merged <- merged[order(merged$Channel), ]
  labels <- ifelse(is.na(merged$marker), merged$metal, merged$marker)

  # Assign channel names to the multi-channel image
  channelNames(images) <- labels
  marker = plotPixels(
    image = images,
    colour_by = "Total HH3"
  )
  ```

  ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-2-1.png)

Code

``` r
# Plot the mask
plotCells(masks)
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-3-1.png)

Code

``` r
# Plot overlay 
combine = plotPixels(
  image    = images,
  mask     = masks,
  img_id   = "sample_id",      # column name in mcols(images)/mcols(masks)
  colour_by = "Total HH3"     
)
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-4-1.png)

## Part 1: Exploring the data

### 1.1: Initial data exploration

At the start of any analysis pipeline, it is often good to explore the
data to get a sense of the structure and its complexity. Let’s explore
the data to answer the questions below:

**Questions**

1.  How many features and observations are there in the data and what do
    the features represent?
2.  What covariates are in our data?
3.  Given our question of interest, what variable would be our outcome
    variable?

  

- [SCE object](#tabset-3-1)
- [Data overview](#tabset-3-2)
- [Covariates](#tabset-3-3)

&nbsp;

- First, we take a quick look at the structure of the SingleCell
  Experiment object.

  Code
  ``` r
  # Load data

  load(system.file("extdata", "breastCancer.RData",
                          package = "scdneySpatialABACBS2025"))
  data_sce = IMC

  # Structure of our data
  data_sce
  ```

      class: SingleCellExperiment
      dim: 38 76307
      metadata(0):
      assays(2): counts logcounts
      rownames(38): HH3_total CK19 ... H3K27me3 CK5
      rowData names(0):
      colnames(76307): MB-0002:345:93 MB-0002:345:107 ... MB-0663:394:577
        MB-0663:394:578
      colData names(17): file_id metabricId ... x_cord y_cord
      reducedDimNames(1): UMAP
      mainExpName: NULL
      altExpNames(0):

  Code
  ``` r
  # Metadata
  #DT::datatable(data.frame(colData(data_sce)))
  head(data.frame(colData(data_sce)))
  ```

                           file_id metabricId core_id ImageNumber ObjectNumber
      MB-0002:345:93  MB0002_1_345    MB-0002       1         345           93
      MB-0002:345:107 MB0002_1_345    MB-0002       1         345          107
      MB-0002:345:113 MB0002_1_345    MB-0002       1         345          113
      MB-0002:345:114 MB0002_1_345    MB-0002       1         345          114
      MB-0002:345:125 MB0002_1_345    MB-0002       1         345          125
      MB-0002:345:127 MB0002_1_345    MB-0002       1         345          127
                      Fibronectin Location_Center_X Location_Center_Y SOM_nodes
      MB-0002:345:93    0.4590055          179.7088          27.90110       130
      MB-0002:345:107   0.1176471          193.0588          29.29412       131
      MB-0002:345:113   0.2512903          162.1935          28.83871        83
      MB-0002:345:114   0.4788444          198.3111          29.15556       128
      MB-0002:345:125   6.1901112          262.8889          31.26667         7
      MB-0002:345:127   0.6897432          136.0946          32.98649       130
                      pg_cluster    description   region      ki67
      MB-0002:345:93          54       HR- CK7- region_1  low_ki67
      MB-0002:345:107         54       HR- CK7- region_4  low_ki67
      MB-0002:345:113         28    HRlow CKlow region_1 high_ki67
      MB-0002:345:114         54       HR- CK7- region_4 high_ki67
      MB-0002:345:125         11 Myofibroblasts region_2  low_ki67
      MB-0002:345:127         54       HR- CK7- region_1  low_ki67
                                    celltype  sample   x_cord   y_cord
      MB-0002:345:93        HR- CK7--region1 MB-0002 179.7088 27.90110
      MB-0002:345:107       HR- CK7--region4 MB-0002 193.0588 29.29412
      MB-0002:345:113    HRlow CKlow-region1 MB-0002 162.1935 28.83871
      MB-0002:345:114       HR- CK7--region4 MB-0002 198.3111 29.15556
      MB-0002:345:125 Myofibroblasts-region2 MB-0002 262.8889 31.26667
      MB-0002:345:127       HR- CK7--region1 MB-0002 136.0946 32.98649

Here, we take a quick look at what our rows and columns represent and
the dimensions of our data.

Code

``` r
# Glimpse the first few observations 
head(colnames(data_sce))
```

    [1] "MB-0002:345:93"  "MB-0002:345:107" "MB-0002:345:113" "MB-0002:345:114"
    [5] "MB-0002:345:125" "MB-0002:345:127"

Code

``` r
# Glimpse the first few features
head(rownames(data_sce))
```

    [1] "HH3_total" "CK19"      "CK8_18"    "Twist"     "CD68"      "CK14"     

Code

``` r
# How many features and observations are in our dataset?
dim(data_sce)
```

    [1]    38 76307

In addition to the proteomics data, it’s important to understand what
other covariates, such as clinical variables, are in our dataset. These
can help us in answering our question(s) of interest or formulate new
questions. For example, we can’t explore the association between smoking
status and breast cancer severity if the variable for smoking status
doesn’t exist in our data.

Code

``` r
# Explore covariates

# DT::datatable(clinical)
colnames(colData(data_sce))
```

     [1] "file_id"           "metabricId"        "core_id"
     [4] "ImageNumber"       "ObjectNumber"      "Fibronectin"
     [7] "Location_Center_X" "Location_Center_Y" "SOM_nodes"
    [10] "pg_cluster"        "description"       "region"
    [13] "ki67"              "celltype"          "sample"
    [16] "x_cord"            "y_cord"           

### 1.2: Is this a complex dataset?

Now that we have a basic idea of what our data looks like, we can look
at it in more detail. While initial data exploration reveals fundamental
patterns, deeper examination is very helpful. As it serves two critical
purposes: first, to detect anomalies or biases requiring remediation,
and second, to inform our choice of analytical methods tailored to the
biological questions at hand.

- [Missing values](#tabset-4-1)
- [Imbalance](#tabset-4-2)
- [Relationships](#tabset-4-3)
- [Visualise cells](#tabset-4-4)

&nbsp;

- Examining the proportion of missing values in a dataset is crucial for
  ensuring the accuracy and validity of data analysis. Missing values
  can significantly impact the reliability of statistical results,
  potentially leading to biased conclusions and reduced statistical
  power. Here, we plot the proportion of missing values for each of our
  clinical variables.

  Code
  ``` r
  # What is the proportion of missing values in our clinical data?
  gg_miss_var(clinical, show_pct = TRUE) + theme_bw()
  ```

  ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-7-1.png)

Looking at imbalance in data is crucial because it can lead to biased
models and inaccurate predictions, especially in classification tasks.
Imbalance occurs when one class is significantly underrepresented
compared to others. This can cause models to be overly influenced by the
majority class, leading to poor performance on the minority class. Here,
we tabularise some of the relevant clinical variables and plot the
distribution of the `Time to Recurrence-Free Survival` variable.

Code

``` r
# Estrogen receptor status
ggplot(clinical, aes(x = ER.Status)) + 
  geom_bar(fill = "#0099B4", color = "white", alpha = 0.8) +
  labs(y = "Count\n", 
       x = "\nER status") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-8-1.png)

Code

``` r
# Histological type
ggplot(clinical, aes(x = Histological.Type)) + 
  geom_bar(fill = "#fc6203", color = "white", alpha = 0.8) +
  labs(y = "Count\n", 
       x = "\nHistological Type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-8-2.png)

Code

``` r
# eventRFS: "Event in Recurrence-Free Survival."It indicates whether the event has occurred.#
ggplot(clinical, aes(x = factor(eventRFS))) + 
  geom_bar(fill = "#40dbb2", color = "white", alpha = 0.8) +
  labs(y = "Count\n", 
       x = "\nEvent in Recurrence-Free Survival") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-8-3.png)

Code

``` r
# timeRFS: "Time to Recurrence-Free Survival." It is the time period until recurrence occurs. 
ggplot(clinical, aes(x = timeRFS)) + 
  geom_histogram(fill = "#de9921", color = "white", alpha = 0.8, bins=20) +
  labs(y = "Frequency\n", 
       x = "\ntimeRFS") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-8-4.png)

Here, we explore the distribution of the outcomes and variables in the
meta-data. We use cross-tabulation to examine the following variables:
Surgery vs death, ER status, and Grade.

Code

``` r
# Stage and death
table(clinical$Breast.Surgery, clinical$Death, useNA = "ifany") 
```

                         0  1
      BREAST CONSERVING 38 14
      MASTECTOMY        16  6
      <NA>               2  1

Code

``` r
# ER status and grade
table(clinical$ER.Status, clinical$Grade)
```

           1  2  3
      neg  0  0 11
      pos 13 30 17

Code

``` r
# "Number of individuals based on Grade
table(clinical$Grade, clinical$Death)
```

         0  1
      1 12  2
      2 25  5
      3 17 11

To assess potential batch effects and sample-specific clustering
patterns, we visualize the cells in our data colored by the sample of
origin. This qualitative inspection provides an intuitive first
assessment of any batch effects in our data integration. If there is
batch effects then we’d need to apply batch correction to resolve this
issue.

Code

``` r
#data_sce <- runUMAP(data_sce, scale=TRUE)

# With the UMAP function we can highlight by meta data of interest
# Here we highlight the UMAP by sample ID
plotUMAP(data_sce, colour_by = "metabricId") + theme(legend.position = "none") +
  coord_equal()
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-10-1.png)

*Discussion*

- Are there any anomalies or biases that might need to be corrected to
  ensure our analysis is robust?
- Given our question of interest and the characteristics of our data,
  are there any particular analytic techniques that would be
  appropriate? Are there any that would not be appropriate?

  

## Part 2: Quality control

Here, we explore key approaches for evaluating the quality of IMC data.
A robust assessment should consider multiple factors, including the
density of marker expression, marker correlations, and their
co-expression. Something we might want to think about is what
characeristics might indicate whether a sample(s) is low/high quality?

Code

``` r
cell_counts <- IMC@colData |>  as.data.frame() |>
  dplyr::count(metabricId, name = "cell_count") |>  # Count rows per metabricId
  dplyr::arrange(desc(cell_count)) 

reducedDim(IMC, "spatialCoords") <- IMC@colData[,c("Location_Center_X","Location_Center_Y")]


cell_type_mapping <- c(
  "B cells" = "B cell",
  "T cells" = "T cell",
  "Macrophages Vim+ CD45low" = "Macrophage",
  "Macrophages Vim+ Slug-" = "Macrophage",
  "Macrophages Vim+ Slug+" = "Macrophage",

  "Endothelial" = "Endothelial",

  "Fibroblasts" = "Fibroblast",
  "Fibroblasts CD68+" = "Fibroblast",
  "Myofibroblasts" = "Fibroblast",
  "Vascular SMA+" = "Fibroblast",

  "Myoepithelial" = "Myoepithelial",

  "HR+ CK7-" = "Tumor HR+",
  "HR+ CK7- Ki67+" = "Tumor HR+",
  "HR+ CK7- Slug+" = "Tumor HR+",

  "HR- CK7+" = "Tumor HR-",
  "HR- CK7-" = "Tumor HR-",
  "HR- Ki67+" = "Tumor HR-",
  "HR- CKlow CK5+" = "Tumor HR-",

  "HRlow CKlow" = "Tumor HR-low/mixed",
  "HER2+" = "Tumor HER2+",
  "Basal CKlow" = "Tumor Basal-like",
  "Hypoxia" = "Tumor Hypoxic"
)

# Convert IMC description to higher-level categories
IMC$high_level_category <- recode(IMC$description, !!!cell_type_mapping)


imc.spe=IMC[,IMC$metabricId%in%cell_counts[1:20,]$metabricId]

# Convert IMC description to higher-level categories
imc.spe$high_level_category <- recode(imc.spe$description, !!!cell_type_mapping)


plots <- plot_marker_densities(imc.spe, sample_id_col = "metabricId",
                                markers_to_plot = c("CD20","CD3","CD68","vWF_CD31","Vimentin","HER2"),
                               assay_name = "logcounts")
```

### 2.1: QC for the panel / markers

Here, we want to assess the quality of the markers: it is desirable to
see at least two peaks, which would indicate the existence of cell type
specific markers in our data.

**Questions**

1.  Can we see any cell type specific markers, which are they?
2.  Are there any non-cell type specific markers?
3.  What can you notice about the distribution of markers across images?

We assess whether the marker intensities require transformation or
normalisation. This step is important for two main reasons:

- **Skewed distributions:** Marker intensities are often right-skewed,
  which can distort downstream analyses such as clustering or
  dimensionality reduction.

- **Inconsistent scales across images:** Due to technical variation, the
  same marker may show very different intensity ranges across images.
  This can shift what’s considered “positive” or “negative” expression,
  making it difficult to label cells accurately.

By applying transformation and normalisation, we aim to stabilise
variance and bring the data onto a more comparable scale across images.

  

- [CD20](#tabset-5-1)
- [CD3](#tabset-5-2)
- [CD68](#tabset-5-3)
- [vWF_CD31](#tabset-5-4)
- [Vimentin](#tabset-5-5)
- [HER2](#tabset-5-6)

&nbsp;

- ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-12-1.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-12-2.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-12-3.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-12-4.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-12-5.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-12-6.png)

#### Normalisation of marker expression

> **Tip**
>
> **What we’re looking for**
>
> 1.  Do the CD31+ and CD31- peaks clearly separate out in the density
>     plot? To ensure that downstream analysis goes smoothly, we want
>     our cell type specific markers to show 2 distinct peaks
>     representing our CD31+ and CD31- cells. If we see 3 or more peaks
>     where we don’t expect, this might be an indicator that further
>     normalisation is required.
> 2.  Are our CD31+ and CD31- peaks consistent across our images? We
>     want to make sure that our density plots for CD3 are largely the
>     same across images so that a CD3+ cell in one image is equivalent
>     to a CD3+ cell in another image.

Code

``` r
# leave out the nuclei markers from our normalisation process
useMarkers <- rownames(data_sce)[!rownames(data_sce) %in% c("DNA1", "DNA2", "HH3", "HH3_total", "HH3_ph")]

# transform and normalise the marker expression of each cell type
data_sce <- normalizeCells(data_sce,
                        markers = useMarkers,
                        transformation = NULL,
                        method = c("trim99", "minMax", "PC1"),
                        assayIn = "counts",
                        imageID = "metabricId")

selected_patients = which(colData(data_sce)$metabricId %in% cell_counts[1:10,]$metabricId)

norm = plot_marker_densities(data_sce[,selected_patients], sample_id_col = "metabricId", markers_to_plot = "vWF_CD31", assay_name = "norm")$vWF_CD31 + theme(legend.position = "none") + ggtitle("Normalized Counts - CD31")

log = plot_marker_densities(data_sce[,selected_patients], sample_id_col = "metabricId", markers_to_plot = c("vWF_CD31"), assay_name = "logcounts")$vWF_CD31 + theme(legend.position = "none") + ggtitle("Log Counts - CD31")

ggarrange(log, norm, ncol=2)
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-13-1.png)

In the plot above, the normalised data appeas more bimodal. We can
observe one clear CD31- peak at around 0.50, and a CD31+ peak at 1.00.
Image-level batch effects also appear to have been mitigated, since most
peaks occur at around the same CD31 intensity.

*Discussion*

- Do we have a sufficient amount of cell type specific markers?
- Is the panel of genes sufficient for our study?

  

### 2.2: QC of individual samples

Here we plot the correlation of different markers for a subset of
samples. This helps us determine whether a sample has low/high
contamination of markers.

**Questions**

1.  Which pairs of markers do we expect to be correlated/not correlated.
2.  Is there any evidence of marker contamination for a given sample?

  

Code

``` r
coexp_df <- readRDS(system.file("extdata", "coexp_df.rds",
                        package = "scdneySpatialABACBS2025"))


marker_list <- c("CD20", "CD3", "CD68", "vWF_CD31", "Vimentin", "HER2")
 
# Proportion heatmaps
heatmaps_prop <- plot_pairwise_heatmaps_per_sample(coexp_df, marker_list, stat = "proportion")
```

- [MB-0002](#tabset-6-1)
- [MB-0064](#tabset-6-2)
- [MB-0128](#tabset-6-3)
- [MB-0130](#tabset-6-4)
- [MB-0132](#tabset-6-5)
- [MB-0136](#tabset-6-6)
- [MB-0138](#tabset-6-7)
- [MB-0142](#tabset-6-8)
- [MB-0145](#tabset-6-9)
- [MB-0154](#tabset-6-10)
- [MB-0168](#tabset-6-11)
- [MB-0175](#tabset-6-12)
- [MB-0180](#tabset-6-13)
- [MB-0192](#tabset-6-14)
- [MB-0201](#tabset-6-15)
- [MB-0208](#tabset-6-16)
- [MB-0225](#tabset-6-17)
- [MB-0227](#tabset-6-18)
- [MB-0232](#tabset-6-19)
- [MB-0240](#tabset-6-20)
- [MB-0245](#tabset-6-21)
- [MB-0246](#tabset-6-22)
- [MB-0249](#tabset-6-23)
- [MB-0252](#tabset-6-24)
- [MB-0255](#tabset-6-25)
- [MB-0256](#tabset-6-26)
- [MB-0258](#tabset-6-27)
- [MB-0260](#tabset-6-28)
- [MB-0263](#tabset-6-29)
- [MB-0275](#tabset-6-30)
- [MB-0308](#tabset-6-31)
- [MB-0318](#tabset-6-32)
- [MB-0320](#tabset-6-33)
- [MB-0347](#tabset-6-34)
- [MB-0351](#tabset-6-35)
- [MB-0367](#tabset-6-36)
- [MB-0377](#tabset-6-37)
- [MB-0383](#tabset-6-38)
- [MB-0394](#tabset-6-39)
- [MB-0397](#tabset-6-40)
- [MB-0400](#tabset-6-41)
- [MB-0401](#tabset-6-42)
- [MB-0405](#tabset-6-43)
- [MB-0410](#tabset-6-44)
- [MB-0420](#tabset-6-45)
- [MB-0439](#tabset-6-46)
- [MB-0451](#tabset-6-47)
- [MB-0454](#tabset-6-48)
- [MB-0455](#tabset-6-49)
- [MB-0461](#tabset-6-50)
- [MB-0463](#tabset-6-51)
- [MB-0467](#tabset-6-52)
- [MB-0470](#tabset-6-53)
- [MB-0475](#tabset-6-54)
- [MB-0480](#tabset-6-55)
- [MB-0481](#tabset-6-56)
- [MB-0487](#tabset-6-57)
- [MB-0498](#tabset-6-58)
- [MB-0508](#tabset-6-59)
- [MB-0516](#tabset-6-60)
- [MB-0518](#tabset-6-61)
- [MB-0520](#tabset-6-62)
- [MB-0531](#tabset-6-63)
- [MB-0536](#tabset-6-64)
- [MB-0571](#tabset-6-65)
- [MB-0584](#tabset-6-66)
- [MB-0591](#tabset-6-67)
- [MB-0596](#tabset-6-68)
- [MB-0604](#tabset-6-69)
- [MB-0605](#tabset-6-70)
- [MB-0628](#tabset-6-71)
- [MB-0635](#tabset-6-72)
- [MB-0636](#tabset-6-73)
- [MB-0642](#tabset-6-74)
- [MB-0661](#tabset-6-75)
- [MB-0662](#tabset-6-76)
- [MB-0663](#tabset-6-77)

&nbsp;

- ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-1.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-2.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-3.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-4.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-5.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-6.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-7.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-8.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-9.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-10.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-11.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-12.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-13.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-14.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-15.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-16.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-17.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-18.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-19.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-20.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-21.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-22.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-23.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-24.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-25.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-26.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-27.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-28.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-29.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-30.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-31.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-32.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-33.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-34.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-35.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-36.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-37.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-38.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-39.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-40.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-41.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-42.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-43.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-44.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-45.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-46.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-47.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-48.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-49.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-50.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-51.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-52.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-53.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-54.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-55.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-56.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-57.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-58.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-59.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-60.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-61.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-62.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-63.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-64.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-65.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-66.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-67.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-68.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-69.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-70.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-71.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-72.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-73.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-74.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-75.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-76.png)

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-15-77.png)

### 2.3: QC of study as a whole

**Questions**

1.  Are there any samples you would like to remove from the data?
2.  Which samples would you like to keep?

  

Code

``` r
cell_counts <- IMC@colData |>  as.data.frame() |>
  dplyr::count(metabricId, name = "cell_count") |>  # Count rows per metabricId
  dplyr::arrange(desc(cell_count)) 


mean_coexp_per_sample <- coexp_df %>%
  dplyr::select(-cell_id, -cell_type) %>%  # remove columns not related to prob values
  group_by(metabricId) %>%
  summarise(mean_coexp_prob = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup()


merged_df <- left_join(mean_coexp_per_sample, cell_counts, by = "metabricId")

# Scatter plot
ggplot(merged_df, aes(x = cell_count, y = mean_coexp_prob)) +
  geom_point(size = 3, color = "steelblue") + 
  labs(
    x = "Cell Count per Sample",
    y = "Mean Co-expression",
    title = "Mean Co-expression vs. Cell Count"
  ) +
  theme_minimal(base_size = 14)
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-16-1.png)

### 2.4: Does my data agree with literature?

Another way to QC our data is to explore whether the DE genes or
associations we find matches current literature. For example, we have a
breast cancer cohort. Breast cancer is very well studied and so it
should be very easy to perform a DE analysis after pseudobulking our
data to see if it makes findings in previous bulk RNA-seq studies of
breast cancer cohorts. We can also do this at the single-cell level and
compare with previous single-cell proteomic/RNA-seq studies of breast
cancer cohorts. We won’t be running the code below. But here is an
example of how one might perform a basic DE analysis using package
[limma](https://academic.oup.com/nar/article/43/7/e47/2414268).

Code

``` r
# Note: The same approach can be used for pseudobulk or single-cell level data. 
# Just change `logcounts(data_sce)`

# We want low ki67 to be the reference level (baseline level) 
data_sce$ki67 <- factor(data_sce$ki67, levels=c("low_ki67", "high_ki67"))

# Here we specify the design matrix. We are specifying that Y are the expression values and X is ki67 status (but it could be any other variable - such as "good" or "bad" prognosis)
# Y~X: Expression as a function of ki67 status. 
design_matrix <- model.matrix(~data_sce$ki67)

# Fit model
fit <- lmFit(assay(data_sce, "norm"), design = design_matrix)

# Estimate variance and SE of coefficients
fit <- eBayes(fit)

tt <- topTable(fit, coef = 2, n = 5, adjust.method = "BH", sort.by = "p")
tt
```

                     logFC   AveExpr         t       P.Value     adj.P.Val
    HH3_total   2.99807594 5.3783014 126.61566  0.000000e+00  0.000000e+00
    HH3_ph      0.22158959 0.3549916  67.81304  0.000000e+00  0.000000e+00
    Ki67        0.07903174 0.1129605  65.26755  0.000000e+00  0.000000e+00
    E_cadherin  0.03507204 0.5389372  31.99330 4.085283e-223 3.881019e-222
    H3K27me3   -0.03479189 0.4609833 -30.33584 6.094169e-201 4.631568e-200
                       B
    HH3_total  7265.7780
    HH3_ph     2222.8176
    Ki67       2062.7347
    E_cadherin  498.5455
    H3K27me3    447.5419

### 2.5: Cell Type Classification

Here, we can also check whether the cell type annotations we have make
sense. One way to do this is to check whether the marker expression
levels agree with the cell type annotations. For example, we expect CD3
to be highly expressed in T cells, CD20 to be highly expressed in B
cells, CD68 to be highly expressed in Macrophages, vWF_CD31 to be highly
expressed in Endothelial cells, Vimentin to be highly expressed in
Fibroblasts, and HER2 to be highly expressed in Tumor HER2+ cells.

Code

``` r
# Calculate average expression per cell type
avg_expr <- aggregate(t(assay(data_sce, "norm")), 
                      by = list(CellType = IMC$high_level_category), 
                      FUN = mean) %>%
  tibble::column_to_rownames(var = "CellType")

avg_expr <- avg_expr[, c("CD20", "CD3", "CD68", "CD45", 
                         "vWF_CD31", "Vimentin", "HER2")]  # select only our mark

# Plot heatmap of marker by celltype
pheatmap(avg_expr, scale = "column", main = "Average Marker Expression per Cell Type\n", 
         cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(10000))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-18-1.png)

*Discussion*

- How might be check whether these results make sense? What potential
  issues might we find?

  

## Part 3: Exploring spatial data

This section assumes pre-existing cell type annotations. For unannotated
data, two primary annotation strategies are available:

- Unsupervised approach: cluster cells and identify marker genes for
  manual annotation.
- Supervised approach: transfer labels from reference datasets using
  supervised learning/classification tools.

For supervised annotation, we recommend
[scClassify](https://pmc.ncbi.nlm.nih.gov/articles/PMC7306901/), a
robust framework for cell-type classification. Below we visualise the
data with and without spatial information before exploring the data
further.

#### Cell types

Following initial assessment of potential batch effects through
sample-origin visualization, we now examine cell-type-specific
clustering patterns. Distinct, biologically meaningful clusters should
emerge for each annotated cell type. The presence of heterogeneous
clusters containing unrelated cell types may indicate incomplete or
inaccurate cell-type annotations.

Code

``` r
plotUMAP(data_sce, colour_by = "description")
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-19-1.png)

#### Spatial structure

The advantage with spatial omics is that we can examine the organisation
of the cell-types as it occurs on the tissue slide. Here, we visualise
one of the slides from a patient. We select a particular patient
“MB-0002” and visualise the tissue sample from this patient using
ggplot. Do we see any spatial patterning or does it look randomly
distributed?

Code

``` r
# obtaining the meta data for this patient 
one_sample <- data_sce[, data_sce$metabricId  == "MB-0002"]
one_sample  <- data.frame(colData(one_sample))
 
ggplot(one_sample, aes(x = Location_Center_X, y = Location_Center_Y, colour = description)) + 
  geom_point(alpha=0.7) + 
  theme(panel.background=element_blank(),
        axis.line=element_line(color="black")) +
  ggtitle("Original slide")
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-20-1.png)

**Questions**

1.  What kinds of spatial information might be of interest given the
    question we’d like to answer?
2.  How might we try to capture these spatial relationships?

  

- [Spatial regions](#tabset-7-1)
- [Spatial regions detection across multiple individuals
  (Banksy)](#tabset-7-2)
- [Across individuals (Optional)](#tabset-7-3)
- [Specific regions (Optional)](#tabset-7-4)

&nbsp;

- Background: One of the most common questions or analyses in spatial
  data is spatial domain detection. “Spatial domains” are regions within
  a tissue where cells share similar gene expression profiles and are
  physically clustered together. Here, we will use the terminology
  “spatial domain” and “regions” interchangeably. Most common analytical
  strategies involve spatial clustering, where different methods use
  different levels of information, such as gene expression data only,
  cell type information, and spatial coordinates.

  The strategy we demonstrate here uses the `lisaClust` function, which
  uses cell type information to cluster cells into five distinct spatial
  regions. As a case study, we will compare individuals with good or
  poor prognosis and examine them graphically to see if any regions
  appear to be different between good or poor prognosis. We define: -
  Good prognosis as individuals with \> 10 years recurrence-free
  survival and - Poor prognosis as individuals with \< 5 years
  recurrence-free survival.

  Below we visualise the spatial domain (regions) detection result based
  on one individual. Here we will use the terminology “spatial domain”
  and “regions” interchangeably.

  Depending on the number of regions, it may be more useful to visualise
  the spatial regions either collectively in a single graph or
  separately in multiple graphs. To visualise it in a single graph, the
  `hatchingPlot()` function is used to create hatching patterns for
  representating spatial regions and cell-types.

  Code
  ``` r
  ## Extract time to recurrence-free survival
  clinical$survivalgroup <- "neither"

  ## Define poor and good prognosis
  clinical$survivalgroup[which( clinical$timeRFS  < 5* 365) ] <- "poor" 
  clinical$survivalgroup[which( clinical$timeRFS  > 10* 365) ] <-  "good"

  ## To visualise it in a single graph
  hatchingPlot(
    data_sce,
    region = "region",
    imageID = "metabricId",
    cellType = "description",
    spatialCoords = c("Location_Center_X", "Location_Center_Y") ) 
  ```

  ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-21-1.png)

    

  We have written a small function `draw_region_clustering_result` to
  visualise the data separately in multiple graphs for the individual
  `MB-0002`.

  Code
  ``` r
  draw_region_clustering_result(data_sce , 
                                sample = "metabricId" , 
                                selected_sample =  "MB-0002" )
  ```

      [[1]]

  ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-22-1.png)

      [[2]]

  ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-22-2.png)

  *Discussion*

  - If we’d like to examine multiple individuals, is there a better
    visualise this information?

    

To perform spatial region detection across multiple individuals, we will
use the `Banksy` package. The `Banksy` method identifies spatial regions
by integrating both gene expression and spatial information. Below is
the code to perform spatial region detection using `Banksy` and
visualise the results.

Code

``` r
library(Banksy)
library(SpatialExperiment)
library(cowplot)

### The code below takes a while to run, so we have saved the output as an RDS file that can be loaded directly.

# spe_base <- SpatialExperiment(
#   assays = assays(data_sce),
#   rowData = rowData(data_sce),
#   colData = colData(data_sce),
#   spatialCoords = as.matrix(colData(data_sce)[, c("Location_Center_X", "Location_Center_Y")])
# )
# assay_to_use <- if ("counts" %in% assayNames(spe_base)) "counts" else assayNames(spe_base)[1]
# sample_ids <- unique(spe_base$metabricId)
# spe_list <- lapply(sample_ids, function(id) {
#   spe_subset <- spe_base[, spe_base$metabricId == id]
#   Banksy::computeBanksy(
#     spe_subset,
#     assay_name = assay_to_use,
#     k_geom = 18
#   )
# })
# spe_joint <- do.call(cbind, spe_list)
# spe_joint <- Banksy::runBanksyPCA(
#   spe_joint,
#   lambda = 0.8,
#   npcs = 30,
#   group = "metabricId"
# )
# spe_joint <- Banksy::clusterBanksy(
#   spe_joint,
#   lambda = 0.8,
#   npcs = 30,
#   algo = "kmeans",
#   kmeans.centers = 5
# )
# clust_col <- grep("^clust_", colnames(colData(spe_joint)), value = TRUE)[1]
# data_sce$banksy_region <- colData(spe_joint)[colnames(data_sce), clust_col]

banksy_ouput_spe_joint <- readRDS(system.file("extdata", 
                                              "banksy_ouput_spe_joint.rds",
                                              package = "scdneySpatialABACBS2025"))

colData(banksy_ouput_spe_joint)[,c("Location_Center_X", "Location_Center_Y")] <- NULL

hatchingPlot(
  banksy_ouput_spe_joint,
  region = "banksy_region",
  imageID = "metabricId",
  cellType = "description",
  spatialCoords = c("Location_Center_X", "Location_Center_Y")
)
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-23-1.png)

We can also visualise the spatial region detection results across
multiple individuals side by side for comparison. Below is an example of
visualising two individuals, `MB-0002` and `MB-0064`. Notice the
different cell type compositions in the different regions between the
two individuals.

Code

``` r
samples_to_plot <- c("MB-0002", "MB-0064")
plot_list <- lapply(samples_to_plot, function(id) {
  idx <- banksy_ouput_spe_joint$metabricId == id
  if (!any(idx)) return(NULL)
  hatchingPlot(
    banksy_ouput_spe_joint[, idx],
    region = "banksy_region",
    imageID = "metabricId",
    cellType = "description",
    spatialCoords = c("Location_Center_X", "Location_Center_Y")
  ) + ggtitle(id)
})
plot_list <- plot_list[!vapply(plot_list, is.null, logical(1))]
cowplot::plot_grid(plotlist = plot_list, ncol = 2)
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-24-1.png)

We can better interpret the region output by summarising the proportion
of each cell-type in a region across the individuals. We look at the
composition of cell-types in each region and compare between prognostic
outcomes.

Code

``` r
draw_dotplot(data_sce,  
             sample = "metabricId" , 
             celltype = "description" ,  
             group= "survivalgroup" , 
             group_of_interest =  c("poor" , "good"))
```

    `summarise()` has grouped output by 'Var1', 'Var2'. You can override using the
    `.groups` argument.

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-25-1.png)

The number of sub-cell types increase considerably when we want to add
spatial domain (region) information. To enhance clarity and facilitate
understanding, it may be helpful to choose a predetermined region. The
code generates a set of boxplots that enable the comparison of cell-type
proportions between individuals with good and poor prognosis in
`region_5`.

Code

``` r
draw_selected_region_boxplot(data_sce, 
                             sample = "metabricId" , 
                             celltype ="description" , 
                             group  = "survivalgroup", 
                             group_of_interest =  c("poor" , "good"), 
                             select_region = "region_5")
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-26-1.png)

## Part 4: Feature engineering with scFeatures

Here, we use scFeatures to generate molecular features for each
individual using the `features x cells` matrices. These features are
interpretable and can be used for downstream analyses. In general,
`scFeatures` generates features across six categories:

1.  Cell-type proportions.
2.  Cell-type specific gene expressions.
3.  Cell-type specific pathway expressions.
4.  Cell-type interaction scores.
5.  Aggregated gene expressions.
6.  Spatial metrics: Nearest neighbour’s correlation, L statistics, and
    Moran’s I.

The different types of features constructed enable a more comprehensive
multi-view understanding of each individual. By default, the function
will generate a total of 13 different types of features and are stored
as 13 `samples x features` matrices, one for each feature type.

In this section, we will examine spatial information from two
perspectives. Utilising spatial domain detection described in part 3, we
select a specific region of interest and create molecular
representations of that region for each individual (Section 4.1).
Second, we will utilise spatial statistics to capture spatial
relationships within the region of interest, such as Moran’s I (Section
4.2)

Code

``` r
region <- data_sce$region

# Define a series of sub-cell-types that is regional specific
data_sce$celltype <- paste0( data_sce$description , "-" , region)
```

### 4.1 How to create molecular representations of individuals for an ROI?

Here, we can consider regional information and the following code allows
us to create cell-type specific features for each region. We use the
function “paste0” to construct **region-specific sub cell-types** and
name it as `celltype` in the R object `data_sce`. For simplicity, in
this workshop, the variable `celltype` in the R object `data_sce` refers
to region-specific sub-cell-types.

There are a total of 13 different types of features (feature_types) that
you can choose to generate. The argument type refers the type of input
data we have. This is either `scrna` (single-cell RNA-sequencing data),
`spatial_p` (spatial proteomics data), or `spatial_t` (single-cell
spatial data). In this section, we will ignore spatial information and
generate non-spatial features, such as pseudobulking at the sample /cell
type levels, overall expression, cell type proportions etc…

- [Generate features](#tabset-8-1)
- [ROI](#tabset-8-2)
- [Association Study](#tabset-8-3)
- [scFeatures code (Optional)](#tabset-8-4)
- [scFeatures output (Optional)](#tabset-8-5)

&nbsp;

- Suppose that we are interested in determining the proportion of
  cell-types within each region for each individual. It is necessary to
  specify `type = spatial_p` to reflect that we have spatial proteomics
  data and `feature_types = proportion_raw` to indicate we intend to
  calculate cell-type proportions for each of the region-specific
  cell-types. Suppose we are only interested in the molecular
  representation of `HR- CK7+` within individuals within region 5.

  Code
  ``` r
  ## [A] The next few lines extract specific information from data_sce as input to scFeatures. 
  ## Extract the expression matrix from data_sce
  IMCmatrix <- assay(data_sce, "norm")

  ## Extract the sample information 
  ## append the condition to the individuals so we can easily retrieve the individuals condition 
  sample <- data_sce$metabricId 
  cond  <- clinical[match(sample, clinical$metabricId), ]$survivalgroup
  sample <- paste0(sample, "_cond_", cond ) 

  ## Extract the region-specific sub-cell-types
  celltype <- data_sce$celltype

  ## Extract the spatial coordinates
  spatialCoords <- list(colData(data_sce)$Location_Center_X, colData(data_sce)$Location_Center_Y)

  ### [B] Running scFeatures
  scfeatures_result <- scFeatures(IMCmatrix, 
                                  sample = sample, 
                                  celltype = celltype, 
                                  spatialCoords = spatialCoords,
                                  feature_types = "proportion_raw", 
                                  type = "spatial_p" )
  ```

**Question**

1.  Are there any regions that are associated with “good” and “poor”
    prognosis?
2.  Is this the right way to visualise the results?

  

Code

``` r
## [A] The next few lines extract specific information from data_sce as input to scFeatures. 
## Select the HR- CK7+-region sub-cell-type 

# There are different ways you can use `scFeatures` to generate molecular representations for individuals and it requires the following information for spatial data.
# 
# data,\
# sample,
# X coordinates,
# Y coordinates,
# feature_types, and
# type

index <- grep("HR- CK7+-region" , data_sce$celltype, fixed=TRUE)
selected_data <- IMCmatrix[, index]
selected_sample <- sample[index]
selected_celltype <- data_sce$celltype[ index] 
selected_spatialCoords <- list(colData(data_sce)$Location_Center_X[index], 
                               colData(data_sce)$Location_Center_Y[index])

### [B] Running scFeatures
scfeatures_result <- scFeatures( selected_data, 
                                 sample = selected_sample, 
                                 celltype = selected_celltype,
                                 spatialCoords = selected_spatialCoords,
                                 feature_types = "proportion_raw", type = "spatial_p" )

### [C] Visualize the regional composition makeup for each individual for HR- CK7+ and HR- CK7-
feature <- scfeatures_result$proportion_raw
feature <- feature[ grep("poor|good", rownames(feature)),  ]

plot_barplot(feature ) + ggtitle("Proportion raw feature") + labs(y = "proportion\n")
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-29-1.png)

Can we identify “differential expression” for a feature of interest? The
R object scfeatures_result contains a variety of features. A important
question focuses on the identification of features that reflect an
association with the prognostic outcome, specifically distinguishing
between good and poor outcomes. The code provided below demonstrates the
use of the limma() function to fit a linear model for the purpose of
analysing gene_mean_celltype as an illustration feature type. The
feature type known as gene_mean_celltype represents the mean protein
expression for each sub-cell-type specific to a spatial region. It is a
matrix consisting of 77 individuals and 4180 features. It is important
to acknowledge that within the context of our IMC data, the term “gene”
is used to refer to “protein”.

Code

``` r
scfeatures_result <- readRDS(system.file("extdata", 
                                              "scfeatures_result.rds",
                                              package = "scdneySpatialABACBS2025"))
# Extract cell-type specific gene expression for all regions. 
gene_mean_celltype <- scfeatures_result$gene_mean_celltype

# Extract HR+ CK7 cell-type specific gene expression for Region5
index <-  grep("HR+ CK7--region5", colnames(gene_mean_celltype) , fixed= T)
gene_mean_celltype <- gene_mean_celltype [, index] 

# transpose to ensure we have gene by sample matrix
gene_mean_celltype <- t(gene_mean_celltype)
      
# Extract the two conditions of interest - poor prognosis vs good prognosis
condition  <- unlist( lapply( strsplit( colnames(gene_mean_celltype) , "_cond_"), `[`, 2))
condition <- data.frame(sample = colnames(gene_mean_celltype), condition = condition )
select_index <- which( condition$condition %in% c("poor",  "good" ))
condition <- condition[ select_index, ]
gene_mean_celltype<- gene_mean_celltype [ ,  select_index]


# Calculate log fold change each protein using limma
design <- model.matrix(~condition, data = condition)
fit <- lmFit(gene_mean_celltype, design)
fit <- eBayes(fit)
tT <- topTable(fit, n = Inf)
tT$gene <- rownames(tT)
tT[1:6] <- signif(tT[1:6], digits=2)
DT::datatable(tT)
```

We visualise the comparison using a volcano plot and a dotplot for the
cell-type specific expression feature. This is a type of scatter-plot
that is used to quickly identify changes in large datasets and represent
the significance (y-axis) versus effect size or fold-change (x-axis).

Code

``` r
# order the proteins by log fold change 
tT <- tT[ order(tT$logFC, decreasing = T), ]
tT <- tT[1:20, ]
ggplot( tT , aes( y = reorder(gene, logFC) , x = logFC  ) )+
      geom_point(aes(colour=-log10(P.Value)), alpha=2/3, size=4) +
      scale_colour_gradient(low="blue",high="red") + theme_bw() + 
      xlab("logFC") + ylab("region specific cell type specfic features" ) 
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-31-1.png)

The code below enable us to generate all feature types for all
cell-types in a line. Due to limitations with today’s computational
capacity, **Please DO NOT run it in today’s workshop, it will crash your
system**.

Code

``` r
# here, we specify that this is a spatial proteomics data
# scFeatures support parallel computation to speed up the process 
# scfeatures_result <- scFeatures(IMCmatrix, 
#                                 type = "spatial_p",
#                                 sample = sample, 
#                                 celltype = celltype, 
#                                 spatialCoords = spatialCoords,
#                                 ncores = 32)
```

Assuming you have already generated a collection of molecular
representation for individuals, please load the prepared RDS file
`scfeatures_result.rds`. Again, you can remind yourself that all
generated feature types are stored in a matrix of `samples x features`.

Code

``` r
# Upload pre-generated RDS file


scfeatures_result <- readRDS(system.file("extdata", 
                                              "scfeatures_result.rds",
                                              package = "scdneySpatialABACBS2025"))

# What are the features and the dimensions of features matrices that we have generated?
lapply(scfeatures_result, dim)
```

    $proportion_raw
    [1]  77 110

    $proportion_logit
    [1]  77 110

    $proportion_ratio
    [1]   77 5995

    $gene_mean_celltype
    [1]   77 4180

    $gene_prop_celltype
    [1]   77 4180

    $gene_cor_celltype
    [1]    77 27958

    $gene_mean_bulk
    [1] 77 38

    $gene_prop_bulk
    [1] 77 38

    $gene_cor_bulk
    [1]  77 703

    $L_stats
    [1]   77 3842

    $celltype_interaction
    [1]   77 4393

    $morans_I
    [1] 77 38

    $nn_correlation
    [1] 77 38

### 4.2 How can we represent spatial information and relationships for a given individual?

We will now look at the spatial statistic output by `scFeatures` and
qualitively assess whether there is any association between between
these statistics and the “good” and “bad” prognosis groups.

**Questions**

1.  What kind of information would spatial statistics provide?
2.  Are there any differences in the distribution of spatial statistics
    between the “good” and “bad” prognosis groups?

  

- [Nearest neighbour correlation](#tabset-9-1)
- [L statistics](#tabset-9-2)
- [Moran’s I](#tabset-9-3)

&nbsp;

- This metric measures the correlation of proteins/genes between a cell
  and its nearest neighbour cell.

  Code
  ``` r
  feature <- scfeatures_result$nn_correlation
  feature <- feature[ grep("poor|good", rownames(feature )),  ]

  top_features <- feature %>%
    mutate(Group = ifelse(grepl("good", rownames(feature)), "Good Prognosis", "Poor Prognosis")) %>%
    pivot_longer(cols = -Group, names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      cohens_d = (mean(value[Group == "Good Prognosis"]) - mean(value[Group == "Poor Prognosis"])) /sd(value),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(cohens_d))) %>%
    slice_head(n = 5) %>%
    pull(variable)

  feature %>%
    mutate(Group = ifelse(grepl("good", rownames(feature)), "Good Prognosis", "Poor Prognosis")) %>%
    pivot_longer(cols = -Group, names_to = "variable", values_to = "value") %>%
    filter(variable %in% top_features) %>%  # Keep only top 5
    ggplot(aes(y = variable, x = value, fill = Group)) + 
    geom_boxplot() + theme_bw() + labs(x = "Nearest Neighbour Correlation", y = "Marker")
  ```

  ![](NUS_workshop_v2_files/figure-html/unnamed-chunk-34-1.png)

The L function is a spatial statistic used to assess the spatial
distribution of cell-types. It assesses the significance of cell-cell
interactions, by calculating the density of a cell-type with other
cell-types within a certain radius. High values indicate spatial
association (co-localisation), low values indicate spatial avoidance. To
demonstrate the L-function, we will plot a specific patient “MB-0128”
who has a high L value for B cells interacting with Fibroblasts and a
low L value for B cells interacting with HR- CK7- cells.

Code

``` r
tableau_palette <- scale_colour_tableau( palette = "Tableau 20")
color_codes <- tableau_palette$palette(10)

# Create a named color vector
cell_colors <- c(
  "B cells" = color_codes[1],
  "Fibroblasts" = color_codes[2],
  "HR- CK7-" = color_codes[9],
  "others" = color_codes[4]
)

# select one patient 
one_sample  <- data_sce[ , data_sce$metabricId == "MB-0128"  ]
one_sample <- data.frame( colData(one_sample) )

one_sample$celltype <- one_sample$description

# High L-function value plot.
one_sample <- data.frame(colData(data_sce[, data_sce$metabricId == "MB-0128"]))
one_sample$celltype <- one_sample$description
index <- one_sample$celltype %in% c("B cells", "Fibroblasts")
one_sample$celltype[!index] <- "others"
a <- ggplot(one_sample, aes(x = Location_Center_X, y = Location_Center_Y, colour = celltype)) + 
  geom_point() + theme(panel.background=element_blank(),
        axis.line=element_line(color="black")) +
  scale_colour_manual(values = cell_colors) +  # Use named vector
  ggtitle("Patient MB-0128 - high L value with \n B cells interacting Fibroblasts")

# Low L-function value plot.
one_sample$celltype <- one_sample$description
index <- one_sample$celltype %in% c("Fibroblasts", "HR- CK7-")
one_sample$celltype[!index] <- "others"
b <- ggplot(one_sample, aes(x = Location_Center_X, y = Location_Center_Y, colour = celltype)) + 
  geom_point() + theme(panel.background=element_blank(),
        axis.line=element_line(color="black")) +
  scale_colour_manual(values = cell_colors) +  # Use named vector
  ggtitle("Patient MB-0128 - low L value with \n B cells interacting HR_ CK7")

ggarrange(plotlist = list(a, b))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-35-1.png)

Moran’s I is a spatial autocorrelation statistic based on both location
and values. It quantifies whether similar values tend to occur near each
other or are dispersed.

Code

``` r
feature <- scfeatures_result$morans_I
feature <- feature[ grep("poor|good", rownames(feature )),  ]

top_features <- feature %>%
  mutate(Group = ifelse(grepl("good", rownames(feature)), "Good Prognosis", "Poor Prognosis")) %>%
  pivot_longer(cols = -Group, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    cohens_d = (mean(value[Group == "Good Prognosis"]) - mean(value[Group == "Poor Prognosis"])),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(cohens_d))) %>%
  slice_head(n = 5) %>%
  pull(variable)

feature %>%
  mutate(Group = ifelse(grepl("good", rownames(feature)), "Good Prognosis", "Poor Prognosis")) %>%
  pivot_longer(cols = -Group, names_to = "variable", values_to = "value") %>%
  filter(variable %in% c(top_features, "Ki67")) %>%  # Keep only top 5
  ggplot(aes(y = reorder(variable,value), x = value, fill = Group)) + 
  geom_boxplot() + theme_bw() + labs(x = "Moran's I", y = "Marker")
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-36-1.png)

## Part 5: Disease classification with ClassifyR \[Presentation\]

Recurrence risk estimation is a fundamental concern in medical research,
particularly in the context of patient survival analysis. In this
section, we will estimate recurrence risk using the molecular
representation of individuals generated from `scFeatures` to build a
survival model. We will use classifyR to build the survival model. The
patient outcome is time-to-event, so, by default, ClassifyR will use Cox
proportional hazards ranking to choose a set of features and also Cox
proportional hazards to predict risk scores. We will also demonstrate
other available models in `ClassifyR`.

**Questions**

1.  Are spatial features globally informative in predicting survival?
2.  If not, for which individuals is it important in predicting
    survival?

  

- [Building a survival model (Optional)](#tabset-10-1)
- [Model performance (Optional)](#tabset-10-2)
- [Does spatial information matter?](#tabset-10-3)

&nbsp;

- Recall in the previous section that we have stored the 13 matrices of
  different feature types in a list. Instead of manually retrieving each
  matrix from the list to build separate models, classifyR can directly
  take a list of matrices as an input and run a repeated
  cross-validation model on each matrix individually. Below, we run 5
  repeats of 5-fold cross-validation. A high score indicates prognosis
  of a worse outcome than a lower risk score. Although we have provided
  the code below, to save time, **just load the prepared RDS file
  `classifyr_result_IMC.rds`** and we will focus on the interpretation
  in this workshop.

      # We use the following variables:
      # timeRFS: "Time to Recurrence-Free Survival." It is the time period until recurrence occurs.
      # eventRFS: "Event in Recurrence-Free Survival."It indicates whether the event has occurred.
      # Breast.Tumour.Laterality: Laterality of tumors, eg, whether the tumor is located in left or right.
      # ER.Status: Whether the tumor is ER positive or ER negative.
      # Inferred.Menopausal.State: of the patient.
      # Grade: of the tumor.
      # Size: of the tumor.

      usefulFeatures <- c("Breast.Tumour.Laterality", "ER.Status", "Inferred.Menopausal.State", "Grade", "Size")
      nFeatures <- append(list(clinical = 1:3), lapply(scfeatures_result, function(metaFeature) 1:5))
      clinicalAndOmics <- append(list(clinical = clinical), scfeatures_result)

      ### generate classfyr result
      classifyr_result_IMC <- crossValidate(clinicalAndOmics, c("timeRFS", "eventRFS"),
                          extraParams = list(prepare = list(useFeatures = list(clinical = usefulFeatures))),
                          nFeatures = nFeatures, nFolds = 5, nRepeats = 5, nCores = 5)

  Code
  ``` r
  classifyr_result_IMC <- readRDS(system.file("extdata", 
                                                "classifyr_result_IMC.rds",
                                                package = "scdneySpatialABACBS2025"))
  ```

  Cox proportional hazards is a classical statistical method, as opposed
  to machine learning methods like Random survival forest. These machine
  learning methods can build remarkably complex relationships between
  features, however their running time can be much longer than Cox
  proportional hazards. We use feature selection to limit the number of
  features considered to at most 100 per metafeature and to save time,
  **you can just load the prepared RDS file.** We will compare the
  predictive performance between these methods.

      nFeatures <- append(list(clinical = 1:3), lapply(scfeatures_result[2:length(scfeatures_result)], function(metaFeature) min(100, ncol(metaFeature))))
      survForestCV <- crossValidate(clinicalAndOmics, outcome, nFeatures = nFeatures,
                      classifier = "randomForest",
                      nFolds = 5, nRepeats = 5, nCores = 5)

  Code
  ``` r
  survForestCV <- readRDS(system.file("extdata", 
                                                "survForestCV.rds",
                                                package = "scdneySpatialABACBS2025"))
  ```

To examine the distribution of prognostic performance, use
`performancePlot`. Currently, the only metric for time-to-event data is
C-index and that will automatically be used because the predictive model
type is tracked inside of the result objects.

Code

``` r
## Make axis label 45 degree to improve readiability 
tilt <- theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

## Putting two sets of cross-validation results together
multiresults <- append(classifyr_result_IMC, survForestCV)
ordering <- c("clinical", names(scfeatures_result))
performancePlot(multiresults,
                characteristicsList = list(x = "Assay Name", 
                                           row = "Classifier Name"),
                orderingList = list("Assay Name" = ordering)) + 
                tilt
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-39-1.png)

Note how the resultant plot is a `ggplot2` object and can be further
modified. The same code could be used for a categorical classifier
because the random forest implementation provided by the `ranger`
package has the same interface for both. We will examine feature
selection stability with `selectionPlot`.

Code

``` r
selectionPlot(multiresults,
                characteristicsList = list(x = "Assay Name", row = "Classifier Name"),
                orderingList = list("Assay Name" = ordering)) + tilt
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-40-1.png)

Code

``` r
distribution(classifyr_result_IMC[[1]], plot = FALSE)
```

         assay                   feature proportion
    1 clinical Inferred.Menopausal.State          1

Does each individual require a different collection of features? Using
`samplesMetricMap` compare the per-sample C-index for Cox models for all
kinds of metafeatures.

Code

``` r
library(grid)
samplesMetricMap(classifyr_result_IMC)
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-41-1.png)

    TableGrob (2 x 1) "arrange": 2 grobs
      z     cells    name                 grob
    1 1 (2-2,1-1) arrange       gtable[layout]
    2 2 (1-1,1-1) arrange text[GRID.text.8777]

## Appendix

### Explanation of spatial features

- L function:

The L function is a spatial statistic used to assess the spatial
distribution of cell-types. It assesses the significance of cell-cell
interactions, by calculating the density of a cell-type with other
cell-types within a certain radius. High values indicate spatial
association (co-localisation), low values indicate spatial avoidance.

Code

``` r
tableau_palette <- scale_colour_tableau( palette = "Tableau 20")
color_codes <- tableau_palette$palette(10)
# select one patient 
one_sample  <- data_sce[ , data_sce$metabricId == "MB-0128"  ]
one_sample <- data.frame( colData(one_sample) )

one_sample$celltype <- one_sample$description

# select certain cell types to examine the interaction 
index <-  one_sample$celltype  %in% c("B cells", "Fibroblasts")
one_sample$celltype[!index] <- "others"
a <-ggplot( one_sample, aes(x = Location_Center_X , y = Location_Center_Y, colour = celltype )) + geom_point()  + scale_colour_manual(values = color_codes)  + ggtitle( "Patient MB-0128 - high L value with \n B cells interacting Fibroblasts")
 

one_sample$celltype <- one_sample$description
index <-  one_sample$celltype  %in% c("Fibroblasts", "HR- CK7-")
one_sample$celltype[!index] <- "others"
b <- ggplot( one_sample, aes(x = Location_Center_X , y = Location_Center_Y, colour = celltype )) + geom_point()  + scale_colour_manual(values = color_codes)  + ggtitle( "Patient MB-0128 - low L value with \n B cells interacting HR_ CK7")
 
ggarrange(plotlist = list(a,b))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-42-1.png)

- Cell type interaction composition:

We calculate the nearest neighbours of each cell and then calculate the
pairs of cell-types based on the nearest neighbour. This allows us to
summarise it into a cell-type interaction composition.

Code

``` r
 one_sample  <- data_sce[ , data_sce$metabricId == "MB-0263"  ]
one_sample <- data.frame( colData(one_sample) )

one_sample$celltype <- one_sample$description
 
a <-ggplot( one_sample, aes(x = Location_Center_X , y = Location_Center_Y, colour = celltype )) + geom_point()  + scale_colour_manual(values = color_codes)  + ggtitle("Patient MB-0263")

feature  <- scfeatures_result$celltype_interaction
to_plot <- data.frame( t( feature[ "MB-0263_cond_poor" , ])  )
to_plot$feature <- rownames(to_plot) 
colnames(to_plot) <- c("value", "celltype interaction composition")
 
to_plot <- to_plot[ order(to_plot$value, decreasing = T), ]
to_plot <- to_plot[1:10, ]
to_plot$`celltype interaction composition` <- factor(to_plot$`celltype interaction composition`, levels = to_plot$`celltype interaction composition`)

b <- ggplot(to_plot, aes(x =  `celltype interaction composition`  ,  y = value, fill=`celltype interaction composition`)) + geom_bar(stat="identity" ) + ylab("Major cell-type interactions")  +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggarrange(plotlist = list(a,b))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-43-1.png)

- Moran’s I:

Moran’s I is a spatial autocorrelation statistic based on both location
and values. It quantifies whether similar values tend to occur near each
other or are dispersed.

Code

``` r
high  <- data_sce["Ki67", data_sce$metabricId == "MB-0132"  ]
high_meta <- data.frame( colData(high) ) 
high_meta$expression <- as.vector(logcounts( high)) 

low  <- data_sce["Ki67",  data_sce$metabricId == "MB-0249" ]
low_meta <- data.frame( colData(low) )
low_meta$expression <- as.vector(logcounts(low))


a <- ggplot(high_meta, aes(x = Location_Center_X , y = Location_Center_Y, colour =expression)) + geom_point(alpha=0.5) + scale_colour_viridis_c() + ggtitle("Patient MB-0132 - high Moran's I in Ki67")

b <- ggplot(low_meta, aes(x = Location_Center_X , y = Location_Center_Y, colour =expression)) + geom_point(alpha=0.5) + scale_colour_viridis_c() + ggtitle("Patient MB-0249 - low Moran's I in Ki67")

ggarrange(plotlist = list(a,b))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-44-1.png)

- Nearest Neighbor Correlation:

This metric measures the correlation of proteins/genes between a cell
and its nearest neighbour cell.

Code

``` r
plot_nncorrelation <- function(thissample , thisprotein){
   
       sample_name <- thissample
       thissample <- data_sce[, data_sce$metabricId ==     sample_name]
    
      
      exprsMat <- logcounts(thissample)
     
    # calculate NN correlation 
    cell_points_cts <- spatstat.geom::ppp(
            x = as.numeric(thissample$Location_Center_X ), y = as.numeric(thissample$Location_Center_Y),
            check = FALSE,
            xrange = c(
                min(as.numeric(thissample$Location_Center_X)),
                max(as.numeric(thissample$Location_Center_X))
            ),
            yrange = c(
                min(as.numeric(thissample$Location_Center_Y)),
                max(as.numeric(thissample$Location_Center_Y))
            ),
            marks = t(as.matrix(exprsMat))
        )
    
     value <-  spatstat.explore::nncorr(cell_points_cts)["correlation", ]
      value <-  value[  thisprotein]
     
    # Find the indices of the two nearest neighbors for each cell
    nn_indices <- nnwhich(cell_points_cts, k = 1)
    
    protein <-  thisprotein
    df <- data.frame(thiscell_exprs  = exprsMat[protein, ] , exprs =  exprsMat[protein,nn_indices ])
    
   p <-  ggplot(df, aes( x =thiscell_exprs ,  y = exprs , colour =  exprs  )) +
      geom_point(alpha = 0.3) + ggtitle(paste0( "Patient ", sample_name ,  " nn_corr = " ,  round(value, 2)  )) + scale_colour_viridis_c() + xlab("This cell expression") + ylab("Neighbouring cell expression")
   
   return (p ) 

}

    
p1 <- plot_nncorrelation("MB-0605",  "HER2")
p2 <- plot_nncorrelation("MB-0258",  "HER2")

ggarrange(plotlist = list(p1, p2))
```

![](NUS_workshop_v2_files/figure-html/unnamed-chunk-45-1.png)

### SessionInfo

Code

``` r
sessionInfo()
```

    R version 4.5.2 (2025-10-31)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.3 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

    time zone: UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] grid      stats4    stats     graphics  grDevices utils     datasets
    [8] methods   base

    other attached packages:
     [1] cowplot_1.2.0               SpatialExperiment_1.20.0
     [3] Banksy_1.6.0                cytomapper_1.22.0
     [5] EBImage_4.52.0              purrr_1.2.0
     [7] simpleSeg_1.12.0            naniar_1.1.0
     [9] reshape_0.8.10              scran_1.38.0
    [11] scater_1.38.0               scuttle_1.20.0
    [13] spatstat_3.4-1              spatstat.linnet_3.3-2
    [15] spatstat.model_3.4-2        rpart_4.1.24
    [17] spatstat.explore_3.5-3      nlme_3.1-168
    [19] spatstat.random_3.4-2       spatstat.geom_3.6-0
    [21] spatstat.univar_3.1-5       spatstat.data_3.1-9
    [23] survminer_0.5.1             ggpubr_0.6.2
    [25] tidyr_1.3.1                 scattermore_1.2
    [27] plotly_4.11.0               limma_3.66.0
    [29] dplyr_1.1.4                 spicyR_1.22.0
    [31] pheatmap_1.0.13             ggthemes_5.1.0
    [33] lisaClust_1.18.0            ClassifyR_3.14.0
    [35] survival_3.8-3              BiocParallel_1.44.0
    [37] MultiAssayExperiment_1.36.0 scFeatures_1.3.4
    [39] ggplot2_4.0.1               SingleCellExperiment_1.32.0
    [41] SummarizedExperiment_1.40.0 Biobase_2.70.0
    [43] GenomicRanges_1.62.0        Seqinfo_1.0.0
    [45] IRanges_2.44.0              S4Vectors_0.48.0
    [47] BiocGenerics_0.56.0         generics_0.1.4
    [49] MatrixGenerics_1.22.0       matrixStats_1.5.0

    loaded via a namespace (and not attached):
      [1] igraph_2.2.1               graph_1.88.0
      [3] Formula_1.2-5              BiocBaseUtils_1.12.0
      [5] tidyselect_1.2.1           bit_4.6.0
      [7] doParallel_1.0.17          clue_0.3-66
      [9] lattice_0.22-7             rjson_0.2.23
     [11] blob_1.2.4                 stringr_1.6.0
     [13] rngtools_1.5.2             S4Arrays_1.10.0
     [15] parallel_4.5.2             png_0.1-8
     [17] cli_3.6.5                  ggplotify_0.1.3
     [19] ProtGenerics_1.42.0        multtest_2.66.0
     [21] goftest_1.2-3              textshaping_1.0.4
     [23] BiocIO_1.20.0              bluster_1.20.0
     [25] grr_0.9.5                  BiocNeighbors_2.4.0
     [27] uwot_0.2.4                 curl_7.0.0
     [29] mime_0.13                  evaluate_1.0.5
     [31] tidytree_0.4.6             tiff_0.1-12
     [33] ComplexHeatmap_2.26.0      ggh4x_0.3.1
     [35] V8_8.0.1                   stringi_1.8.7
     [37] backports_1.5.0            lmerTest_3.1-3
     [39] XML_3.99-0.20              orthogene_1.16.0
     [41] httpuv_1.6.16              AnnotationDbi_1.72.0
     [43] magrittr_2.0.4             rappdirs_0.3.3
     [45] splines_4.5.2              mclust_6.1.2
     [47] KMsurv_0.1-6               jpeg_0.1-11
     [49] doRNG_1.8.6.2              DT_0.34.0
     [51] ggbeeswarm_0.7.2           DBI_1.2.3
     [53] terra_1.8-80               HDF5Array_1.38.0
     [55] jquerylib_0.1.4            withr_3.0.2
     [57] reformulas_0.4.2           class_7.3-23
     [59] matrixTests_0.2.3.1        systemfonts_1.3.1
     [61] ggnewscale_0.5.2           GSEABase_1.72.0
     [63] bdsmatrix_1.3-7            rtracklayer_1.70.0
     [65] htmlwidgets_1.6.4          fs_1.6.6
     [67] ggrepel_0.9.6              labeling_0.4.3
     [69] SparseArray_1.10.1         SingleCellSignalR_2.0.1
     [71] h5mread_1.2.0              annotate_1.88.0
     [73] zoo_1.8-14                 raster_3.6-32
     [75] XVector_0.50.0             knitr_1.50
     [77] nnls_1.6                   UCSC.utils_1.6.0
     [79] AUCell_1.32.0              foreach_1.5.2
     [81] dcanr_1.26.0               patchwork_1.3.2
     [83] data.table_1.17.8          ggtree_4.0.1
     [85] rhdf5_2.54.0               R.oo_1.27.1
     [87] ggiraph_0.9.2              irlba_2.3.5.1
     [89] gridGraphics_0.5-1         lazyeval_0.2.2
     [91] yaml_2.3.10                crayon_1.5.3
     [93] RColorBrewer_1.1-3         tweenr_2.0.3
     [95] later_1.4.4                codetools_0.2-20
     [97] GlobalOptions_0.1.2        KEGGREST_1.50.0
     [99] sccore_1.0.6               Rtsne_0.17
    [101] shape_1.4.6.1              gdtools_0.4.4
    [103] Rsamtools_2.26.0           filelock_1.0.3
    [105] leidenAlg_1.1.5            pkgconfig_2.0.3
    [107] GenomicAlignments_1.46.0   aplot_0.2.9
    [109] spatstat.sparse_3.1-0      ape_5.8-1
    [111] viridisLite_0.4.2          xtable_1.8-4
    [113] car_3.1-3                  plyr_1.8.9
    [115] httr_1.4.7                 rbibutils_2.4
    [117] tools_4.5.2                beeswarm_0.4.0
    [119] broom_1.0.10               dbplyr_2.5.1
    [121] crosstalk_1.2.2            survMisc_0.5.6
    [123] assertthat_0.2.1           lme4_1.1-37
    [125] digest_0.6.38              numDeriv_2016.8-1.1
    [127] Matrix_1.7-4               farver_2.1.2
    [129] AnnotationFilter_1.34.0    reshape2_1.4.5
    [131] yulab.utils_0.2.1          viridis_0.6.5
    [133] glue_1.8.0                 cachem_1.1.0
    [135] BiocFileCache_3.0.0        polyclip_1.10-7
    [137] proxyC_0.5.2               Biostrings_2.78.0
    [139] visdat_0.6.0               ggalluvial_0.12.5
    [141] statmod_1.5.1              concaveman_1.2.0
    [143] ScaledMatrix_1.18.0        fontBitstreamVera_0.1.1
    [145] carData_3.0-5              minqa_1.2.8
    [147] httr2_1.2.1                glmnet_4.1-10
    [149] dqrng_0.4.1                gtools_3.9.5
    [151] ggsignif_0.6.4             gridExtra_2.3
    [153] shiny_1.11.1               GSVA_2.4.1
    [155] BulkSignalR_1.2.1          R.utils_2.13.0
    [157] rhdf5filters_1.22.0        RCurl_1.98-1.17
    [159] memoise_2.0.1              rmarkdown_2.30
    [161] scales_1.4.0               R.methodsS3_1.8.2
    [163] stabledist_0.7-2           svglite_2.2.2
    [165] RANN_2.6.2                 fontLiberation_0.1.0
    [167] km.ci_0.5-6                EnsDb.Mmusculus.v79_2.99.0
    [169] cluster_2.1.8.1            msigdbr_25.1.1
    [171] spatstat.utils_3.2-0       coxme_2.2-22
    [173] scam_1.2-20                colorspace_2.1-2
    [175] rlang_1.1.6                EnsDb.Hsapiens.v79_2.99.0
    [177] GenomeInfoDb_1.46.0        DelayedMatrixStats_1.32.0
    [179] sparseMatrixStats_1.22.0   shinydashboard_0.7.3
    [181] aricode_1.0.3              ggforce_0.5.0
    [183] homologene_1.4.68.19.3.27  circlize_0.4.16
    [185] dbscan_1.2.3               mgcv_1.9-3
    [187] xfun_0.54                  iterators_1.0.14
    [189] abind_1.4-8                tibble_3.3.0
    [191] treeio_1.34.0              Rhdf5lib_1.32.0
    [193] bitops_1.0-9               Rdpack_2.6.4
    [195] fftwtools_0.9-11           promises_1.5.0
    [197] RSQLite_2.4.4              DelayedArray_0.36.0
    [199] compiler_4.5.2             boot_1.3-32
    [201] beachmat_2.26.0            RcppHungarian_0.3
    [203] Rcpp_1.1.0                 fontquiver_0.2.1
    [205] edgeR_4.8.0                BiocSingular_1.26.0
    [207] tensor_1.5.1               MASS_7.3-65
    [209] ggupset_0.4.1              babelgene_22.9
    [211] R6_2.6.1                   fastmap_1.2.0
    [213] rstatix_0.7.3              vipor_0.4.7
    [215] ensembldb_2.34.0           rsvd_1.0.5
    [217] gtable_0.3.6               deldir_2.0-4
    [219] htmltools_0.5.8.1          bit64_4.6.0-1
    [221] lifecycle_1.0.4            S7_0.2.1
    [223] nloptr_2.2.1               restfulr_0.0.16
    [225] sass_0.4.10                vctrs_0.6.5
    [227] ggfun_0.2.0                sp_2.2-0
    [229] bslib_0.9.0                gprofiler2_0.2.4
    [231] pillar_1.11.1              GenomicFeatures_1.62.0
    [233] magick_2.9.0               metapod_1.18.0
    [235] locfit_1.5-9.12            otel_0.2.0
    [237] jsonlite_2.0.0             svgPanZoom_0.3.4
    [239] cigarillo_1.0.0            GetoptLong_1.0.5          

### Acknowledgments

The authors thank all their colleagues, particularly at The University
of Sydney, Sydney Precision Data Science and Charles Perkins Centre for
their support and intellectual engagement. Special thanks to Ellis
Patrick, Shila Ghazanfar, Andy Tran, Helen, and Daniel for guiding and
supporting the building of this workshop.
