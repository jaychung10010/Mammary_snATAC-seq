snATAC-seq of mouse mammary cells
================
Jay Chung |
Feb 2019

Raw files can be downloaded here:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125523>

``` r
library(GenomicRanges)
library(IRanges)
library(Matrix)
library(irlba)
library(Rtsne)
library(scales)
```

## Feature (regions) selection

### Read files:

``` r
peaks.df <- read.table("Mammary_fetal+adult_merged_macs2_summits.bed") # these are regions called by macs2 from aggregate snATAC-seq
nrow(peaks.df)
```

    ## [1] 437903

``` r
head(peaks.df)
```

    ##     V1      V2      V3
    ## 1 chr1 3100660 3100661
    ## 2 chr1 3203924 3203925
    ## 3 chr1 3203940 3203941
    ## 4 chr1 3212845 3212846
    ## 5 chr1 3212897 3212898
    ## 6 chr1 3212923 3212924

### Resize regions to 1 kb:

``` r
peaks.gr <- GRanges(seqnames = peaks.df[,1], 
                    ranges = IRanges(peaks.df[,2], peaks.df[,3]))
peaks.rs.gr <- resize(reduce(resize(peaks.gr, 1000, fix = "center")), 1000, fix = "center")
all(width(peaks.rs.gr) == 1000)
```

    ## [1] TRUE

``` r
peaks.rs.df <- as.data.frame(peaks.rs.gr)[,1:3]
nrow(peaks.rs.df)
```

    ## [1] 164657

### Separate promoters from enhancers:

``` r
proms.df <- read.table("mm10_protein_coding_genes_3kbTSS.bed")
proms.gr <- GRanges(seqnames = proms.df[,1], 
                    ranges = IRanges(proms.df[,2], proms.df[,3]))
peaks.np.gr <- peaks.rs.gr[-queryHits(findOverlaps(peaks.rs.gr, proms.gr))]
peaks.p.gr <- peaks.rs.gr[queryHits(findOverlaps(peaks.rs.gr, proms.gr))]
all(width(peaks.np.gr) == 1000)
```

    ## [1] TRUE

``` r
all(width(peaks.p.gr) == 1000)
```

    ## [1] TRUE

``` r
peaks.np.ex.df <- as.data.frame(peaks.np.gr)[,1:3]
peaks.p.ex.df <- as.data.frame(peaks.p.gr)[,1:3]
nrow(peaks.np.ex.df)
```

    ## [1] 145453

``` r
nrow(peaks.p.ex.df)
```

    ## [1] 21179

``` r
write.table(peaks.np.ex.df, "Mammary_v1+2_fetal+adult_distal.ygi", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(peaks.p.ex.df, "Mammary_v1+2_fetal+adult_promoter.ygi", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
```

## Data transformation and dimension reduction

``` r
load("./fMaSC+adult_v1+2_rep1+2.distal.sparse_mat.RData") # load sparse bmat
dim(all_distal_bmat_sm)
```

    ## [1]   7846 145453

### Remove bottom 10% usage regions

``` r
dcount <- as.vector(Matrix::colSums(all_distal_bmat_sm))
didx_keep <- which(dcount > quantile(dcount, probs = 0.1))
all_distal_bmat_sm <- all_distal_bmat_sm[,didx_keep]
all_distal_bmat_sm_t <- t(all_distal_bmat_sm) # transpose to row:regions, column:cells
```

### LSI of distal regions

``` r
nfreqs <- t(t(all_distal_bmat_sm_t) / Matrix::colSums(all_distal_bmat_sm_t))
idf <- as(log(1 + ncol(all_distal_bmat_sm_t) / Matrix::rowSums(all_distal_bmat_sm_t)), "sparseVector")
tf_idf_counts_distal <- as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs
dim(tf_idf_counts_distal)
```

    ## [1] 130899   7846

``` r
set.seed(524)
bmat_svd <- irlba(tf_idf_counts_distal, 50, 50)
d_diagtsne <- matrix(0, 50, 50)
diag(d_diagtsne) <- bmat_svd$d
bmat_svd_vd <- t(d_diagtsne %*% t(bmat_svd$v))
dim(bmat_svd_vd)
```

    ## [1] 7846   50

### t-SNE of distal regions

Tune t-SNE by selecting lowest KL divergence

    tsne_svd <- vector("list", 9)
    KL <- numeric()
    par(mfrow=c(3,3))
    set.seed(524)
    for (i in 1:9){
      mat <- Rtsne(bmat_svd_vd, dims = 2, perplexity = 30, verbose = TRUE, 
                   max_iter = 1000, check_duplicates = FALSE, is_distance = FALSE, 
                   theta = 0.5, pca = FALSE, exaggeration_factor = 12)
      tsne_svd[[i]] <- mat$Y
      plot(tsne_svd[[i]], cex = 0.5, pch = 16, xlab = "t-SNE1", ylab = "t-SNE2", 
           main=paste0(i, "-TFIDF SVD t-SNE v1+2 distal"), col = alpha(col, 0.3))
      KL[i] <- round(min(mat$itercosts), 4)
      legend("topleft", paste0("KL = ", KL[i]), text.col = "black", bty = "n")
    }
    best_tsne <- tsne_svd[[which.min(KL)]] # lowest KL
    write.csv(best_tsne, "./Bmat_tfidf_svd_tsne_best_v1+2_distal.csv")

``` r
best_tsne <- read.csv("./Bmat_tfidf_svd_tsne_best_v1+2_distal.csv", row.names = 1)
col <- c(rep("darkgreen", 2577), rep("orange", 5269))
plot(best_tsne, cex = 0.8, pch = 16, xlab = "t-SNE1", ylab = "t-SNE2", 
     main = "t-SNE of snATAC-seq (7846 cells)", col = alpha(col, 0.5))
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Mojave 10.14
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] scales_1.0.0         Rtsne_0.13           irlba_2.3.2         
    ## [4] Matrix_1.2-15        GenomicRanges_1.32.7 GenomeInfoDb_1.16.0 
    ## [7] IRanges_2.14.12      S4Vectors_0.18.3     BiocGenerics_0.26.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.0             knitr_1.20             XVector_0.20.0        
    ##  [4] magrittr_1.5           zlibbioc_1.26.0        munsell_0.5.0         
    ##  [7] colorspace_1.3-2       lattice_0.20-38        stringr_1.3.1         
    ## [10] tools_3.5.2            grid_3.5.2             htmltools_0.3.6       
    ## [13] yaml_2.2.0             digest_0.6.18          GenomeInfoDbData_1.1.0
    ## [16] bitops_1.0-6           RCurl_1.95-4.11        evaluate_0.12         
    ## [19] rmarkdown_1.11         stringi_1.2.4          compiler_3.5.2
