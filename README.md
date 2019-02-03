snATAC-seq of mouse mammary cells
================
Jay Chung |
Feb 2019

Raw files can be downloaded here:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125523>

## Feature (regions) selection

``` r
library(GenomicRanges)
library(IRanges)
```

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
    ## [1] GenomicRanges_1.32.7 GenomeInfoDb_1.16.0  IRanges_2.14.12     
    ## [4] S4Vectors_0.18.3     BiocGenerics_0.26.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.0             digest_0.6.18          bitops_1.0-6          
    ##  [4] magrittr_1.5           evaluate_0.12          zlibbioc_1.26.0       
    ##  [7] stringi_1.2.4          XVector_0.20.0         rmarkdown_1.11        
    ## [10] tools_3.5.2            stringr_1.3.1          RCurl_1.95-4.11       
    ## [13] yaml_2.2.0             compiler_3.5.2         htmltools_0.3.6       
    ## [16] knitr_1.20             GenomeInfoDbData_1.1.0
