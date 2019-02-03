snATAC-seq of mouse mammary cells
================
Jay Chung
2/2/2019

## Feature (regions) selection

### Read files:

``` r
library(GenomicRanges)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colMeans, colnames, colSums, dirname, do.call, duplicated,
    ##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
    ##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
    ##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
    ##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

``` r
library(IRanges)
peaks.df <- read.table("Mammary_fetal+adult_merged_macs2_summits.bed") # these are regions called by macs2 from aggregate snATAC-seq
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
peaks.rs.gr <- resize(reduce(resize(peaks.gr, 1000, fix = "center")), 1000, fix = "center") # 164657
all(width(peaks.rs.gr) == 1000)
```

    ## [1] TRUE

``` r
peaks.rs.df <- as.data.frame(peaks.rs.gr)[,1:3]
```

### Separate promoters from enhancers:

``` r
proms.df <- read.table("mm10_protein_coding_genes_3kbTSS.bed")
proms.gr <- GRanges(seqnames = proms.df[,1], 
                    ranges = IRanges(proms.df[,2], proms.df[,3]))
peaks.np.gr <- peaks.rs.gr[-queryHits(findOverlaps(peaks.rs.gr, proms.gr))] # 145453
peaks.p.gr <- peaks.rs.gr[queryHits(findOverlaps(peaks.rs.gr, proms.gr))] # 21179
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
write.table(peaks.np.ex.df, "Mammary_v1+2_fetal+adult_distal.ygi", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(peaks.p.ex.df, "Mammary_v1+2_fetal+adult_promoter.ygi", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
```
