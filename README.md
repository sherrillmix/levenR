# Calculate Levenshtein distance between strings

## Introduction

`violinPointR` provides a few functions for simple Levenshtein alignment and distance calculation with multiple threads, ends-free and reduced homopolymer gap costs.

## Installation


```r
devtools::install_github("sherrillmix/levenR")
```

## Examples

### Generate Levenshtein distance matrix

An example of calculating the Levenshtein distance between several strings to make a distance matrix:

```r
library(levenR)
seqs <- c("AAATA", "AATA", "AAAT", "ACCTA")
leven(seqs)
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    0    1    1    2
## [2,]    1    0    2    2
## [3,]    1    2    0    3
## [4,]    2    2    3    0
```

### Compare many to one

An example of calculating the Levenshtein distance between several strings against a longer reference sequence:

```r
library(levenR)
seqs <- c("AAATA", "AATA", "AAAT", "ACCTA")
ref <- "CCAAATACCGACC"
leven(seqs, ref, substring2 = TRUE)
```

```
##      [,1]
## [1,]    0
## [2,]    0
## [3,]    0
## [4,]    1
```

### Find the best reference

An example of calculating the Levenshtein distance between several strings against two longer reference sequences and determining the best match for each read:

```r
library(levenR)
seqs <- c("AAATA", "AATA", "AAAT", "ACCTA")
refs <- c("CCATAATACCGACC", "GGAAATACCTA")
dist <- leven(seqs, refs, substring2 = TRUE)
apply(dist, 1, which.min)
```

```
## [1] 2 1 2 2
```

### Alignment

An example of aligning strings against a longer reference:

```r
library(levenR)
seqs <- c("AAATA", "AATA", "AAAT", "ACCTA")
ref <- "CCAAATACCGACC"
levenAlign(seqs, ref, substring2 = TRUE)
```

```
## $ref
## [1] "CCAAATACCGACC"
## 
## $align
## [1] "--AAATA------" "---AATA------" "--AAAT-------" "------ACCTA--"
```



