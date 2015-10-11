# Calculate Levenshtein distance between strings

## Introduction

`violinPointR` provides a few functions for simple Levenshtein alignment and distance calculation with multiple threads, ends-free and reduced homopolymer gap costs.

## Installation


```r
devtools::install_github("sherrillmix/levenR")
```

## Examples

### Basic Levenshtein distance

We use the provided function `offsetX` to generate the x-offsets for plotting.

```r
library(levenR)
seqs<-c('AAATA','AATA','AAAT','ACCTA')
leven(seqs)
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    0    1    1    2
## [2,]    1    0    2    2
## [3,]    1    2    0    3
## [4,]    2    2    3    0
```

