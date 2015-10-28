# LevenR --- Calculate Levenshtein distance between strings

[![Build Status](https://travis-ci.org/sherrillmix/levenR.svg?branch=master)](https://travis-ci.org/sherrillmix/levenR)

## Introduction

`violinPointR` provides a few functions for simple Levenshtein alignment and distance calculation with multiple threads, ends-free and reduced homopolymer gap costs.

## Installation
To install directly from github, use the [<code>devtools</code>](https://github.com/hadley/devtools) library and run:

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

### Ignoring indels in homopolymers
An example of calculating the Levenshtein distance between several strings to make a distance matrix while ignoring indels in long homopolymers (an error type common in 454 and IonTorrent sequencing):

```r
library(levenR)
seqs <- c("AAAAATA", "AAATTTTTA", "AAAAATTTA")
leven(seqs, homoLimit = 3)
```

```
##      [,1] [,2] [,3]
## [1,]    0    2    2
## [2,]    2    0    0
## [3,]    2    0    0
```

### Using multiple threads 
An example of calculating the Levenshtein distance between several strings using multiple threads:

```r
library(levenR)
seqs <- replicate(50, paste(sample(letters, 100, TRUE), collapse = ""))
system.time(leven(seqs))
```

```
##    user  system elapsed 
##   0.566   0.090   0.645
```

```r
system.time(leven(seqs, nThreads = 4))
```

```
##    user  system elapsed 
##   0.159   0.000   0.061
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



