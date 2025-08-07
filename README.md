# moltenNMF

`moltenNMF` is an R package for Nonnegative Matrix Factorization for tidy data format.

## Reference

Ko ABE and Teppei SHIMAMURA (2023) UNMF: A unified non-negative matrix factorization for multi-dimensional omics data. Briefings in Bioinformatics. [https://doi.org/10.1093/bib/bbad253](https://doi.org/10.1093/bib/bbad253)

## Installation

```R
devtools::install_github("abikoushi/moltenNMF")
```

## Example

1. [A simple example of 'moltenNMF'](https://github.com/abikoushi/moltenNMF/blob/main/Rmd/moltenNMF.md)
2. [missing data imputation](https://github.com/abikoushi/moltenNMF/blob/main/Rmd/gehag.md)


## Note

Version 2 will be released soon.

It will be possible to handle variational Bayesian methods in a way that is suitable for sparse matrix format.

Reference : 
Ko Abe & Teppei Shimamura “Scalable Procedures for Nonnegative Matrix Factorization with Skip-Zeros and its Extension” (in preparation)

New Features : 

- `NMF2D_vb`
- `NMF2D_svb`
- `mNMF_svb`
