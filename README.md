[![PyPI version](https://badge.fury.io/py/dnctree.svg)](https://badge.fury.io/py/dnctree)
# dnctree: Randomized divide and conquer algorithm for phylogenetic trees

See the [biorXiv](https://doi.org/10.1101/2023.10.11.561902) paper about the method.


## Example usage

``` shell
dnctree -f phylip sample_data/L500_sd0.0_WAG_n0.phylip
dnctree --first_triple H32 H14 H0  --verbose --base_case_size 12 -f phylip sample_data/L500_sd0.0_WAG_n0.phylip
```
