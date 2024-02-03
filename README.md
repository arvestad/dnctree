[![PyPI version](https://badge.fury.io/py/dnctree.svg)](https://badge.fury.io/py/dnctree)
# dnctree: Randomized divide and conquer algorithm for phylogenetic trees

This is a distance-based method, inspired by Neighbor-Joining, for inferring phylogenies.
Its main feature is that it scales to very large input datasets by avoiding estimating
a pairwise distance matrix. 

The input is a multiple sequence alignment and the output is a tree in standard Newick
format. Option `--json-output` wraps the Newick-formatted tree in a JSON file, with some
additional data describing the computation.

The implementation is (currently) 100% pure Python, but you can handle very large datasets
in reasonable time anyways. See the [preprint](https://doi.org/10.1101/2023.10.11.561902)
for some examples!

## Algorithms

There are currently two algorithm versions implemented in dnctree. 

- The default algorithm, 'core tree', seems to be as accurate as Neighbor-Joining in our experiments
  so far, but scales much better than Neighbor-Joining. We have submitted a paper on this algorithm.
- The 'simple' algorithms (see option `--simple`) is much faster, but has much worse accuracy than
  the core tree algorithm. The simple algorithm is described in a
  [biorXiv](https://doi.org/10.1101/2023.10.11.561902) preprint.

Both algorithms use Divide-and-Conquer. Problem instances with fewer sequences than
what is given by the "base-case size", 100 by default, are handled by Neighbor-Joining
and larger instances are partitioned and handled recursively. 


## Input formats

The formats Fasta, Phylip, Clustal, Nexus, and Stockholm (Pfam) are currently accepted
input formats. This is determined by what the BioPython package accepts.


## Example usage

```shell
dnctree testdata/s83_L500.phylip
dnctree -f phylip testdata/s83_L500.phylip   # Making it very clear input is a Phylip file
dnctree --simple  testdata/s83_L500.phylip   # Using the faster "simple" algorithm
dnctree --base-case-size 10 testdata/s83_L500.phylip  # Divide and conquer on larger inputs
```

Examples with output:
```shell
$ dnctree testdata/s83_L500.phylip
((((((L26,L27),(L24,L25)),((L28,L29),(L30,L31))),(((L16,L17),(L18,L19)),((L22,L23),(L20,L21)))),((((L8,L9),(L10,L11)),((L12,L13),(L14,L15))),(((L2,L3),(L0,L1)),((L6,L7),(L4,L5))))),(((((L86,L87),(L84,L85)),((L80,L81),(L82,L83))),(((L94,L95),(L92,L93)),((L90,L91),(L88,L89)))),((((L74,L75),(L72,L73)),((L78,L79),(L76,L77))),(((L64,L65),(L66,L67)),((L68,L69),(L70,L71))))),(((((L44,L45),(L46,L47)),((L42,L43),(L40,L41))),(((L34,L35),(L32,L33)),((L38,L39),(L36,L37)))),((((L50,L51),(L48,L49)),((L54,L55),(L52,L53))),(((L58,L59),(L56,L57)),((L62,L63),(L60,L61))))));

$ dnctree --json-output --base-case-size 10 testdata/s83_L500.phylip
{
    "version": "dnctree 1.0",
    "tree": "((((((L45,L44),(L47,L46)),((L43,L42),(L40,L41))),((((L54,L55),(L52,L53)),((L50,L51),(L48,L49))),(((L58,L59),(L57,L56)),((L63,L62),(L61,L60))))),(((L35,L34),(L33,L32)),((L39,L38),(L36,L37)))),(((((L91,L90),(L89,L88)),((L94,L95),(L92,L93))),(((L81,L80),(L82,L83)),((L87,L86),(L84,L85)))),((((L78,L79),(L77,L76)),(((L69,L68),(L70,L71)),((L65,L64),(L66,L67)))),((L74,L75),(L72,L73)))),(((((L7,L6),(L4,L5)),((L3,L2),(L0,L1))),((((L21,L20),(L23,L22)),((L16,L17),(L19,L18))),(((L30,L31),(L29,L28)),((L26,L27),(L24,L25))))),(((L9,L8),(L10,L11)),((L12,L13),(L14,L15)))));",
    "infile": "testdata/s83_L500.phylip",
    "aligned": true,
    "base-case-size": 10,
    "distances-computed": 1711,
    "fraction-computed-distances": 0.375,
    "n-taxa": 96,
    "comment": "Computed 1711 distances for 96 taxa. A full distance matrix would contain 4560 pairs. Savings: 62.5 %",
    "model-name": "WAG",
    "description": "AA alignment",
    "msa-width": 500,
    "computing-time": 0.907991542
}
```

# Credits

* Amy Lee Jalsenius developed and implemented the "core tree" algorithm which is now the default.
* Mazen Mardini added PaHMM code (see `https://github.com/marbogusz/paHMM-Tree`), enabling experiments.