# Covering design records

Explicit covering designs that improve the best known size listed by the [Covering Repository](https://www.coveringrepository.com). A (v,k,t) covering design is a family of k-subsets (blocks) of a v-set such that every t-subset lies in at least one block.

| Parameters (v,k,t) | Previous best | New size | File |
|---|---:|---:|---|
| (89,8,3) | 2611 (Bluskov and Sidorenko, product construction, 22 February 2026) | **2606** | `designs/C_89_8_3_2606.txt` |

Each file lists one block per line, on the points 1 to v.

## Method

Post-optimisation of the published design. `generators/cd_postop.c` removes the block that covers the fewest triples alone, then repairs the uncovered triples by simulated annealing, and repeats while it succeeds. The program makes no call to a language model at run time.

    cc -O3 -o cd_postop generators/cd_postop.c -lm
    ./cd_postop 89 8 3 300 1 < start.txt > out.txt

The starting blocks are the 2611-block design published in the repository. The run is randomised and time-limited, so a rerun gives a valid design of similar size, possibly a different one.

## Verification

    python3 verify.py designs/C_89_8_3_2606.txt 89 8 3

The script checks each block and enumerates the triples they cover: all 113564 are covered.

## Prior art

Searched on 31 October 2026: the Covering Repository, the most recent arXiv papers and Zenodo records matching "covering design". Nothing found on these parameters. A smaller unpublished design remains possible.

## Use of a language model

The programs were generated with Claude (Anthropic), which was also used as a writing tool for this text.

## Licence

The designs are under CC BY 4.0, the code under the MIT licence.
