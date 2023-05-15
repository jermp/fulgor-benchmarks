Benchmarks for the Fulgor index
===============================

This repository contains benchmarking code and evaluation scripts for the [Fulgor](https://github.com/jermp/fulgor) index.

  - Snakefile - Rules for generating the simulated reads and evaluating pseudo-alignment output.
  - sim_reads.sh - Script using mason to generate simulated reads.
  - evaluate_mappings.py - Given the pseudo-alignment mappings output by Fulgor, generate the json file that contain the different mapping statistics.
  - map_contigs.sh - Script mapping references to ids.
  
**More scripts to be populated soon.**

**Packages** <br />
Mason - 2.0.9 <br />
Seqkit - 2.4.0
