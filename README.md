# plasmid-model-cli

This tool runs the host-plasmid and plasmid-plasmid dynamics model in [Matlock and MacLean (2026)](https://doi.org/10.64898/2026.02.19.706745). If you use the tool, please cite:

```
@article{matlock2026conjugation,
  title={Conjugation structures plasmid populations through host-lineage restriction},
  author={Matlock, William and MacLean, R Craig},
  journal={bioRxiv},
  pages={2026--02},
  year={2026},
  publisher={Cold Spring Harbor Laboratory}
}
```

All you need is a host chromosomal tree (as a `.nwk`), and a presence/absence matrix of plasmid groups (formatted like `test_data/test_matrix.csv`). Please refer to the paper for interpretation of the outputs.

> Of course, you could use any tree and any binary traits!

## Quick start with Docker

From the repository root:

```bash
docker run --rm -it \
  -v "$PWD:/work" \
  wtmatlock/plasmid-model-cli:latest \
  --matrix /work/test_data/test_matrix.csv \
  --tree /work/test_data/test_tree.nwk \
  --out-dir /work/results \
  --k-dims 10 \
  --warmup 1000 \
  --samples 1000
```

The tool will provide information about your run:

```
Preparing data...

  Tips: 100
  Plasmids: 10
  Retained phylogenetic dimensions: 10
  Variance explained: 40.41%
  Phylogenetic dimensions needed to explain at least 99% of variance: 86

Running CmdStan sampling...

  Model: /app/stan/model
  Chains: 4
  Warmup: 1000
  Sampling: 1000
  Total iterations: 2000
  Seed: 123

  [00:00:33] ######################################## 8000/8000 (100%) sampling complete                                                                                                                                                                                      

Writing raw outputs...

Writing model summaries...

Completed successfully. Results saved to /work/results
```

The results are written to:

- `results/inputs/`, which contains the input data converted to JSON
- `results/raw/`, which contains the MCMC draws and other logs
- `results/summary/`, which contains the following summary tables:
  - `results/summary/plasmid_summary.csv` gives $k$ and $\sigma_k$ medians with 95% credible intervals
  - `results/summary/phi.csv` gives the signed tail probabilities for $\Phi$
  - `results/summary/r.csv` gives the signed tail probabilities for $R$
 
## Tips and best practice

- The model uses the host tree primarily as a measure of host relatedness. However, your input tree should be rooted, so it can be converted into a phylogenetic covariance matrix. In practice, this isn't always possible, and in my experience, it's fine to midpoint root.
- The number of phylogenetic dimensions retained (`--k-dims`) should balance faithful representation of the tree-derived covariance structure against computational cost. The tool reports the cumulative variance explained by the retained dimensions in the terminal. For publishable analyses, we recommend retaining enough dimensions to explain at least 99% of the variance (this value is also printed to the terminal).
- To that end, we also recommend `--warmup 5000 --samples 5000` as a minimum. You can inspect MCMC convergence using the files in `results/raw/`.

## Confessional 
Rust is a new (and exciting) language for me. I used GitHub Copilot to help cobble together an initial CLI, but after several unassisted iterations, I am confident it's doing the job. If you spot any weirdness, please let me know. My preference is to run Stan models in R, as I did for the paper.
