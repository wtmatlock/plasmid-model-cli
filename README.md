# plasmid-model-cli

This tool runs the plasmid-host dynamics model in [Matlock and MacLean (2026)](https://doi.org/10.64898/2026.02.19.706745). If you use the tool, please cite:

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

## Quick start with Docker

From repo root:

```bash
docker build -t plasmid-model-cli --build-arg STAN_MODEL=stan/model.stan .
```

This will compile:

- CmdStan
- Stan model executable (`stan/model.stan` to `stan/model`)
- Rust CLI binary

You can then test with the bundled data:

```bash
docker run --rm -it \
  -v "$PWD:/work" \
  plasmid-model-cli \
  --matrix /work/test_data/test_matrix.csv \
  --tree /work/test_data/test_tree.nwk \
  --out-dir /work/results \
  --k-dims 10 \
  --warmup 1000 \
  --samples 1000
```

The results are written to:

- `results/inputs/`, which contains the input data converted to JSON
- `results/raw/`, which contains the MCMC draws and other logs
- `results/summary/`, which contains the following summary tables:
  - `results/summary/plasmid_summary.csv` gives $k$ and $\sigma_k$ medians with 95% credible intervals
  - `results/summary/phi.csv` gives the signed tail probabilities for $\Phi$
  - `results/summary/r.csv` gives the signed tail probabilities for $R$
 
## Best practices

- The number of phylogenetic dimensions retained (--k-dims) should balance faithful representation of the tree-derived covariance structure against computational cost. The tool reports the cumulative variance explained by the retained dimensions in the terminal. For publishable analyses, we recommend retaining enough dimensions to explain at least 99% of the variance.
- To that end, we also recommend `--warmup 5000 --samples 5000` as a minimum. You can inspect MCMC convergence using the files in `results/raw/`.
