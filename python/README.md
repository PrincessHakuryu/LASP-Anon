# AMM°‰Q3 Reference Structures and NNP®CShengBTE Workflow

This folder contains the ten reference AMM°‰Q3 structures used in the paper.

## Required environment

The Python scripts were tested with the following environment:

```text
ase       3.28.0
numpy     2.4.6
PyYAML    6.0.3
phonopy   3.5.1
torch     2.7.0+cu126
jax       0.10.1
lxml      6.1.0
```

The `thirdorder` package is also required. Install it from its GitHub repository and provide the path to `thirdorder_vasp.py` when running the script.

## Optional calculator packages

Install these packages only if they are selected or configured in `nnp-kappa.py`:

```text
nequip
nequix
mace
```

## Input structure

Each structure should contain a `POSCAR` file:

```text
article-10/
©¿©§©§ BaCuScSe3/
©¶   ©∏©§©§ POSCAR
©¿©§©§ BaCuScTe3/
©¶   ©∏©§©§ POSCAR
©¿©§©§ CsZnYS3/
©¶   ©∏©§©§ POSCAR
©∏©§©§ ...
```

## Configuring the neural-network potential

Edit the calculator configuration in `nnp-kappa.py`. For example, for Nequix:

```python
CALCULATOR_KIND = "nequix"
NEQUIX_MODEL_NAME = "nequix-oam-1-pft"
NEQUIX_BACKEND = "jax"
NEQUIX_USE_KERNEL = True
```

If you use a different calculator, modify the `make_calc()` function accordingly.

## Running the workflow

For a single structure:

```bash
python nnp-kappa.py \
    --poscar BaCuScSe3/POSCAR \
    --workdir BaCuScSe3-BTE \
    --no-relax \
    --no-relax-cell \
    --fc2-dim 2 2 1 \ll
    
    --fc3-dim 2 2 1 \
    --fc3-cutoff -3 \
    --qgrid 12 12 12 \
    --thirdorder /path/to/thirdorder_vasp.py \
    --overwrite
```

The options `--no-relax` and `--no-relax-cell` disable ionic and cell optimization, respectively. The existing `POSCAR` is used directly.

## Generated results

The output directory will contain:

```text
article-10/
©¿©§©§ BaCuScSe3-BTE/
©¶   ©¿©§©§ POSCAR
©¶   ©¿©§©§ CONTROL
©¶   ©¿©§©§ FORCE_CONSTANTS_2ND
©¶   ©¿©§©§ FORCE_CONSTANTS_3RD
©¶   ©∏©§©§ NNP_BTE_REPORT.yaml
©∏©§©§ ...
```

The script generates the second- and third-order force constants and a `CONTROL` file for ShengBTE. ShengBTE must be run separately after this workflow.
