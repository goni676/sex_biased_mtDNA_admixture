# Sex-Biased mtDNA Admixture

This project compares theoretical mathematical predictions and msprime simulations for the probability of observing extreme mtDNA ancestry proportions after an admixture event.

## Guiding Questions

### Scenario 1 — No Initial Bias
>
>**Given a population of $N$ mtDNA haplotypes, 50:50 from two populations $A$ and $B$, what are the chances of observing 80% $A$-derived mtDNA after $X$ generations?**
>
#### Theoretical model

let $K_t$ denote the number of $A$-derived mtDNA haplotypes after $t$ generations, and let $p_t = K_t/N$ denote their frequency.

The drift process is modeled as:

$$
K_t \sim \mathrm{Bin}\left(N, \frac{K_{t-1}}{N}\right)
$$

The probability of interest is:

$$
P\left(K_X \geq 0.8N \mid K_0 = \frac{N}{2}\right)
$$


### Scenario 2 - Initial mtDNA Bias
>
>**Given a population of $N$ mtDNA haplotypes originated from populations $A$ and $B$, where the initial mitochondrial contribution is biased toward population $A$, what are the chances of observing 80% A-derived mtDNA after $X$ generations?**
>


Let $K_t$ denote the number of A-derived mtDNA haplotypes after $t$ generations.

Let $M_A$ denote the proportion of mothers from population A.

$$
K_0 = M_A \cdot N
$$

Then, in each generation $t$:

$$
K_t \sim \mathrm{Bin}\left(N,\frac{K_{t-1}}{N}\right)
$$

The required probability is therefore:

$$
P\left(K_X \ge 0.8N \\middle|\ K_0 = M_A N\right)
$$

### msprime simulations

A demographic model is constructed using **msprime** with four populations: two source populations `A` and `B`, an admixed population `ADMIX`, and a common ancestral population `ANC`.

The admixture event is modeled by assigning each mtDNA lineage from ADMIX to either `A` or `B` at the admixture time, which corresponds to the parameter $X$ in the mathematical formulation. Setting `record_migrations=True` makes it possible to recover, for each sampled mtDNA lineage, whether its ancestry traces back to population `A` or population `B`.

The simulations are repeated 1000 times with different random seeds, and I estimate the probability of observing at least 80% A-derived mtDNA.

## Repository structure

```text
project/
│
├── notebooks/
│   ├── 00_exploring_msprime.ipynb
│   ├── 01_drift_prob_results.ipynb
│   └── 02_sex_bias_results.ipynb
│
├── src/
│   ├── __init__.py
│   ├── mtdna_msprime.py
│   └── mtdna_theory.py
│
└── results/
```

### notebooks/

#### `00_exploring_msprime.ipynb`
Exploring the demographic model, the admixture event, and the msprime data structures used to document ancestry.

#### `01_drift_prob_results.ipynb`
Comparing the theoretical model and the msprime simulations for **Scenario 1**.

#### `02_sex_bias_results.ipynb`
Comparing the theoretical model and the msprime simulations for **Scenario 2**.

### src/

#### `mtdna_msprime.py`

- demographic model construction using msprime
- ancestry simulation using msprime
- estimation of probabilities from repeated simulations

#### `mtdna_theory.py`

- implementation of the theoretical drift model
- calculation of theoretical probabilities
