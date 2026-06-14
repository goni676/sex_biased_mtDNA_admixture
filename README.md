# Sex-Biased mtDNA Admixture

## Scenario 1 — No Initial Bias

The guiding question is:

**Given a population of $N$ mtDNA haplotypes, 50:50 from two populations $A$ and $B$, what are the chances of observing 80% $A$-derived mtDNA after $X$ generations?**

Mathematically, let $K_t$ denote the number of $A$-derived mtDNA haplotypes after $t$ generations, and let $p_t = K_t/N$ denote their frequency.

The drift process is modeled as:

$$
K_t \sim \mathrm{Bin}\left(N, \frac{K_{t-1}}{N}\right)
$$

The probability of interest is:

$$
P\left(K_X \geq 0.8N \mid K_0 = \frac{N}{2}\right)
$$


### Scenario 2 - Initial mtDNA Bias

**Given a population of N mtDNA haplotypes originated from populations A and B, where the initial mitochondrial contribution is biased toward population A, what are the chances of observing 80% A-derived mtDNA after X generations?**

Let $K_t$ denote the number of A-derived mtDNA haplotypes after $t$ generations.

Let $M_A$ denote the proportion of mothers from population A.

$$
K_0 = M_A \cdot N
$$

Then, in each generation $t$,

$$
K_t \sim \mathrm{Bin}\left(N,\frac{K_{t-1}}{N}\right)
$$

The required probability is therefore

$$
P\left(K_X \ge 0.8N \\middle|\ K_0 = M_A N\right)
$$

To my understanding, if we model sex bias alone, without drift, then we would not randomly sample $N$ haplotypes in each generation. Instead, the mtDNA proportion would remain constant, so that for every generation $t$,

$$
p_t = M_A.
$$
