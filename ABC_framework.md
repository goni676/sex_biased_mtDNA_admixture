## ABC framework

Generate 1000 parameter vectors, denoted by $\theta_i$. Each parameter in $\theta_i$ is sampled from a uniform distribution over its range. Thus:

$$p(\theta_i) = \prod_{x \in \theta_i} p(x)$$

In ABC terminology, this product defines the **prior** distribution of parameter vector $i$. Note that $p(\theta_i) = p(\theta_j)$ for every $i,j \in \[1000\]$.

Under the same framework, the **likelihood** expresses the probability of observing approximately 10% of EUR_mtDNA. We need to define an acceptance range around 10%, for example 8%–12%, and then $D = \{0.08 \le EUR_{mtDNA} \le 0.12\}$, and the likelihood is $p(D|\theta_i)$. From an implementation standpoint, we simulate the process for 100 independent repetitions and calculate the proportion of repetitions that achieve this goal:

$$p(D|\theta_i) = \frac{1}{100} \sum_{k=1}^{100} \text{INDICATOR}(mtDNA_{ik} \in [0.08, 0.12])$$

For the **posterior**, we use the rule:

$$p(\theta_i|D) = \frac{p(D|\theta_i) \cdot p(\theta_i)}{p(D)}$$

$$p(D) = \sum_{i=1} p(D|\theta_i) \cdot p(\theta_i) = p(\theta_i) \sum_{i} p(D \mid \theta_i)$$


Therefore:

$$p(\theta_i \mid D) = \frac{p(D \mid \theta_i)}{\sum_{i} p(D \mid \theta_i)}$$
