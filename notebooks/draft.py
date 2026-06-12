import numpy as np
import matplotlib.pyplot as plt

# ----- Your assumptions / parameters -----

N = 100                 # number of mtDNA haplotypes
M_A = 0.6               # initial maternal contribution from population A
generations = 100       # number of generations to simulate
num_simulations = 1000  # number of repeated simulations
threshold = 0.8         # 80% A-derived mtDNA

# Initial number of A-derived mtDNA haplotypes
K0 = int(M_A * N)

# p_values[simulation, generation]
p_values = np.zeros((num_simulations, generations + 1))

# ----- Simulation -----

for sim in range(num_simulations):
    K = K0
    p_values[sim, 0] = K / N

    for t in range(1, generations + 1):
        p = K / N
        
        # drift: sample N mtDNA haplotypes from previous generation
        K = np.random.binomial(N, p)
        
        p_values[sim, t] = K / N

# ----- Probability graph -----

prob_above_threshold = np.mean(p_values >= threshold, axis=0)


plt.figure(figsize=(8, 5))

plt.plot(range(generations + 1), prob_above_threshold)

plt.xlabel("Generation")
plt.ylabel("P(A-derived mtDNA >= 80%)")
plt.title("Probability of reaching 80% A-derived mtDNA over generations")
plt.ylim(0, 1)

plt.show()