import msprime
import numpy as np
import pandas as pd

"""
========================================
Demographic model and ancestry tracing
========================================
Defines the demographic model and runs an mtDNA ancestry simulation.
'count_origin' counts how many samples originate from a specified population.

"""

def define_demography_AFR_EUR(theta, seed):
    N_ANC = theta["N_ANC"]
    N_AFR = theta["N_AFR"]
    N_EUR = theta["N_EUR"]
    N_ADMIX = theta["N_ADMIX"]
    ADMIX_TIME = theta["ADMIX_TIME"]
    SPLIT_TIME = theta["SPLIT_TIME"]
    P_AFR = theta["p_AFR"]

    demography = msprime.Demography()
    demography.add_population(name="AFR", initial_size=N_AFR)
    demography.add_population(name="EUR", initial_size=N_EUR)
    demography.add_population(name="ADMIX", initial_size=N_ADMIX)
    demography.add_population(name="ANC", initial_size=N_ANC)
    
    demography.add_population_split(time=SPLIT_TIME, derived=["AFR", "EUR"], ancestral="ANC")

    demography.add_admixture(time=ADMIX_TIME, derived="ADMIX", ancestral=["AFR", "EUR"], proportions=[P_AFR, 1 - P_AFR]) 
    
    demography.sort_events()
    
    ts = msprime.sim_ancestry(
        samples={"ADMIX": 1000},
        demography=demography,
        ploidy=1,
        sequence_length=16569,
        recombination_rate=0,
        random_seed=seed,
        record_migrations=True,
    )
    return ts


def count_origin(curr_origin, ts, ADMIX_TIME):
    cnt = 0
    tree = ts.first()
    for sample in ts.samples():
        origin = None
        for m in ts.migrations():
            if m.time == ADMIX_TIME:
                if tree.is_descendant(sample, m.node):
                    origin = ts.population(m.dest).metadata["name"]
                    if origin == curr_origin:
                        cnt += 1
                    break
    
    return cnt


"""
=================
Theta creators
=================
Functions for sampling a theta parameter vector from the specified parameter ranges.
The sampling method is selected according to prior_type (uniform/normal)

"""

def create_theta_normal(parameter_vector):
    rng = np.random.default_rng()
    theta = {}

    for parameter, values in parameter_vector.items():
        if parameter == "SPLIT_TIME":
            theta[parameter] = values
            continue
        lower, mean, upper = values
        sd = (upper-lower) / (2 * 1.96)
        sampled_value = rng.normal(
            loc=mean,
            scale=sd,
        )

        if parameter == "p_AFR":
            theta[parameter] = float(sampled_value)
        else:
            theta[parameter] = int(round(sampled_value))

    return theta

def create_theta_uniform(parameter_vector):
    rng = np.random.default_rng()
    theta = {}

    for parameter, values in parameter_vector.items():
        if parameter == "SPLIT_TIME":
            theta[parameter] = values
            continue
        lower, upper = values
        if parameter == "p_AFR":
            theta[parameter] = rng.uniform(lower, upper)
        else:
            theta[parameter] = rng.integers(lower, upper)
    return theta

"""
========================================
Likelihood estimation for a theta vector
========================================
Samples one theta vector from the selected prior and estimates its likelihood.

"""

def simulate_theta_replicates(parameter_vectors, n_replicates, prior_type):
    success = 0
    parameter_vector = parameter_vectors[prior_type]


    if prior_type == "uniform":
        theta = create_theta_uniform(parameter_vector)
    else:
        theta = create_theta_normal(parameter_vector)

    for seed in range(1, n_replicates+1):
        ts = define_demography_AFR_EUR(theta, seed)
        cnt_EUR = count_origin("EUR", ts, theta["ADMIX_TIME"])/1000
        if 0.08 <= cnt_EUR <= 0.12:
            success += 1

    return (theta, success/n_replicates)



parameter_vectors = {
    "uniform": {
        "N_ANC": (7000, 20001),
        "N_AFR": (10000, 50001),
        "N_EUR": (5000, 30001),
        "N_ADMIX": (5000, 50001),
        "ADMIX_TIME": (8, 17),
        "SPLIT_TIME": 2000,
        "p_AFR": (0.7, 0.9),
    },  
    "normal": {
        "N_ANC": (5000, 10000, 15000),
        "N_AFR": (10_000, 20000, 30000),
        "N_EUR": (5000, 10000, 15000),
        "N_ADMIX": (5000, 20000, 35000),
        "ADMIX_TIME": (5, 10, 15),
        "p_AFR": (0.70, 0.80, 0.90),
        "SPLIT_TIME": 2000,
    },
}

prior_type="normal"
n_replicates = 100
n_vectors = 1
results = []
for vector_id in range(0, n_vectors):
    print(vector_id)
    tmp = simulate_theta_replicates(parameter_vectors, n_replicates, prior_type=prior_type)
    results.append({
        "vector_id": vector_id,
        **tmp[0],
        "p_EUR": 1 - tmp[0]['p_AFR'],
        "success_rate | likelihood": tmp[1],
        })
    
results_df = pd.DataFrame(results)
results_df.columns = results_df.columns.str.replace('parameters.', '', regex=False)
total_sum = results_df['success_rate | likelihood'].sum()
results_df['posterior'] = results_df['success_rate | likelihood'] / (total_sum)
output_file = f"output_{n_vectors}_{prior_type}.csv"

results_df.to_csv(
    output_file,
    index=False,
)









    