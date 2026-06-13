import msprime
import pandas as pd

def define_demography(N, ADMIX_TIME, SPLIT_TIME, p_A, seed):
    
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N)
    demography.add_population(name="B", initial_size=N)
    demography.add_population(name="ADMIX", initial_size=N)
    demography.add_population(name="ANC", initial_size=N)
    
    demography.add_population_split(time=SPLIT_TIME, derived=["A", "B"], ancestral="ANC")

    demography.add_admixture(time=ADMIX_TIME, derived="ADMIX", ancestral=["A", "B"], proportions=[p_A, 1 - p_A]) 
    
    demography.sort_events()
    
    ts = msprime.sim_ancestry(
        samples={"ADMIX": N},
        demography=demography,
        ploidy=1,
        sequence_length=16569,
        recombination_rate=0,
        random_seed=seed,
        record_migrations=True,
    )
    return ts

def count_origin_A(ts, ADMIX_TIME):
    cnt_A = 0
    tree = ts.first()
    for sample in ts.samples():
        origin = None
        for m in ts.migrations():
            if m.time == ADMIX_TIME:
                if tree.is_descendant(sample, m.node):
                    origin = ts.population(m.dest).metadata["name"]
                    if origin == "A":
                        cnt_A += 1
                    break
    

    return cnt_A

def prob_A_above_threshold(times, N, reps, threshold, SPLIT_TIME, p_A):
    results = {}
    for t in times:
        success = 0
        for seed in range(1, reps+1):
            ts = define_demography(N, t, SPLIT_TIME, p_A, seed)
            f_A = count_origin_A(ts, t)/N
            if f_A >= threshold:
                success += 1
        results[t] = success/reps
        
     
    return results


            
            
            
