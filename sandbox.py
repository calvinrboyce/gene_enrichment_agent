import pickle
import numpy as np
from scipy import stats

with open("barcoding/results_with_descriptions.pkl", "rb") as f:
    results = pickle.load(f)

term_hallucination_rates = [t['enrichment_results']['hallucination_metrics']['term_hallucination_rate'] for t in results if 'enrichment_results' in t]
pvalue_hallucination_rates = [t['enrichment_results']['hallucination_metrics']['pvalue_hallucination_rate'] for t in results if 'enrichment_results' in t]
runtimes = [t['enrichment_results']['runtime'] for t in results if 'enrichment_results' in t]
total_tokens = [t['enrichment_results']['token_usage']['total_tokens'] for t in results if 'enrichment_results' in t]

def print_summary_statistics(data):
    # 1. Calculate the Mean
    mean_rate = np.mean(data)

    # 2. Calculate the Median (50th percentile)
    median_rate = np.median(data)

    # 3. Calculate the 95th Percentile (tail risk)
    p95_rate = np.percentile(data, 95)

    # 4. Calculate the Bootstrapped 95% Confidence Interval for the mean
    # scipy.stats.bootstrap expects the data in a sequence, so we pass it as a tuple: (rates,)
    bootstrap_ci = stats.bootstrap(
        (data,), 
        statistic=np.mean, 
        confidence_level=0.95,
        method='BCa'  # Bias-corrected and accelerated method (great for skewed data)
    )

    ci_lower = bootstrap_ci.confidence_interval.low
    ci_upper = bootstrap_ci.confidence_interval.high

    # Print the formatted results
    print(f"Mean:            {mean_rate}")
    print(f"95% CI (Mean):   [{ci_lower}, {ci_upper}]")
    print(f"Median:          {median_rate}")
    print(f"95th Percentile: {p95_rate}")

# print('Term Statistics')
# print_summary_statistics(term_hallucination_rates)
# print()
# print('P-value Statistics')
# print_summary_statistics(pvalue_hallucination_rates)
# print()
# print('Runtime Statistics')
# print_summary_statistics(runtimes)
# print()
# print('Token Usage Statistics')
# print_summary_statistics(total_tokens)
# print()

with open('go/results_with_descriptions.pkl', 'rb') as f:
    go_results = pickle.load(f)

with open('msigdb/results_with_descriptions.pkl', 'rb') as f:
    msigdb_results = pickle.load(f)

with open('reactome/results_with_descriptions.pkl', 'rb') as f:
    reactome_results = pickle.load(f)

with open('panglao/results_with_descriptions.pkl', 'rb') as f:
    panglao_results = pickle.load(f)

with open('shortlist/results_with_descriptions.pkl', 'rb') as f:
    shortlist_results = pickle.load(f)

control_runtimes = []
control_runtimes.extend([t['enrichment_results']['runtime'] for t in go_results if 'enrichment_results' in t])
control_runtimes.extend([t['enrichment_results']['runtime'] for t in msigdb_results if 'enrichment_results' in t])
control_runtimes.extend([t['enrichment_results']['runtime'] for t in reactome_results if 'enrichment_results' in t])
control_runtimes.extend([t['enrichment_results']['runtime'] for t in panglao_results if 'enrichment_results' in t])
control_runtimes.extend([t['enrichment_results']['runtime'] for t in shortlist_results if 'enrichment_results' in t])

# print('Control Runtime Statistics')
# print_summary_statistics(control_runtimes)
# print()

control_input_token_usage = []
control_input_token_usage.extend([t['enrichment_results']['token_usage']['input_tokens'] for t in go_results if 'enrichment_results' in t])
control_input_token_usage.extend([t['enrichment_results']['token_usage']['input_tokens'] for t in msigdb_results if 'enrichment_results' in t])
control_input_token_usage.extend([t['enrichment_results']['token_usage']['input_tokens'] for t in reactome_results if 'enrichment_results' in t])
control_input_token_usage.extend([t['enrichment_results']['token_usage']['input_tokens'] for t in panglao_results if 'enrichment_results' in t])
control_input_token_usage.extend([t['enrichment_results']['token_usage']['input_tokens'] for t in shortlist_results if 'enrichment_results' in t])

# print('Control Input Token Usage Statistics')
# print_summary_statistics(control_input_token_usage)
# print()

control_output_token_usage = []
control_output_token_usage.extend([t['enrichment_results']['token_usage']['output_tokens'] for t in go_results if 'enrichment_results' in t])
control_output_token_usage.extend([t['enrichment_results']['token_usage']['output_tokens'] for t in msigdb_results if 'enrichment_results' in t])
control_output_token_usage.extend([t['enrichment_results']['token_usage']['output_tokens'] for t in reactome_results if 'enrichment_results' in t])
control_output_token_usage.extend([t['enrichment_results']['token_usage']['output_tokens'] for t in panglao_results if 'enrichment_results' in t])
control_output_token_usage.extend([t['enrichment_results']['token_usage']['output_tokens'] for t in shortlist_results if 'enrichment_results' in t])

# print('Control Output Token Usage Statistics')
# print_summary_statistics(control_output_token_usage)
# print()

control_total_token_usage = []
control_total_token_usage.extend([t['enrichment_results']['token_usage']['total_tokens'] for t in go_results if 'enrichment_results' in t])
control_total_token_usage.extend([t['enrichment_results']['token_usage']['total_tokens'] for t in msigdb_results if 'enrichment_results' in t])
control_total_token_usage.extend([t['enrichment_results']['token_usage']['total_tokens'] for t in reactome_results if 'enrichment_results' in t])
control_total_token_usage.extend([t['enrichment_results']['token_usage']['total_tokens'] for t in panglao_results if 'enrichment_results' in t])
control_total_token_usage.extend([t['enrichment_results']['token_usage']['total_tokens'] for t in shortlist_results if 'enrichment_results' in t])

# print('Control Total Token Usage Statistics')
# print(len(control_total_token_usage))
# print_summary_statistics(control_total_token_usage)
# print()

def print_results(results):
    sem_sims = [t['sem_sim'] for t in results if 'sem_sim' in t]
    percentiles = [t['percentile'] for t in results if 'percentile' in t]
    control_sem_sims = [t['control_sem_sim'] for t in results if 'control_sem_sim' in t]
    control_percentiles = [t['control_percentile'] for t in results if 'control_percentile' in t]
    ontological_distances = [t['ontological_distance'] for t in results if 'ontological_distance' in t]

    print('Sem Sims')
    print_summary_statistics(sem_sims)
    print()
    print('Percentiles')
    print_summary_statistics(percentiles)
    print()
    print('Control Sem Sims')
    print_summary_statistics(control_sem_sims)
    print()
    print('Control Percentiles')
    print_summary_statistics(control_percentiles)
    print()
    print('Ontological Distances')
    if ontological_distances:
        print_summary_statistics(ontological_distances)
    else:
        print('No ontological distances available')
    print()
    recovered_terms = (np.array(percentiles)>=.95).sum()
    total_terms_tested = len(percentiles)

    # Calculate the exact binomial test (defaults to 95% confidence level)
    result = stats.binomtest(k=recovered_terms, n=total_terms_tested)

    # Extract the Clopper-Pearson confidence interval
    ci_lower, ci_upper = result.proportion_ci(confidence_level=0.95)

    print(f"Recovery Rate: {recovered_terms / total_terms_tested:.1%}")
    print(f"95% CI:        [{ci_lower:.1%}, {ci_upper:.1%}]")

# print('Shortlist')
# print_results(shortlist_results)


global_percentiles = []
# global_percentiles.extend([t['percentile'] for t in go_results if 'percentile' in t])
# global_percentiles.extend([t['percentile'] for t in msigdb_results if 'percentile' in t])
# global_percentiles.extend([t['percentile'] for t in reactome_results if 'percentile' in t])
# global_percentiles.extend([t['percentile'] for t in panglao_results if 'percentile' in t])
global_percentiles.extend([t['percentile'] for t in shortlist_results if 'percentile' in t])

print('Global Percentiles')
print(len(global_percentiles))
print_summary_statistics(global_percentiles)
print()
recovered_terms = (np.array(global_percentiles)>=.95).sum()
total_terms_tested = len(global_percentiles)

# Calculate the exact binomial test (defaults to 95% confidence level)
result = stats.binomtest(k=recovered_terms, n=total_terms_tested)

# Extract the Clopper-Pearson confidence interval
ci_lower, ci_upper = result.proportion_ci(confidence_level=0.95)

print(f"Recovery Rate: {recovered_terms / total_terms_tested:.1%}")
print(f"95% CI:        [{ci_lower:.1%}, {ci_upper:.1%}]")
