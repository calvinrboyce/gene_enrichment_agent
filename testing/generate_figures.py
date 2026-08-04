import pickle
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from scipy.stats import pearsonr, spearmanr


def percentile_distribution(results_dir, ax):
    # plot the distribution of percentiles in the test cases in a violin plot
    with open(results_dir + "/results_with_descriptions.pkl", "rb") as f:
        results = pickle.load(f)

    test_sem_sim = [float(test_case["sem_sim"]) for test_case in results if "sem_sim" in test_case]
    test_percentiles = [float(test_case["percentile"]) for test_case in results if "percentile" in test_case]
    control_sem_sim = [float(test_case["control_sem_sim"]) for test_case in results if "control_sem_sim" in test_case]
    control_percentiles = [float(test_case["control_percentile"]) for test_case in results if "control_percentile" in test_case]

    data_to_plot = [test_sem_sim, control_sem_sim, test_percentiles, control_percentiles]

    parts = ax.violinplot(data_to_plot, showmeans=True)
    
    test_color = 'tab:blue'
    control_color = 'tab:green'
    
    for i, pc in enumerate(parts['bodies']):
        if i % 2 == 0:
            pc.set_facecolor(test_color)
        else:
            pc.set_facecolor(control_color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)
        
    for partname in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
        vp = parts[partname]
        vp.set_edgecolor('black')
        vp.set_linewidth(1)

    import matplotlib.patches as mpatches
    test_patch = mpatches.Patch(color=test_color, alpha=0.7, label='Test')
    control_patch = mpatches.Patch(color=control_color, alpha=0.7, label='Control')
    ax.legend(handles=[test_patch, control_patch], loc='best')

    # # add a table of statistics to the plot
    # table_data = [
    #     ["Mean", "Median", "Std Dev"],
    #     [f"{np.mean(test_percentiles):.2f}", f"{np.median(test_percentiles):.2f}", f"{np.std(test_percentiles):.2f}"]
    # ]
    # table = ax.table(
    #     cellText=table_data,
    #     cellLoc='center',
    #     # bbox=[.55, .153, .4, .2]  # msigdb
    #     # bbox=[.55, .2155, .4, .2]  # go
    #     # bbox=[.55, .106, .4, .2]  # reactome
    #     bbox=[.55, .17, .4, .2]  # panglao
    # )
    # for key, cell in table.get_celld().items():
    #     cell.set_edgecolor('lightgrey')  # black borders for better visibility
    #     # cell.set_linewidth(1.0)
    #     cell.set_facecolor('white')  # white background for opacity
    #     # cell.set_alpha(0.9)  # slight transparency

    ax.set_ylabel('Scores', fontsize=11)
    ax.set_title('Distribution of Test Case Percentiles', fontsize=12)
    ax.set_xticks([1.5, 3.5])
    ax.set_xticklabels(["Semantic Similarity", "Percentile"], fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.grid(True, alpha=0.3, axis='y')


def gene_list_histogram(results_dir, ax):
    with open(results_dir + "/results.pkl", "rb") as f:
        results = pickle.load(f)

    gene_list_lengths = [len(result["genes"]) for result in results if "enrichment_results" in result]
    ax.hist(gene_list_lengths, log=True, bins=20, alpha=0.7, edgecolor='black')
    ax.set_xlabel("Number of genes")
    ax.set_ylabel("Number of terms (log scale)")
    ax.grid(True, alpha=0.3)


def cosine_similarity_distributions(results_dir, ax):
    with open(results_dir + "/results_with_descriptions.pkl", "rb") as f:
        results = pickle.load(f)

    with open(results_dir + "/name_to_index.pkl", "rb") as f:
        name_to_index = pickle.load(f)

    data_to_plot = np.array([result["cosine_similarities"] for result in results if "cosine_similarities" in result])
    
    scatter_handles = []
    # randomly sample 10 indices
    indices = np.random.choice(len(data_to_plot), 10, replace=False)
    ax.violinplot([data_to_plot[idx] for idx in indices], showmedians=True)
    for i, idx in enumerate(indices):
        result = results[idx]
        similarities = result["cosine_similarities"]
        print(result["name"])
        # print(result["most_similar_theme"]["theme"])
        # print(result["most_similar_theme"]["description"])
        # print(result["percentile"])
        # print()
        # Only add the label to the first scatter, so only one legend entry
        if i == 0:
            handle = ax.scatter(i+1, similarities[name_to_index[result["name"]]], color='red', zorder=5, 
                               label="Proposed Name")
            scatter_handles.append(handle)
        else:
            ax.scatter(i+1, similarities[name_to_index[result["name"]]], color='red', zorder=5)
    print('-'*100)
    ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    ax.set_xticklabels(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], fontsize=10)
    ax.set_xlabel("Test Case", fontsize=11)
    ax.set_ylabel("Semantic Similarity Score", fontsize=11)
    ax.set_title("Distribution of Similarity Scores", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=10)
    if scatter_handles:
        ax.legend(handles=scatter_handles, fontsize=10, frameon=True, fancybox=True, shadow=True)


def recovery_arc(results_dir, ax):
    with open(results_dir + "/results_with_descriptions.pkl", "rb") as f:
        results = pickle.load(f)
    
    percentiles = np.array([result["percentile"]*100 for result in results if "percentile" in result])
    terms_recovered = [np.mean(percentiles >= i) for i in range(101)]
    ax.plot(terms_recovered)

    # Find the value at x=95
    x_val = 95
    y_val = np.mean(percentiles >= x_val)
    print(results[0]["cosine_similarities"].shape)
    print(len(percentiles))
    print(np.sum(percentiles >= x_val))

    # Dotted line from x axis up to the curve at x=95
    ax.plot([x_val, x_val], [0, y_val], color='red', linestyle='--')
    ax.text(91, 0.01, f"{x_val}", color='red')

    # Dotted line from y axis over to the curve at y=y_val
    ax.plot([0, x_val], [y_val, y_val], color='red', linestyle='--')
    ax.text(0.5, y_val+0.015, f"{y_val:.2f}", color='red')

    ax.set_xlim(0, 101)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Similarity Percentile", fontsize=11)
    ax.set_ylabel("Percentage of Terms Recovered", fontsize=11)
    ax.set_title("Recovery Curve", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=10)


def percentile_vs_term_length(results_dir, ax):
    with open(results_dir + "/results_with_descriptions.pkl", "rb") as f:
        results = pickle.load(f)
    
    term_lengths = [len(result["genes"]) for result in results if "percentile" in result]
    percentiles = [result["percentile"] for result in results if "percentile" in result]
    ax.scatter(term_lengths, percentiles, alpha=0.6, s=30, color='#1f77b4', linewidth=0.5)
    ax.set_xlabel("Number of Genes", fontsize=11)
    ax.set_ylabel("Similarity Percentile", fontsize=11)
    ax.set_title("Term Length vs. Percentile", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=10)


def confidence_distribution(ax):
    with open("panglao/results_with_descriptions.pkl", "rb") as f:
        panglao_results = pickle.load(f)
    with open("reactome/results_with_descriptions.pkl", "rb") as f:
        reactome_results = pickle.load(f)
    with open("go/results_with_descriptions.pkl", "rb") as f:
        go_results = pickle.load(f)
    with open("msigdb/results_with_descriptions.pkl", "rb") as f:
        msigdb_results = pickle.load(f)
    with open("shortlist/results_with_descriptions.pkl", "rb") as f:
        shortlist_results = pickle.load(f)
    
    panglao_confidences = [theme["confidence"] for result in panglao_results if "enrichment_results" in result for theme in result["enrichment_results"]["themes"]]
    reactome_confidences = [theme["confidence"] for result in reactome_results if "enrichment_results" in result for theme in result["enrichment_results"]["themes"]]
    go_confidences = [theme["confidence"] for result in go_results if "enrichment_results" in result for theme in result["enrichment_results"]["themes"]]
    msigdb_confidences = [theme["confidence"] for result in msigdb_results if "enrichment_results" in result for theme in result["enrichment_results"]["themes"]]
    shortlist_confidences = [theme["confidence"] for result in shortlist_results if "enrichment_results" in result for theme in result["enrichment_results"]["themes"]]
    combined_confidences = panglao_confidences + reactome_confidences + go_confidences + msigdb_confidences + shortlist_confidences
    
    ax.violinplot(combined_confidences, orientation='horizontal')
    low_medium_split = .87
    medium_high_split = .93
    ax.axvline(x=low_medium_split, color='r', linestyle='--', linewidth=1.5)
    ax.axvline(x=medium_high_split, color='r', linestyle='--', linewidth=1.5)
    ax.text(low_medium_split - .045, .73, f"{low_medium_split:.2f}", color='r')
    ax.text(medium_high_split + .005, .73, f"{medium_high_split:.2f}", color='r')
    # ax.hist(combined_confidences, bins=20, alpha=0.7, edgecolor='black')
    # ax.hist(combined_confidences, bins=[0, 0.87, 0.93, 1], alpha=0.7, edgecolor='black')
    ax.set_xlabel("Confidence", fontsize=11)
    ax.set_ylabel("Number of terms", fontsize=11)
    ax.set_title("Confidence Score Distribution", fontsize=12)
    ax.set_yticks([])
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=10)

def nonsense_figure(ax):
    with open("reactome/results_with_descriptions.pkl", "rb") as f:
        reactome_results = pickle.load(f)
    
    with open("nonsense/results_with_descriptions.pkl", "rb") as f:
        nonsense_results = pickle.load(f)
    
    confidence_scores = defaultdict(lambda: [0, 0, 0])
    percentile_scores = defaultdict(lambda: [0, 0, 0])
    # add nonsense results
    for result in nonsense_results:
        if result['name'].startswith('Action of antimicrobials'):
            continue
        if "enrichment_results" in result:
            if result["name"].endswith(" (50/50 mix)"):
                name = result["name"].rsplit(" (", 1)[0]
                # insert confidence score in the front of the list
                confidence_scores[name][1] = result["most_similar_theme"]["confidence"]
                percentile_scores[name][1] = result["percentile"]
            elif result["name"].endswith(" (random)"):
                name = result["name"].rsplit(" (", 1)[0]
                confidence_scores[name][2] = result["most_similar_theme"]["confidence"]
                percentile_scores[name][2] = result["percentile"]
    
    # add reactome results
    for result in reactome_results:
        try:
            if result["name"] in confidence_scores:
                confidence_scores[result["name"]][0] = result["most_similar_theme"]["confidence"]
                percentile_scores[result["name"]][0] = result["percentile"]
        except Exception as e:
            print(result["name"])
            print(result.keys())
            raise e
    
    # bin the confidence scores
    def confidence_bin(score):
        if score <= 0.87:
            return 0
        elif score <= 0.93:
            return 1
        else:
            return 2
    
    # confidence bins
    pure_confidence_bins = np.bincount([confidence_bin(value[0]) for value in confidence_scores.values()])
    pure_confidence_bins = pure_confidence_bins / np.sum(pure_confidence_bins)
    print(np.sum(pure_confidence_bins[1:]))
    mixed_confidence_bins = np.bincount([confidence_bin(value[1]) for value in confidence_scores.values()])
    mixed_confidence_bins = mixed_confidence_bins / np.sum(mixed_confidence_bins)
    print(np.sum(mixed_confidence_bins[1:]))
    random_confidence_bins = np.bincount([confidence_bin(value[2]) for value in confidence_scores.values()])
    random_confidence_bins = random_confidence_bins / np.sum(random_confidence_bins)
    print(np.sum(random_confidence_bins[1:]))
    if len(random_confidence_bins) == 2:
        random_confidence_bins = np.append(random_confidence_bins, 0)

    # recovery scores
    pure_percentiles = np.array([value[0] for value in percentile_scores.values()])
    pure_recovery = np.mean(pure_percentiles >= 0.95)
    mixed_percentiles = np.array([value[1] for value in percentile_scores.values()])
    mixed_recovery = np.mean(mixed_percentiles >= 0.95)
    random_percentiles = np.array([value[2] for value in percentile_scores.values()])
    random_recovery = np.mean(random_percentiles >= 0.95)

    # Stack the bars
    labels = [f'Pure Test Cases', f'Mixed Test Cases', f'Random Test Cases']
    categories = ['High', 'Medium', 'Low']
    
    # Professional color scheme: light gray, light blue, dark blue
    colors = ['#1F77B4', '#9BCBEE', '#E0E0E0']  # Light gray, light blue, dark blue
    
    bottoms = np.zeros(len(categories))

    for i in range(3):
        ax.bar(labels, [pure_confidence_bins[::-1][i], mixed_confidence_bins[::-1][i], random_confidence_bins[::-1][i]], 
               bottom=bottoms, label=categories[i], color=colors[i])
        bottoms += [pure_confidence_bins[::-1][i], mixed_confidence_bins[::-1][i], random_confidence_bins[::-1][i]]

    # Add red line for recovery
    ax.plot([0, 1, 2], [pure_recovery, mixed_recovery, random_recovery], color='red', marker='o', linewidth=2, label='Recovery')

    ax.set_ylabel('Proportion', fontsize=11)
    ax.set_title('Confidence Evaluation with Noise', fontsize=12)
    ax.legend(title='Confidence Level', title_fontsize=10, fontsize=9, frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Improve tick labels
    ax.tick_params(axis='both', which='major', labelsize=10)

def calibration_analysis(ax, confidences, accuracies, num_bins=10):
    """
    Computes Adaptive ECE, correlation metrics, and plots a quantile-based reliability diagram.
    
    Parameters:
    confidences (array-like): Model confidence scores [0, 1]
    accuracies (array-like): Ground truth accuracy or similarity scores [0, 1]
    num_bins (int): Number of equal-mass bins (quantiles)
    """
    # Convert to numpy arrays for safety
    confidences = np.array(confidences)
    accuracies = np.array(accuracies)
    
    # 1. Sort both arrays based on confidence scores
    sorted_indices = np.argsort(confidences)
    sorted_confs = confidences[sorted_indices]
    sorted_accs = accuracies[sorted_indices]
    
    # 2. Split into equal-mass bins (quantiles)
    # np.array_split handles arrays that don't divide perfectly by num_bins
    conf_bins = np.array_split(sorted_confs, num_bins)
    acc_bins = np.array_split(sorted_accs, num_bins)
    
    # 3. Calculate mean confidence and mean accuracy per bin
    mean_confs = np.array([np.mean(b) for b in conf_bins])
    mean_accs = np.array([np.mean(b) for b in acc_bins])
    bin_sizes = np.array([len(b) for b in conf_bins])
    
    # 4. Compute Adaptive Expected Calibration Error (ACE)
    n = len(confidences)
    ace = np.sum((bin_sizes / n) * np.abs(mean_accs - mean_confs))
    
    # 5. Compute Correlation Metrics
    # Pearson: Linear relationship (Do they scale together proportionally?)
    pearson_corr, _ = pearsonr(confidences, accuracies)
    # Spearman: Monotonic relationship (Does higher confidence reliably mean higher accuracy, regardless of scale?)
    spearman_corr, _ = spearmanr(confidences, accuracies)
    
    # Perfect calibration reference line
    ax.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Perfect Calibration')
    
    # Model calibration curve
    ax.plot(mean_confs, mean_accs, marker='o', color='tab:blue', 
            linewidth=2, markersize=8, label='Model Calibration')
    
    # Formatting the plot
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_xlabel('Mean Confidence', fontsize=12)
    ax.set_ylabel('Mean Similarity Percentile', fontsize=12)
    ax.set_title('Reliability Diagram', fontsize=12)
    
    # Add metrics text box to the plot
    metrics_text = (
        f"Adaptive ECE: {ace:.4f}\n"
        f"Pearson r: {pearson_corr:.3f}\n"
        f"Spearman ρ: {spearman_corr:.3f}"
    )
    ax.text(0.05, 0.95, metrics_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(facecolor='white', alpha=0.9, edgecolor='lightgray'))
    
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)

def percentile_vs_confidence(ax, confidences, percentiles):
    ax.scatter(confidences, percentiles, alpha=0.6, s=30, color='#1f77b4', linewidth=0.5)
    ax.set_xlabel("Confidence", fontsize=11)
    ax.set_ylabel("Similarity Percentile", fontsize=11)
    ax.set_title("Confidence vs. Percentile", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=10)


def characterize_results(results_dir, best=True):
    with open(results_dir + "/results.pkl", "rb") as f:
        results = pickle.load(f)
    
    results = sorted([result for result in results if "percentile" in result], key=lambda x: x["percentile"], reverse=best)

    for result in results[:10]:
        print("Name: ", result["name"])
        print("Number of genes: ", len(result["genes"]))
        print()
        # print("Full summary: ", result["enrichment_results"]["summary"])
        # print()
        print("Proposed themes:\n" + '\n'.join([theme["theme"] for theme in result["enrichment_results"]["themes"]]))
        print()
        print("Most similar theme: ", result["most_similar_theme"]["theme"])
        print("Description: ", result["most_similar_theme"]["description"])
        print("Percentile: ", result["percentile"])
        print()
        print('-'*100)


def generate_plots(dataset):
    """
    Generate a comprehensive figure for a given dataset suitable for scientific papers.
    
    Args:
        dataset (str): Name of the dataset directory (e.g., "msigdb", "reactome", "panglao")
    
    Returns:
        matplotlib.figure.Figure: The generated figure
    """
    # Set up the figure with professional formatting
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    dir_to_name = {
        'go': 'GO Biological Processes',
        'msigdb': 'MSigDB Hallmark Gene Sets',
        'reactome': 'Reactome Pathways',
        'panglao': 'PanglaoDB Cell Markers',
        'shortlist': 'Shortlist from Reactome Pathways'
    }
    print()
    print("Generating plots for dataset: ", dataset)
    fig.suptitle(dir_to_name[dataset], 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Generate all plots
    percentile_distribution(dataset, axs[0, 0])
    percentile_vs_term_length(dataset, axs[0, 1])
    cosine_similarity_distributions(dataset, axs[1, 0])
    recovery_arc(dataset, axs[1, 1])
    
    # Add subplot labels (A, B, C, D)
    subplot_labels = ['a', 'b', 'c', 'd']
    for i, ax in enumerate(axs.flat):
        ax.text(-0.1, 1.05, subplot_labels[i], transform=ax.transAxes, 
                fontsize=14, fontweight='bold', va='top')
    
    # Professional formatting
    plt.tight_layout()
    plt.subplots_adjust(top=0.92, hspace=0.3, wspace=0.3)
    plt.savefig(f"{dataset}/plots_with_descriptions.png", dpi=300, bbox_inches='tight')
    
    return fig



# Generate plots for different datasets
# fig = generate_plots("go")
# fig = generate_plots("msigdb")
# fig = generate_plots("reactome")
# fig = generate_plots("panglao")
fig = generate_plots("shortlist")


# # plot confidence distribution and nonsense 

# # pull data
# with open("panglao/results_with_descriptions.pkl", "rb") as f:
#     panglao_results = pickle.load(f)
# with open("reactome/results_with_descriptions.pkl", "rb") as f:
#     reactome_results = pickle.load(f)
# with open("go/results_with_descriptions.pkl", "rb") as f:
#     go_results = pickle.load(f)
# with open("msigdb/results_with_descriptions.pkl", "rb") as f:
#     msigdb_results = pickle.load(f)
# with open("shortlist/results_with_descriptions.pkl", "rb") as f:
#     shortlist_results = pickle.load(f)
# with open("nonsense/results_with_descriptions.pkl", "rb") as f:
#     nonsense_results = pickle.load(f)

# # confidences
# panglao_confidences = [result["most_similar_theme"]["confidence"] for result in panglao_results if "most_similar_theme" in result]
# reactome_confidences = [result["most_similar_theme"]["confidence"] for result in reactome_results if "most_similar_theme" in result]
# go_confidences = [result["most_similar_theme"]["confidence"] for result in go_results if "most_similar_theme" in result]
# msigdb_confidences = [result["most_similar_theme"]["confidence"] for result in msigdb_results if "most_similar_theme" in result]
# shortlist_confidences = [result["most_similar_theme"]["confidence"] for result in shortlist_results if "most_similar_theme" in result]
# nonsense_confidences = [result["most_similar_theme"]["confidence"] for result in nonsense_results if "most_similar_theme" in result]
# combined_confidences = panglao_confidences + reactome_confidences + go_confidences + msigdb_confidences + shortlist_confidences + nonsense_confidences

# # percentiles
# panglao_percentiles = [result["percentile"] for result in panglao_results if "percentile" in result]
# reactome_percentiles = [result["percentile"] for result in reactome_results if "percentile" in result]
# go_percentiles = [result["percentile"] for result in go_results if "percentile" in result]
# msigdb_percentiles = [result["percentile"] for result in msigdb_results if "percentile" in result]
# shortlist_percentiles = [result["percentile"] for result in shortlist_results if "percentile" in result]
# nonsense_percentiles = [result["percentile"] for result in nonsense_results if "percentile" in result]
# combined_percentiles = panglao_percentiles + reactome_percentiles + go_percentiles + msigdb_percentiles + shortlist_percentiles + nonsense_percentiles



# fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# # Add subplot labels (a, b)
# subplot_labels = ['a', 'b'  , 'c', 'd']
# for i, ax in enumerate(axs.flat):
#     ax.text(-0.1, 1.05, subplot_labels[i], transform=ax.transAxes, 
#             fontsize=14, fontweight='bold', va='top')

# confidence_distribution(axs[0, 0])
# nonsense_figure(axs[0, 1])
# calibration_analysis(axs[1,0], combined_confidences, combined_percentiles)
# percentile_vs_confidence(axs[1,1], combined_confidences, combined_percentiles)
# plt.tight_layout()
# plt.subplots_adjust(top=0.92, hspace=0.3)
# fig.suptitle("Confidence Calibration", fontsize=16, fontweight='bold', y=0.98)
# plt.savefig("testing/confidence_distribution.png", dpi=300, bbox_inches='tight')
