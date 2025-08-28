import pickle
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict


def percentile_distribution(results_dir, ax):
    # plot the distribution of percentiles in the test cases in a violin plot
    with open(results_dir + "/results.pkl", "rb") as f:
        results = pickle.load(f)

    test_percentiles = [float(test_case["percentile"])*100 for test_case in results if "percentile" in test_case]

    data_to_plot = [test_percentiles]

    ax.violinplot(data_to_plot)

    # add a table of statistics to the plot
    table_data = [
        ["Mean", "Median", "Std Dev"],
        [f"{np.mean(test_percentiles):.2f}", f"{np.median(test_percentiles):.2f}", f"{np.std(test_percentiles):.2f}"]
    ]
    table = ax.table(
        cellText=table_data,
        cellLoc='center',
        # bbox=[.55, .153, .4, .2]  # msigdb
        # bbox=[.55, .2155, .4, .2]  # go
        # bbox=[.55, .106, .4, .2]  # reactome
        bbox=[.55, .17, .4, .2]  # panglao
    )
    for key, cell in table.get_celld().items():
        cell.set_edgecolor('lightgrey')  # black borders for better visibility
        # cell.set_linewidth(1.0)
        cell.set_facecolor('white')  # white background for opacity
        # cell.set_alpha(0.9)  # slight transparency

    ax.set_ylabel('Percentile Scores', fontsize=11)
    ax.set_title('Distribution of Test Case Percentiles', fontsize=12)
    ax.set_xticks([1])
    ax.set_xticklabels(["PanglaoDB"], fontsize=10)
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
    with open(results_dir + "/results.pkl", "rb") as f:
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
    with open(results_dir + "/results.pkl", "rb") as f:
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
    with open(results_dir + "/results.pkl", "rb") as f:
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
    with open("panglao/results.pkl", "rb") as f:
        panglao_results = pickle.load(f)
    with open("reactome/results.pkl", "rb") as f:
        reactome_results = pickle.load(f)
    with open("go/results.pkl", "rb") as f:
        go_results = pickle.load(f)
    with open("msigdb/results.pkl", "rb") as f:
        msigdb_results = pickle.load(f)
    
    panglao_confidences = [[theme["confidence"] for theme in result["enrichment_results"]["themes"]] for result in panglao_results if "enrichment_results" in result]
    reactome_confidences = [[theme["confidence"] for theme in result["enrichment_results"]["themes"]] for result in reactome_results if "enrichment_results" in result]
    go_confidences = [[theme["confidence"] for theme in result["enrichment_results"]["themes"]] for result in go_results if "enrichment_results" in result]
    msigdb_confidences = [[theme["confidence"] for theme in result["enrichment_results"]["themes"]] for result in msigdb_results if "enrichment_results" in result]
    flattened_panglao_confidences = [item for sublist in panglao_confidences for item in sublist]
    flattened_reactome_confidences = [item for sublist in reactome_confidences for item in sublist]
    flattened_go_confidences = [item for sublist in go_confidences for item in sublist]
    flattened_msigdb_confidences = [item for sublist in msigdb_confidences for item in sublist]
    combined_confidences = flattened_panglao_confidences + flattened_reactome_confidences + flattened_go_confidences + flattened_msigdb_confidences
    
    ax.violinplot(combined_confidences, orientation='horizontal')
    low_medium_split = .87
    medium_high_split = .93
    ax.axvline(x=low_medium_split, color='r', linestyle='--', linewidth=1.5)
    ax.axvline(x=medium_high_split, color='r', linestyle='--', linewidth=1.5)
    ax.text(low_medium_split - .035, .73, f"{low_medium_split:.2f}", color='r')
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
    with open("reactome/results.pkl", "rb") as f:
        reactome_results = pickle.load(f)
    
    with open("nonsense/results.pkl", "rb") as f:
        nonsense_results = pickle.load(f)
    
    confidence_scores = defaultdict(lambda: [0, 0, 0])
    percentile_scores = defaultdict(lambda: [0, 0, 0])
    # add nonsense results
    for result in nonsense_results:
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
        if result["name"] in confidence_scores:
            confidence_scores[result["name"]][0] = result["most_similar_theme"]["confidence"]
            percentile_scores[result["name"]][0] = result["percentile"]
    
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

    # recovery scores
    pure_percentiles = np.array([value[0] for value in percentile_scores.values()])
    pure_recovery = np.mean(pure_percentiles >= 0.95)
    mixed_percentiles = np.array([value[1] for value in percentile_scores.values()])
    mixed_recovery = np.mean(mixed_percentiles >= 0.95)
    random_percentiles = np.array([value[2] for value in percentile_scores.values()])
    random_recovery = np.mean(random_percentiles >= 0.95)

    # Stack the bars
    labels = [f'Pure Test Cases\nRecovery: {pure_recovery:.2f}', f'Mixed Test Cases\nRecovery: {mixed_recovery:.2f}', f'Random Test Cases\nRecovery: {random_recovery:.2f}']
    categories = ['High', 'Medium', 'Low']
    
    # Professional color scheme: light gray, light blue, dark blue
    colors = ['#1F77B4', '#9BCBEE', '#E0E0E0']  # Light gray, light blue, dark blue
    
    bottoms = np.zeros(len(categories))

    for i in range(3):
        ax.bar(labels, [pure_confidence_bins[::-1][i], mixed_confidence_bins[::-1][i], random_confidence_bins[::-1][i]], 
               bottom=bottoms, label=categories[i], color=colors[i])
        bottoms += [pure_confidence_bins[::-1][i], mixed_confidence_bins[::-1][i], random_confidence_bins[::-1][i]]

    ax.set_ylabel('Proportion', fontsize=11)
    ax.set_title('Confidence Evaluation with Noise', fontsize=12)
    ax.legend(title='Confidence Level', title_fontsize=10, fontsize=9, frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Improve tick labels
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
    # fig.suptitle(f'Gene Enrichment Analysis Results: {dataset.upper()}', 
    #              fontsize=16, fontweight='bold', y=0.98)
    
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
    plt.savefig(f"testing/{dataset}_plots.png", dpi=300, bbox_inches='tight')
    
    return fig



# Generate plots for different datasets
fig = generate_plots("panglao")

# # plot confidence distribution and nonsense figure
# fig, axs = plt.subplots(2, 1, figsize=(10, 12))

# # Add subplot labels (a, b)
# subplot_labels = ['a', 'b']
# for i, ax in enumerate(axs):
#     ax.text(-0.1, 1.05, subplot_labels[i], transform=ax.transAxes, 
#             fontsize=14, fontweight='bold', va='top')

# confidence_distribution(axs[0])
# nonsense_figure(axs[1])
# plt.tight_layout()
# plt.subplots_adjust(top=0.92, hspace=0.3)
# plt.savefig("testing/confidence_distribution.png", dpi=300, bbox_inches='tight')
