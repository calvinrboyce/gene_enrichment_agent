from testing.go import get_go_terms
from testing.embedding import embed_terms
from gene_enrichment_agent import GeneEnrichmentAgent
import os
import dotenv
from tqdm import tqdm
import numpy as np
import pickle
import json
import random
import time
import networkx as nx
import gc

def validation_scheme(test_cases, database_terms, holdout, validation_name, go_graph=None):
    # we changed this function for each of the experiments to accomodate the different databases
    error_log = []
    # with open(validation_name + "/error_log.txt", "r") as f:
    #     error_log = f.read().strip().split('\n')

    # # make directory for results
    # os.makedirs(validation_name, exist_ok=True)

    # with open(validation_name + "/test_cases.json", "w") as f:
    #     json.dump(test_cases, f)

    # embed database terms
    print("Embedding database terms")
    db_name_to_index = {term["name"]: i for i, term in enumerate(database_terms)}
    with open(validation_name + "/name_to_index.pkl", "wb") as f:
        pickle.dump(db_name_to_index, f)
    db_names = [term["name"] for term in database_terms]
    with open(validation_name + "/db_names.pkl", "wb") as f:
        pickle.dump(db_names, f)
    db_embeddings = embed_terms(db_names)

    with open(validation_name + "/db_embeddings.pkl", "wb") as f:
        pickle.dump(db_embeddings, f)
    # with open(validation_name + "/db_embeddings.pkl", "rb") as f:
    #     db_embeddings = pickle.load(f)
    # with open(validation_name + "/name_to_index.pkl", "rb") as f:
    #     db_name_to_index = pickle.load(f)

    # run enrichment_agent on test_cases
    print("Running enrichment agent on test cases")
    dotenv.load_dotenv()
    open_ai_api_key = os.getenv("OPENAI_API_KEY")
    entrez_api_key = os.getenv("ENTREZ_API_KEY")

    gea = GeneEnrichmentAgent(open_ai_api_key, entrez_api_key)
    email = "cboyce3@mgh.harvard.edu"
    i = -1
    for test_case in tqdm(test_cases):
        # # error filling
        # if test_case['name'] not in error_log:
        #     continue
        # #
        if "enrichment_results" not in test_case:
            continue

        i += 1
        try:
            context = "None"
            use_barcodes = True
            if validation_name == 'panglao':
                context = "Include a theme for likely cell type"
            if validation_name == 'barcoding':
                use_barcodes = False
            if validation_name == 'go':
                holdout = test_case['id']
            # test_case["enrichment_results"] = gea.run_analysis(test_case["genes"], email, context=context, ranked=False, save_results=False, holdout=holdout, use_barcodes=use_barcodes)
            test_case["theme_embeddings"] = embed_terms(["Name: " + theme["theme"] + "\nDescription: " + theme["description"] for theme in test_case["enrichment_results"]["themes"]])
            test_case["control_embeddings"] = embed_terms([result["name"] for result in test_case["enrichment_results"]["top_statistical_results"]])
        except Exception as e:
            print(f'Enrichment failed for {test_case["name"]}:', e)
            error_log.append(test_case["name"])
            continue
        # save the test cases
        if not i % 10:
            try:
                with open(validation_name + "/results_with_descriptions.pkl", "wb") as f:
                    pickle.dump(test_cases, f)
            except Exception as e:
                print("Pickling failed at index " + str(i))
                with open(validation_name + "/error_log.txt", "w") as f:
                    for error in error_log:
                        f.write(error + "\n")
                return
    try:
        with open(validation_name + "/results_with_descriptions.pkl", "wb") as f:
            pickle.dump(test_cases, f)
    except Exception as e:
        print("Pickling failed at index " + str(i))
        with open(validation_name + "/error_log.txt", "w") as f:
            for error in error_log:
                f.write(error + "\n")
        return
    
    # save the error log
    print("Saving error log")
    with open(validation_name + "/error_log.txt", "w") as f:
        for error in error_log:
            f.write(error + "\n")
    
    # compute cosine similarities
    print("Computing cosine similarities and percentiles")
    for i, test_case in enumerate(test_cases):
        if "theme_embeddings" not in test_case:
            continue

        ys = db_embeddings
        if test_case["name"].endswith(" (50/50 mix)") or test_case["name"].endswith(" (random)"):
            y_index = db_name_to_index[test_case["name"].rsplit(" (", 1)[0]]
        else:
            y_index = db_name_to_index[test_case["name"]]
        y = ys[y_index]
        yhats = test_case["theme_embeddings"]

        # find the theme most similar to the correct db term
        cosine_similarities = np.dot(yhats, y) / (np.linalg.norm(y) * np.linalg.norm(yhats, axis=1))
        yhat_index = np.argmax(cosine_similarities)
        yhat = yhats[yhat_index]
        test_case["most_similar_theme"] = test_case["enrichment_results"]["themes"][yhat_index]

        # find the percentile of the actual similarity in the list of cosine similarities
        test_case["cosine_similarities"] = np.dot(ys, yhat) / (np.linalg.norm(ys, axis=1) * np.linalg.norm(yhat))
        actual_similarity = test_case["cosine_similarities"][y_index]
        percentile = (np.sum(test_case["cosine_similarities"] <= actual_similarity) - 1) / (len(test_case["cosine_similarities"]) - 1)
        test_case["sem_sim"] = actual_similarity
        test_case["percentile"] = percentile

        # closest GO term
        if validation_name == "go":
            closest_go_idx = np.argmax(test_case["cosine_similarities"])
            closest_go_term = database_terms[closest_go_idx]["id"]
            test_case["closest_go_term"] = closest_go_term
            try:
                distance = nx.shortest_path_length(go_graph, source=closest_go_term, target=test_case["id"])
                test_case["ontological_distance"] = distance
            except Exception as e:
                print(e)

        ## CONTROL SIMILARITY
        # find the raw result most similar to the correct term
        yhats = test_case["control_embeddings"]
        cosine_similarities = np.dot(yhats, y) / (np.linalg.norm(y) * np.linalg.norm(yhats, axis=1))
        yhat_index = np.argmax(cosine_similarities)
        yhat = yhats[yhat_index]
        test_case["most_similar_control"] = test_case["enrichment_results"]["top_statistical_results"][yhat_index]

        # find the percentile of the control similarity
        test_case["control_cosine_similarities"] = np.dot(ys, yhat) / (np.linalg.norm(ys, axis=1) * np.linalg.norm(yhat))
        actual_similarity = test_case["control_cosine_similarities"][y_index]
        percentile = (np.sum(test_case["control_cosine_similarities"] <= actual_similarity) - 1) / (len(test_case["control_cosine_similarities"]) - 1)
        test_case["control_sem_sim"] = actual_similarity
        test_case["control_percentile"] = percentile

        # save as we go
        if not i % 10:
            try:
                with open(validation_name + "/results_with_descriptions.pkl", "wb") as f:
                    pickle.dump(test_cases, f)
            except Exception as e:
                print("Cosine similarity pickling failed at index " + str(i))
                return

    # save the test cases
    try:
        print("Saving test cases")
        with open(validation_name + "/results_with_descriptions.pkl", "wb") as f:
            pickle.dump(test_cases, f)
    except Exception as e:
        print("Final pickling failed")
        return
    
    return test_cases

# # panglaoDB
# with open('testing/panglao_cells.gmt', 'r') as f:
#     panglao_cells = [line.strip().split('\t') for line in f.readlines()]
#     panglao_cells = [{
#         "name": term[0],
#         "genes": term[1:]
#     } for term in panglao_cells]
# with open("panglao/results.pkl", 'rb') as f:
#     panglao_test_cases = pickle.load(f)
# results = validation_scheme(panglao_test_cases, panglao_cells.copy(), None, "panglao")

# # GO
# print("Getting GO terms...")
# go_terms, go_test_cases, go_graph = get_go_terms("testing/go_basic.obo", "biological_process", num_terms=5)
# with open("go/results_with_descriptions.pkl", 'rb') as f:
#     go_test_cases = pickle.load(f)
# results = validation_scheme(go_test_cases, go_terms, None, "go", go_graph)
# breakpoint()
# del go_test_cases
# del go_graph
# del results
# gc.collect()


# # MSigDB
# go_terms, go_test_cases, _ = get_go_terms("testing/go_basic.obo", "biological_process", num_terms=0)
# print("Getting MSigDB terms...")
# with open("testing/msigdb_hallmarks.gmt", "r") as f:
#     msigdb_hallmarks = [line.strip().split('\t') for line in f.readlines()]
#     msigdb_hallmarks = [{
#         "name": " ".join(term[0].split("_")[1:]).lower(),
#         "genes": term[2:]
#     } for term in msigdb_hallmarks]

# background = msigdb_hallmarks + go_terms
# with open("msigdb/results.pkl", 'rb') as f:
#     msigdb_test_cases = pickle.load(f)
# results = validation_scheme(msigdb_test_cases, background, "MSigDB-H", "msigdb")
# del background
# del go_terms
# del msigdb_hallmarks
# del results
# gc.collect()

# # Reactome
# print("Getting Reactome terms...")
# with open("testing/ReactomePathways.gmt", "r") as f:
#     reactome = [line.strip().split('\t') for line in f.readlines()]
#     reactome = [{
#         "name": term[0],
#         "genes": term[2:]
#     } for term in reactome if len(term[2:]) > 5]
# with open("reactome/results.pkl", 'rb') as f:
#     reactome_test_cases = pickle.load(f)
# results = validation_scheme(reactome_test_cases, reactome, "REAC", "reactome")

# # Nonsense
# test_cases = []

# with open("testing/ReactomePathways.gmt", "r") as f:
#     reactome_pathways = [line.strip().split('\t') for line in f.readlines()]
#     reactome_pathways = [{
#         "name": term[0],
#         "genes": term[2:]
#     } for term in reactome_pathways]

# background_genes = set()
# for term in reactome_pathways:
#     background_genes.update(term["genes"])
# background_genes = list(background_genes)

# with open("reactome/test_cases.json", "r") as f:
#     pure_test_cases = json.load(f)

# for test_case in random.sample(pure_test_cases, 200):
#     random_genes = random.sample(background_genes, len(test_case["genes"]))

#     # generate 50/50 mix case
#     mix_case = {
#         "name": test_case["name"] + " (50/50 mix)",
#         "genes": random.sample(test_case["genes"], (len(test_case["genes"])+1)//2) + random.sample(random_genes, len(test_case["genes"])//2)
#     }
#     test_cases.append(mix_case)

#     # generate random case
#     random_case = {
#         "name": test_case["name"] + " (random)",
#         "genes": random_genes
#     }
#     test_cases.append(random_case)

# with open("nonsense/results.pkl", 'rb') as f:
#     nonsense_test_cases = pickle.load(f)
# results = validation_scheme(nonsense_test_cases, reactome_pathways, "REAC", "nonsense")
# del test_cases
# del background_genes
# del results
# gc.collect()
    
# # Barcode ablation
# test_cases = random.sample(pure_test_cases, 200)
# with open("barcoding/results.pkl", 'rb') as f:
#     barcoding_test_cases = pickle.load(f)

# results = validation_scheme(barcoding_test_cases, reactome_pathways, "REAC", "barcoding")

# Shortlist
with open("testing/ReactomePathways.gmt", "r") as f:
    reactome_pathways = [line.strip().split('\t') for line in f.readlines()]
    reactome_pathways = [{
        "name": term[0],
        "genes": term[2:]
    } for term in reactome_pathways]

test_cases = []
with open("reactome/test_cases.json", "r") as f:
    original_test_cases = json.load(f)

random.seed(42)
random.shuffle(original_test_cases)
for test_case in original_test_cases:
    if len(test_cases) >= 200:
        break

    if 5<=len(test_case['genes'])<=500:
        subset = random.sample(test_case['genes'], 5)
        test_case['genes'] = subset
        test_cases.append(test_case)

with open("shortlist/results.pkl", 'rb') as f:
    shortlist_test_cases = pickle.load(f)
breakpoint()

results = validation_scheme(shortlist_test_cases, reactome_pathways, "REAC", "shortlist")