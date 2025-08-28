# from testing.go import get_go_terms
from testing.sapbert import embed_terms
from gene_enrichment_agent import GeneEnrichmentAgent
import os
import dotenv
from tqdm import tqdm
import numpy as np
import pickle
import json
import random
import time

def validation_scheme(test_cases, database_terms, holdout, validation_name):
    # we changed this function for each of the experiments to accomodate the different databases
    error_log = []

    # make directory for results
    os.makedirs(validation_name, exist_ok=True)

    with open(validation_name + "/test_cases.json", "w") as f:
        json.dump(test_cases, f)

    # embed database terms
    print("Embedding database terms")
    db_name_to_index = {term["name"]: i for i, term in enumerate(database_terms)}
    with open(validation_name + "/name_to_index.pkl", "wb") as f:
        pickle.dump(db_name_to_index, f)
    db_names = ["Name: " + term["name"] for term in database_terms]
    # db_names = ["Name: " + term["name"] + ". Description: " + term["definition"] for term in database_terms]
    with open(validation_name + "/db_names.pkl", "wb") as f:
        pickle.dump(db_names, f)
    db_embeddings = embed_terms(db_names)

    with open(validation_name + "/db_embeddings.pkl", "wb") as f:
        pickle.dump(db_embeddings, f)

    # run enrichment_agent on test_cases
    print("Running enrichment agent on test cases")
    dotenv.load_dotenv()
    open_ai_api_key = os.getenv("OPENAI_API_KEY")

    # for msigdb
    # enrichr_sources = {
    #     "CellMarker_2024": "CellMarker",
    #     "ChEA_2022": "ChEA"
    # }

    gea = GeneEnrichmentAgent(open_ai_api_key, num_papers=20)
    email = "cboyce3@mgh.harvard.edu"
    i = -1
    for test_case in tqdm(test_cases):
        i += 1
        try:
            start = time.time()
            # context = "Include a theme for likely cell type"
            test_case["enrichment_results"] = gea.run_analysis(test_case["genes"], email, ranked=False, save_results=False, holdout=holdout)
            test_case["theme_embeddings"] = embed_terms(["Name: " + theme["theme"] + ". Description: " + theme["description"] for theme in test_case["enrichment_results"]["themes"]])
            test_case["inference_time"] = time.time() - start
        except Exception as e:
            error_log.append(test_case["name"])
            continue
        # save the test cases
        if not i % 10:
            try:
                with open(validation_name + "/results.pkl", "wb") as f:
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
        if test_case["name"] in error_log:
            continue
        ys = db_embeddings
        if test_case["name"].endswith(" (50/50 mix)") or test_case["name"].endswith(" (random)"):
            y_index = db_name_to_index[test_case["name"].rsplit(" (", 1)[0]]
        else:
            y_index = db_name_to_index[test_case["name"]]
        y = ys[y_index]
        yhats = test_case["theme_embeddings"]

        # find the theme most similar to the correct GO term
        cosine_similarities = np.dot(yhats, y) / (np.linalg.norm(y) * np.linalg.norm(yhats, axis=1))
        yhat_index = np.argmax(cosine_similarities)
        yhat = yhats[yhat_index]
        test_case["most_similar_theme"] = test_case["enrichment_results"]["themes"][yhat_index]

        # find the percentile of the actual similarity in the list of cosine similarities
        test_case["cosine_similarities"] = np.dot(ys, yhat) / (np.linalg.norm(ys, axis=1) * np.linalg.norm(yhat))
        actual_similarity = test_case["cosine_similarities"][y_index]
        percentile = (np.sum(test_case["cosine_similarities"] <= actual_similarity) - 1) / (len(test_case["cosine_similarities"]) - 1)
        test_case["percentile"] = percentile

        # save as we go
        if not i % 10:
            try:
                with open(validation_name + "/results.pkl", "wb") as f:
                    pickle.dump(test_cases, f)
            except Exception as e:
                print("Cosine similarity pickling failed at index " + str(i))
                return

    # save the test cases
    try:
        print("Saving test cases")
        with open(validation_name + "/results.pkl", "wb") as f:
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

# results = validation_scheme(panglao_cells, panglao_cells.copy(), None, "panglao")

# # GO
# go_terms, go_test_cases = get_go_terms("testing/go_basic.obo", "biological_process", num_terms=500)
# results = validation_scheme(go_test_cases, go_terms, None, "go")

# # MSigDB
# go_terms, go_test_cases = get_go_terms("testing/go_basic.obo", "biological_process", num_terms=0)
# with open("testing/msigdb_hallmarks.gmt", "r") as f:
#     msigdb_hallmarks = [line.strip().split('\t') for line in f.readlines()]
#     msigdb_hallmarks = [{
#         "name": " ".join(term[0].split("_")[1:]).lower(),
#         "genes": term[2:]
#     } for term in msigdb_hallmarks]

# background = msigdb_hallmarks + go_terms

# results = validation_scheme(msigdb_hallmarks, background, 26771021, "msigdb")

# # Reactome
# with open("testing/ReactomePathways.gmt", "r") as f:
#     reactome = [line.strip().split('\t') for line in f.readlines()]
#     reactome = [{
#         "name": term[0],
#         "genes": term[2:]
#     } for term in reactome]

# results = validation_scheme(random.sample(reactome, 500), reactome, "REAC", "reactome")

# Nonsense
test_cases = []
database_terms = []

with open("testing/ReactomePathways.gmt", "r") as f:
    reactome_pathways = [line.strip().split('\t') for line in f.readlines()]
    reactome_pathways = [{
        "name": term[0],
        "genes": term[2:]
    } for term in reactome_pathways]

background_genes = set()
for term in reactome_pathways:
    background_genes.update(term["genes"])
background_genes = list(background_genes)

with open("reactome/test_cases.json", "r") as f:
    pure_test_cases = json.load(f)

for test_case in random.sample(pure_test_cases, 200):
    random_genes = random.sample(background_genes, len(test_case["genes"]))

    # generate 50/50 mix case
    mix_case = {
        "name": test_case["name"] + " (50/50 mix)",
        "genes": random.sample(test_case["genes"], (len(test_case["genes"])+1)//2) + random.sample(random_genes, len(test_case["genes"])//2)
    }
    test_cases.append(mix_case)

    # generate random case
    random_case = {
        "name": test_case["name"] + " (random)",
        "genes": random_genes
    }
    test_cases.append(random_case)

results = validation_scheme(test_cases, reactome_pathways, "REAC", "nonsense")
    

