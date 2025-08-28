import pickle
import matplotlib.pyplot as plt
import json
from testing.sapbert import embed_terms
from gene_enrichment_agent import GeneEnrichmentAgent
import os
import time
import dotenv
import numpy as np
from tqdm import tqdm

def fill_errors(results_dir):
    # we changed this function for each of the experiments to accomodate the different databases
    with open(results_dir + "/error_log1.txt", "r") as f:
        error_log = f.readlines()
    error_log = [error.strip() for error in error_log]
    new_error_log = []

    # load test cases
    with open(results_dir + "/results1.pkl", "rb") as f:
        test_cases = pickle.load(f)

    # embed database terms
    with open(results_dir + "/name_to_index.pkl", "rb") as f:
        db_name_to_index = pickle.load(f)
    with open(results_dir + "/db_embeddings.pkl", "rb") as f:
        db_embeddings = pickle.load(f)

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
    for test_case in tqdm(test_cases):
        if test_case["name"] not in error_log:
            continue
        try:
            start = time.time()
            # context = "Include a theme for likely cell type"
            test_case["enrichment_results"] = gea.run_analysis(test_case["genes"], email, ranked=False, save_results=False, holdout='REAC')
            test_case["theme_embeddings"] = embed_terms(["Name: " + theme["theme"] + ". Description: " + theme["description"] for theme in test_case["enrichment_results"]["themes"]])
            test_case["inference_time"] = time.time() - start
        except Exception as e:
            print(e)
            new_error_log.append(test_case["name"])
            continue
    
    # save the error log
    print("Saving error log")
    with open(results_dir + "/error_log2.txt", "w") as f:
        for error in new_error_log:
            f.write(error + "\n")
    
    # compute cosine similarities
    print("Computing cosine similarities and percentiles")
    for test_case in tqdm(test_cases):
        if test_case["name"] in new_error_log:
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
        percentile = (np.sum(test_case["cosine_similarities"] <= actual_similarity)-1) / (len(test_case["cosine_similarities"]) -1)
        test_case["percentile"] = percentile

    # save the test cases
    try:
        print("Saving test cases")
        with open(results_dir + "/results1.pkl", "wb") as f:
            pickle.dump(test_cases, f)
    except Exception as e:
        print("Final pickling failed")
        return
    
    return test_cases

# fill_errors("nonsense")

with open("nonsense/name_to_index.pkl", "rb") as f:
    name_to_index = pickle.load(f)

print(name_to_index.keys())
