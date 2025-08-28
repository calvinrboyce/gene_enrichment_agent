import random
from collections import defaultdict
import pronto

def get_go_terms(ontology, namespace, num_terms=100, min_genes=5, seed=42):
    # Step 1: Load GO Ontology
    ontology = pronto.Ontology(ontology)

    # Step 2: Filter to GO:BP terms
    terms = [{
        'id': term.id,
        'name': term.name,
        'definition': term.definition
    } for term in ontology.terms() if term.id.startswith("GO:") and term.namespace == namespace]

    # Step 4: Load gene annotations from GAF
    go_to_genes = defaultdict(set)

    with open("testing/goa_human.gaf") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.strip().split("\t")
            gene = parts[2]
            go_id = parts[4]
            aspect = parts[8]
            if aspect == "P":  # Only biological process
                go_to_genes[go_id].add(gene)
    
    # superset all genes
    all_genes = set()
    for genes in go_to_genes.values():
        all_genes.update(genes)
    with open("go/background_genes.txt", "w") as f:
        for gene in all_genes:
            if gene.startswith("URS00"):
                continue
            f.write(gene + "\n")

    # Step 5: Prepare test cases
    test_cases = []

    random.seed(seed)
    random.shuffle(terms)
    for term in terms:
        if len(test_cases) >= num_terms:
            break
        genes = go_to_genes.get(term['id'])
        if genes and 500 >= len(genes) >= min_genes:
            test_cases.append({
                "id": term['id'],
                "name": term['name'],
                "definition": term['definition'],
                "genes": list(genes)
            })

    return terms, test_cases
