import random
from collections import defaultdict
import pronto
import networkx as nx

def get_go_terms(ontology_file, namespace, num_terms=100, min_genes=5, seed=42):
    # Step 1: Load GO Ontology
    ontology = pronto.Ontology(ontology_file)
    G = nx.Graph()

    # Step 2: Load DIRECT gene annotations from GAF
    direct_go_to_genes = defaultdict(set)

    with open("testing/goa_human.gaf") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.strip().split("\t")
            gene = parts[2]
            go_id = parts[4]
            aspect = parts[8]
            
            if aspect == "P":  # Only biological process
                direct_go_to_genes[go_id].add(gene)
    
    # Step 3: Filter terms and propagate annotations
    terms = []
    propagated_go_to_genes = defaultdict(set)
    
    for term in ontology.terms():
        # Only process terms in your target namespace
        if term.id.startswith("GO:") and term.namespace == namespace:
            terms.append({
                'id': term.id,
                'name': term.name,
                'definition': term.definition
            })
            
            # The True Path Rule: Aggregate genes from this term AND all descendants.
            # term.subclasses() yields the term itself and all child terms recursively.
            for child in term.subclasses():
                if child.id in direct_go_to_genes:
                    propagated_go_to_genes[term.id].update(direct_go_to_genes[child.id])
            
            # Build graph
            G.add_node(term.id)
            for parent in term.superclasses(distance=1):
                if parent.id != term.id and parent.id.startswith("GO:") and parent.namespace == namespace:
                    G.add_edge(term.id, parent.id)
            for rel, linked_terms in term.relationships.items():
                for linked_term in linked_terms:
                    if linked_term.id.startswith("GO:") and linked_term.namespace == namespace:
                        G.add_edge(term.id, linked_term.id)

    # Step 4: Superset all background genes
    # We can just flatten the direct dictionary for the background
    all_genes = set(gene for genes in direct_go_to_genes.values() for gene in genes)
    
    with open("go/background_genes.txt", "w") as f:
        for gene in all_genes:
            if not gene.startswith("URS00"):
                f.write(gene + "\n")

    # Step 5: Prepare test cases using the PROPAGATED gene counts
    test_cases = []
    random.seed(seed)
    random.shuffle(terms)
    
    for term in terms:
        if len(test_cases) >= num_terms:
            break
            
        genes = propagated_go_to_genes.get(term['id'], set())
        
        # Check constraints against the fully propagated set
        if 500 >= len(genes) >= min_genes:
            test_cases.append({
                "id": term['id'],
                "name": term['name'],
                "definition": term['definition'],
                "genes": list(genes)
            })

    return terms, test_cases, G
