"""ToppFun enrichment analysis tool integration."""

from collections import defaultdict
from typing import List, Dict, Any
import requests

class ToppFunAnalyzer:
    """Handler for ToppFun enrichment analysis."""
    
    def __init__(self, additional_sources: List[str] = []):
        """Initialize the ToppFun analyzer.
        Args:
            additional_sources: List of ToppFun sources to use in addition to the default sources
        """
        self.category_mapping = {
            'GeneOntologyMolecularFunction': 'GO:MF',
            'GeneOntologyBiologicalProcess': 'GO:BP',
            'GeneOntologyCellularComponent': 'GO:CC',
            'HumanPheno': 'HP',
            'MousePheno': 'MP',
            'Domain': 'DOMAIN',
            'Pathway': 'PATHWAY',
            'Pubmed': 'PUBMED',
            'Interaction': 'PPI',
            'Cytoband': 'CYTOBAND',
            'TFBS': 'TFBS',
            'GeneFamily': 'GENE_FAM',
            'Coexpression': 'COEXP',
            'CoexpressionAtlas': 'COEXP_ATLAS',
            'ToppCell': 'ToppCell Atlas',
            'Computational': 'COMP',
            'MicroRNA': 'MIRNA',
            'Drug': 'DRUG',
            'Disease': 'DISEASE'
        }

        self.additional_sources = [self.category_mapping[source] for source in additional_sources]

    def analyze(self, genes: List[str]) -> Dict[str, Any]:
        """Run ToppFun enrichment analysis and organize results by category.
        
        Args:
            genes: List of gene symbols to analyze
            
        Returns:
            Dict containing:
                - results: Organized results by category
                - summary_stats: Summary statistics of the analysis
                - raw_results: Raw API response data
                
        Raises:
            ValueError: If the API response is invalid or gene lookup fails
            requests.RequestException: If there is an error communicating with the API
        """
        if not genes:
            raise ValueError("Gene list cannot be empty")
        
        entrez_ids = self._lookup_entrez_ids(genes)
        raw_results = self._run_enrichment(entrez_ids)
        organized_results = self._process_results(raw_results)
        
        return organized_results

    def _lookup_entrez_ids(self, genes: List[str]) -> List[int]:
        """Convert gene symbols to Entrez IDs using ToppGene API."""
        try:
            response = requests.post(
                "https://toppgene.cchmc.org/API/lookup", 
                json={'Symbols': genes},
                timeout=30
            )
            response.raise_for_status()
            gene_info = response.json()
            
            if not isinstance(gene_info, dict) or 'Genes' not in gene_info:
                raise ValueError("Invalid response format from ToppGene lookup API")
            
            entrez_ids = [gene['Entrez'] for gene in gene_info['Genes'] if 'Entrez' in gene]
            
            if not entrez_ids:
                raise ValueError("No valid Entrez IDs found for provided genes")
                
            return entrez_ids
            
        except requests.RequestException as e:
            raise ValueError(f"Error communicating with ToppGene API: {str(e)}")

    def _run_enrichment(self, entrez_ids: List[int]) -> List[Dict[str, Any]]:
        """Run ToppFun enrichment analysis with Entrez IDs."""
        try:
            response = requests.post(
                "https://toppgene.cchmc.org/API/enrich", 
                json={'Genes': entrez_ids},
                timeout=30
            )
            response.raise_for_status()
            result_data = response.json()
            
            if not isinstance(result_data, dict) or 'Annotations' not in result_data:
                raise ValueError("Invalid response format from ToppFun enrichment API")
                
            return result_data['Annotations']
            
        except requests.RequestException as e:
            raise ValueError(f"Error communicating with ToppFun API: {str(e)}")

    def _process_results(self, raw_results: List[Dict[str, Any]]) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """Process and organize ToppFun results by category.
        
        Returns:
            Dict where the keys are categories
            Each category has a dictionary where the keys are the IDs and the values are the results
        """
        all_categories = ['GO:BP', 'GO:MF', 'GO:CC', 'PATHWAY', 'PPI'] + self.additional_sources
        organized_results = defaultdict(dict)
        
        for result in raw_results:
            if not isinstance(result, dict) or 'Category' not in result:
                continue
                
            category = self.category_mapping.get(result['Category'])
            if not category or category not in all_categories:
                continue

            result['Genes'] = [gene['Symbol'] for gene in result['Genes']]

            # Okay things start to get tricky here.
            # For results in the ontologies and pathways, we are going to be consolidating the results
            # from all three tools. The ontologies have unique IDs to match between the tools, but the
            # pathways need to use a sanitized version of the term name as the ID. So depending on where
            # the term comes from, we'll add it to the results using the appropriate ID.
            pathways = ['KEGG Legacy Pathways', 'Reactome Pathways', 'WikiPathways']
            if category == 'PATHWAY':
                if result['Source'] in pathways:
                    # These pathway databases use "REACTOME_TERM_NAME" or "KEGG_TERM_NAME" as the name,
                    # so we need to replace underscores with spaces and skip the first word.
                    category = 'KEGG' if result['Source'] == 'KEGG Legacy Pathways' else 'REAC' if result['Source'] == 'Reactome Pathways' else 'WP'
                    ID = ' '.join(result['Name'].lower().split('_')[1:])
                    name = ID
                elif result['Source'] == 'KEGG Medicus Pathways':
                    # These pathway databases use "KEGG_MEDICUS_SOMETHINGELSE_TERM_NAME" as the name,
                    # so we need to skip the first three words.
                    category = 'KEGG'
                    ID = ' '.join(result['Name'].lower().split('_')[3:])
                    name = ID
                else:
                    continue
            elif category == 'PPI':
                ID = result['ID'].split(':')[1].lower()
                name = ID
            else:
                ID = result['ID']
                name = result['Name']

            cleaned_result = {
                'id': result['ID'],
                'name': name,
                'toppfun_p_value': result['QValueFDRBH'],
                'term_size': result['GenesInTerm'],
                'genes': result['Genes']
            }

            
            organized_results[category][ID] = cleaned_result
            
        return organized_results
