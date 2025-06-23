"""ToppFun enrichment analysis tool integration."""

from typing import List, Dict, Any
import requests

class ToppFunAnalyzer:
    """Handler for ToppFun enrichment analysis."""
    
    def __init__(self, sources: Dict[str, str] = {}, terms_per_source: int = 10):
        """Initialize the ToppFun analyzer.
        Args:
            sources: Dictionary of ToppFun sources to use
            terms_per_source: Number of terms to retrieve per source
        """
        self.terms_per_source = terms_per_source
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
            'ToppCell': 'CELL',
            'Computational': 'COMP',
            'MicroRNA': 'MIRNA',
            'Drug': 'DRUG',
            'Disease': 'DISEASE'
        } if not sources else sources

        self.shortlist_categories = ['GO:BP', 'GO:MF', 'GO:CC', 'TFBS', 'CELL']

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
        shortlist = self._generate_shortlist(organized_results)
        
        return {
            "results": organized_results,
            "shortlist": shortlist
        }

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

    def _process_results(self, raw_results: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
        """Process and organize ToppFun results by category."""
        organized_results = {v: [] for v in self.category_mapping.values()}
        
        for result in raw_results:
            if not isinstance(result, dict) or 'Category' not in result:
                continue
                
            category = self.category_mapping.get(result['Category'])
            if not category:
                continue
                
            try:
                cleaned_result = self._clean_result(result)
                organized_results[category].append(cleaned_result)
            except (ValueError, TypeError, KeyError) as e:
                print(f"Warning: Error processing result: {e}")
                continue
        
        for category in organized_results:
            organized_results[category].sort(key=lambda x: x['p_value'])
            
        return organized_results

    def _clean_result(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Clean and transform a single ToppFun result."""
        cleaned_result = {
            'name': result.get('Name', ''),
            'id': result.get('ID', ''),
            'p_value': float(result.get('PValue', 1.0)),
            'q_value': float(result.get('QValueFDRBH', 1.0)),
            'total_genes': int(result.get('TotalGenes', 0)),
            'genes_in_term': int(result.get('GenesInTerm', 0)),
            'genes_in_query': int(result.get('GenesInQuery', 0)),
            'intersection_size': int(result.get('GenesInTermInQuery', 0)),
            'genes': [gene.get('Symbol', '') for gene in result.get('Genes', [])],
            'url': result.get('URL', '').strip() if result.get('URL', '').strip() != ' ' else None
        }

        return cleaned_result
    
    def _generate_shortlist(self, organized_results: Dict[str, List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
        """Generate a shortlist of the names, sources, and p-values of the top 10 results for each category."""
        shortlist = []
        for category in self.shortlist_categories:
            if not organized_results[category]:
                continue
            for result in organized_results[category][:self.terms_per_source]:
                shortlist.append(
                    {
                        'name': result['name'],
                        'genes': result['genes'],
                        'p_value': result['p_value'],
                        'source': category,
                        'tool': 'toppfun'
                    }
                )
        return shortlist
