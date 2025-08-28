"""gProfiler enrichment analysis tool integration."""

from collections import defaultdict
from typing import List, Dict, Any
from gprofiler import GProfiler

class GProfilerAnalyzer:
    """Handler for gProfiler enrichment analysis."""
    
    def __init__(self, additional_sources: List[str] = []):
        """Initialize the gProfiler analyzer.
        Args:
            additional_sources: List of gProfiler sources to use in addition to the default sources
        """
        self.additional_sources = additional_sources

    def analyze(self, genes: List[str], background_genes: List[str] = []) -> Dict[str, Any]:
        """Run gProfiler enrichment analysis and organize results by source.
        
        Args:
            genes: List of gene symbols to analyze
            background_genes: List of background genes to use for the enrichment analysis
        Returns:
            Dict containing:
                - results: Organized results by source (GO:BP, GO:MF, GO:CC, etc.)
                - summary_stats: Summary statistics of the analysis
                - raw_results: Raw gProfiler API response
                
        Raises:
            ValueError: If the gene list is empty or API call fails
        """
        if not genes:
            raise ValueError("Gene list cannot be empty")
            
        raw_results = self._run_query(genes, background_genes)
        organized_results = self._process_results(raw_results)

        return organized_results

    def _run_query(self, genes: List[str], background_genes: List[str] = []) -> List[Dict[str, Any]]:
        """Execute gProfiler enrichment query."""
        try:
            gp = GProfiler()
            if background_genes:
                results = gp.profile(query=genes, background=background_genes)
            else:
                results = gp.profile(query=genes)
            if not results:
                print("No results returned from gProfiler")
                return []
            return results
        except Exception as e:
            raise ValueError(f"Error running gProfiler analysis: {str(e)}")

    def _process_results(self, raw_results: List[Dict[str, Any]]) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """Process and organize gProfiler results by source.
        
        Returns:
            Dict where the keys are categories
            Each category has a dictionary where the keys are the IDs and the values are the results
        """
        all_categories = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP'] + self.additional_sources
        organized_results = defaultdict(dict)
        
        for result in raw_results:
            source = result.get('source')
            if not source or source not in all_categories:
                continue

            # We are going to be consolidating the results from all three tools,
            # so the pathways need to use a sanitized version of the term name as the ID.
            pathways = ['KEGG', 'REAC', 'WP']
            if source in pathways:
                ID = result['name'].lower()
            else:
                ID = result['native']
            
            cleaned_result = {
                'id': result['native'],
                'name': result['name'],
                'gprofiler_p_value': result['p_value'],
                'term_size': result['term_size'],
                'description': result['description']
            }

            organized_results[source][ID] = cleaned_result
            
        return organized_results
