"""gProfiler enrichment analysis tool integration."""

from typing import List, Dict, Any
from gprofiler import GProfiler

class GProfilerAnalyzer:
    """Handler for gProfiler enrichment analysis."""
    
    def __init__(self, sources: Dict[str, str] = {}, terms_per_source: int = 10):
        """Initialize the gProfiler analyzer.
        Args:
            sources: Dictionary of gProfiler sources to use
            terms_per_source: Number of terms to retrieve per source
        """
        self.terms_per_source = terms_per_source
        self.sources = {
            "Biological Process": "GO:BP",
            "Molecular Function": "GO:MF",
            "Cellular Component": "GO:CC",
            "KEGG pathways": "KEGG",
            "Reactome": "REAC",
            "Transcription Factors": "TF",
            "Protein Complexes": "CORUM",
            "Human Phenotype": "HP",
            "Human Protein Atlas": "HPA",
            "WikiPathways": "WP",
            "miRNA targets": "MIRN"
        } if not sources else sources

        self.shortlist_categories = ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'TF']

    def analyze(self, genes: List[str]) -> Dict[str, Any]:
        """Run gProfiler enrichment analysis and organize results by source.
        
        Args:
            genes: List of gene symbols to analyze
            
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
            
        raw_results = self._run_query(genes)
        organized_results = self._process_results(raw_results)
        shortlist = self._generate_shortlist(organized_results)
        
        return {
            "results": organized_results,
            "shortlist": shortlist
        }

    def _run_query(self, genes: List[str]) -> List[Dict[str, Any]]:
        """Execute gProfiler enrichment query."""
        try:
            gp = GProfiler()
            results = gp.profile(query=genes)
            if not results:
                raise ValueError("No results returned from gProfiler")
            return results
        except Exception as e:
            raise ValueError(f"Error running gProfiler analysis: {str(e)}")

    def _process_results(self, raw_results: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
        """Process and organize gProfiler results by source."""
        organized_results = {abbrev: [] for abbrev in self.sources.values()}
        
        for result in raw_results:
            source = result.get('source')
            if source in organized_results:
                organized_results[source].append(result)
        
        for source in organized_results:
            organized_results[source].sort(key=lambda x: x['p_value'])
            
        return organized_results

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
                        'description': result['description'],
                        'p_value': result['p_value'],
                        'source': category,
                        'tool': 'gprofiler'
                    }
                )
        return shortlist