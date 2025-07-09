"""Gene enrichment analysis workflow."""

import os
import json
from typing import List, Dict, Any
from datetime import datetime
import re

from src.enrichment_tools import ToppFunAnalyzer, GProfilerAnalyzer, EnrichrAnalyzer
from src.literature import LiteratureAnalyzer
from src.summarize import SummarizeAnalyzer

class GeneEnrichmentAgent:
    """Main agent for running gene enrichment analysis workflow."""
 
    def __init__(self,
                 open_ai_api_key: str,
                 open_ai_model: str = "gpt-4.1-mini",
                 results_dir: str = "gene_enrichment_analysis_results",
                 enrichr_sources: Dict[str, str] = {},
                 gprofiler_sources: Dict[str, str] = {},
                 toppfun_sources: Dict[str, str] = {},
                 terms_per_source: int = 15,
                 papers_per_gene: int = 2,
                 max_papers: int = 10):
        """Initialize the agent with necessary tools and setup.
        Args:
            open_ai_api_key: OpenAI API key
            open_ai_model: OpenAI model to use
            results_dir: Directory to save results
            enrichr_sources: Dictionary of Enrichr sources to use
            gprofiler_sources: Dictionary of gProfiler sources to use
            toppfun_sources: Dictionary of ToppFun sources to use
            terms_per_source: Number of terms to retrieve per source
            papers_per_gene: Number of papers to retrieve per gene
            max_papers: Maximum number of papers to retrieve for full list
        """
        self.results_dir = results_dir
        self.enrichr = EnrichrAnalyzer(enrichr_sources, terms_per_source)
        self.gprofiler = GProfilerAnalyzer(gprofiler_sources, terms_per_source)
        self.toppfun = ToppFunAnalyzer(toppfun_sources, terms_per_source)
        self.literature = LiteratureAnalyzer(papers_per_gene, max_papers)
        self.summarize = SummarizeAnalyzer(open_ai_api_key, open_ai_model)

    def _save_results(self,
                      genes: List[str],
                      ranked: bool,
                      email: str,
                      search_terms: List[str],
                      context: str,
                      enrichr_results: Dict[str, Any],
                      toppfun_results: Dict[str, Any],
                      gprofiler_results: Dict[str, Any],
                      literature_results: Dict[str, Any],
                      themed_results: Dict[str, Any],
                      analysis_name: str) -> str:
        """Save analysis results to file.
        
        Args:
            genes: List of gene symbols to analyze
            ranked: Whether the genes are ranked by differential expression
            email: Email address for NCBI's reference
            search_terms: List of search terms to use in the literature search
            context: Context of where the genes came from and what you're studying
            enrichr_results: Enrichr analysis results
            toppfun_results: ToppFun analysis results
            gprofiler_results: gProfiler analysis results
            literature_results: Literature search results
            themed_results: Themed results
            analysis_name: Name of the analysis
        Returns:
            Path to the saved results directory
        """
        # Create timestamped directory for this run
        timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        safe_analysis_name = analysis_name if analysis_name else f"run_{timestamp}"
        run_dir = os.path.join(self.results_dir, safe_analysis_name)
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(run_dir, exist_ok=True)

        # Helper to prefix filenames
        def prefix_filename(base):
            return f"{run_dir}/{safe_analysis_name}_{base}.json"

        # Save input parameters
        with open(prefix_filename("input_params"), 'w') as f:
            json.dump({
                "genes": genes,
                "ranked": ranked,
                "email": email,
                "search_terms": search_terms,
                "context": context,
                "timestamp": timestamp
            }, f, indent=2)

        # Save enrichment results
        with open(prefix_filename("enrichr_results"), 'w') as f:
            json.dump(enrichr_results, f, indent=2)
        with open(prefix_filename("toppfun_results"), 'w') as f:
            json.dump(toppfun_results, f, indent=2)
        with open(prefix_filename("gprofiler_results"), 'w') as f:
            json.dump(gprofiler_results, f, indent=2)

        # Save literature results
        with open(prefix_filename("literature_results"), 'w') as f:
            json.dump(literature_results, f, indent=2)

        # Save themed results
        with open(prefix_filename("themed_results"), 'w') as f:
            json.dump(themed_results, f, indent=2)

        # Generate final analysis
        self.summarize.synthesize_analysis(themed_results, genes, email, search_terms, context, analysis_name, run_dir)

        return run_dir

    def run_analysis(self,
                     genes: List[str],
                     email: str,
                     ranked: bool=True,
                     search_terms: List[str] = [],
                     context: str = "None",
                     save_results: bool = True,
                     analysis_name: str = None):
        """Run the complete analysis workflow.
        Args:
            genes: List of gene symbols to analyze
            email: Email address for NCBI's reference
            ranked: Whether the genes are ranked by differential expression
            search_terms: List of search terms to use in the literature search
            context: Context of where the genes came from and what you're studying
            save_results: Whether to save the results to file
            analysis_name: Name of the analysis

        Returns:
            Themed results: A dictionary of themed results with the following keys:
                * themes: A list of themes with the following keys:
                    * theme: The name of the theme
                    * description: A brief description of the function of the theme and why you identified it
                    * confidence: A confidence score for the theme, between 0 and 1
                    * terms: A list of terms from the various tools that were used to identify the theme
                * summary: A summary of the results
        """
        # Sanitize analysis name
        if analysis_name:
            analysis_name = re.sub(r'[ \\/:*?"<>|\'`~!@#$%^&\(\)]', '_', analysis_name)
        
        # Run enrichment analyses
        print("Running enrichment analyses...")
        enrichr_results = self.enrichr.analyze(genes)
        toppfun_results = self.toppfun.analyze(genes)
        gprofiler_results = self.gprofiler.analyze(genes)

        # Search literature
        print("Searching literature...")
        if ranked:
            queries = genes[:10] + [genes]
        else:
            queries = [genes]
        literature_results = self.literature.search_literature(queries, email, search_terms)

        # Group results by theme
        print("Grouping results by theme...")
        themed_results = self.summarize.group_results_by_theme(enrichr_results,
                                                               toppfun_results,
                                                               gprofiler_results,
                                                               literature_results,
                                                               genes,
                                                               ranked,
                                                               context)

        # Save results
        if save_results:
            print("Saving results...")
            self._save_results(genes,
                               ranked,
                               email,
                               search_terms,
                               context,
                               enrichr_results,
                               toppfun_results,
                               gprofiler_results,
                               literature_results,
                               themed_results,
                               analysis_name)

        return themed_results