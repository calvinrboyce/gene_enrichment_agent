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
                 gprofiler_sources: List[str] = [],
                 toppfun_sources: List[str] = ['ToppCell'],
                 terms_per_source: int = 20,
                 num_papers: int = 10):
        """Initialize the agent with necessary tools and setup.
        Args:
            open_ai_api_key: OpenAI API key
            open_ai_model: OpenAI model to use
            results_dir: Directory to save results
            enrichr_sources: Dictionary of Enrichr sources to use
            gprofiler_sources: List of gProfiler sources to use
            toppfun_sources: List of ToppFun sources to use
            terms_per_source: Number of terms to retrieve per source
            num_papers: Number of papers to retrieve from PubMed
        """
        self.results_dir = results_dir
        self.terms_per_source = terms_per_source
        self.enrichr = EnrichrAnalyzer(enrichr_sources)
        self.gprofiler = GProfilerAnalyzer(gprofiler_sources)
        self.toppfun = ToppFunAnalyzer(toppfun_sources)
        self.literature = LiteratureAnalyzer(num_papers)
        self.summarize = SummarizeAnalyzer(open_ai_api_key, open_ai_model)

    def run_analysis(self,
                     genes: List[str],
                     email: str,
                     background_genes: List[str] = [],
                     ranked: bool=True,
                     search_terms: List[str] = [],
                     context: str = "None",
                     save_results: bool = True,
                     analysis_name: str = None,
                     holdout: str = None):
        """Run the complete analysis workflow.
        Args:
            genes: List of gene symbols to analyze
            email: Email address for NCBI's reference
            background_genes: List of background genes to use for the enrichment analysis
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
        # Run enrichment analyses
        print("Running enrichment analyses...")
        enrichr_results = self.enrichr.analyze(genes, background_genes)
        toppfun_results = self.toppfun.analyze(genes)
        gprofiler_results = self.gprofiler.analyze(genes, background_genes)

        # Search literature
        print("Searching literature...")
        literature_results = self.literature.search_literature(genes, email, search_terms)

        # Fetch gene summaries
        if ranked or len(genes) <= 10:
            print("Fetching gene summaries...")
            gene_summaries = self.literature.fetch_gene_summaries(genes[:10])
        else:
            gene_summaries = []

        # Group results by theme
        print("Grouping results by theme...")
        themed_results = self.summarize.group_results_by_theme(enrichr_results,
                                                               toppfun_results,
                                                               gprofiler_results,
                                                               literature_results,
                                                               gene_summaries,
                                                               self.terms_per_source,
                                                               genes,
                                                               ranked,
                                                               context,
                                                               holdout)

        # Save results
        if save_results:
            print("Saving results...")
            if analysis_name:
                analysis_name = analysis_name.replace('\ufeff', '').replace('\u200b', '')
                analysis_name = re.sub(r'[\\/:*?"<>|]', '_', analysis_name)
                analysis_name = analysis_name.strip(' .')
                analysis_name = re.sub(r'_+', '_', analysis_name)
            else:
                timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
                analysis_name = f"run_{timestamp}"
            os.makedirs(self.results_dir, exist_ok=True)
            self.summarize.synthesize_analysis(themed_results, genes, email, search_terms, context, analysis_name, self.results_dir)

        return themed_results