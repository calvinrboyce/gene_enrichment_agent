"""Gene enrichment analysis workflow."""

import os
import json
from typing import List, Dict, Any
from datetime import datetime
import re
import time

from src.enrichment_tools import ToppFunAnalyzer, GProfilerAnalyzer, EnrichrAnalyzer
from src.literature import LiteratureAnalyzer
from src.summarize import SummarizeAnalyzer

class GeneEnrichmentAgent:
    """Main agent for running gene enrichment analysis workflow."""
 
    def __init__(self,
                 open_ai_api_key: str,
                 entrez_api_key: str,
                 open_ai_model: str = "gpt-4.1-mini",
                 results_dir: str = "gene_enrichment_analysis_results",
                 enrichr_sources: Dict[str, str] = {},
                 gprofiler_sources: List[str] = [],
                 toppfun_sources: List[str] = ['ToppCell'],
                 terms_per_source: int = 20,
                 papers_per_gene: int = 2,
                 aggregate_papers: int = 15):
        """Initialize the agent with necessary tools and setup.
        Args:
            open_ai_api_key: OpenAI API key
            entrez_api_key: Entrez API key
            open_ai_model: OpenAI model to use
            results_dir: Directory to save results
            enrichr_sources: Dictionary of Enrichr sources to use
            gprofiler_sources: List of gProfiler sources to use
            toppfun_sources: List of ToppFun sources to use
            terms_per_source: Number of terms to retrieve per source
            papers_per_gene: Number of papers to retrieve for each gene
            aggregate_papers: Number of papers to retrieve mentioning subsets of genes
        """
        self.results_dir = results_dir
        self.terms_per_source = terms_per_source
        self.papers_per_gene = papers_per_gene
        self.aggregate_papers = aggregate_papers
        self.enrichr = EnrichrAnalyzer(enrichr_sources)
        self.gprofiler = GProfilerAnalyzer(gprofiler_sources)
        self.toppfun = ToppFunAnalyzer(toppfun_sources)
        self.literature = LiteratureAnalyzer(entrez_api_key, papers_per_gene, aggregate_papers)
        self.summarize = SummarizeAnalyzer(open_ai_api_key, open_ai_model)
    
    def _sanitize_analysis_name(self, analysis_name: str):
        """Sanitize symbols in analysis name"""
        # Remove bad symbols
        if analysis_name:
            analysis_name = analysis_name.replace('\ufeff', '').replace('\u200b', '')
            analysis_name = re.sub(r'[\\/:*?"<>|]', '_', analysis_name)
            analysis_name = analysis_name.strip(' .')
            analysis_name = re.sub(r'_+', '_', analysis_name)
        else:
            timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
            analysis_name = f"run_{timestamp}"
        
        return analysis_name

    def run_analysis(self,
                     genes: List[str],
                     email: str,
                     background_genes: List[str] = [],
                     ranked: bool=True,
                     search_terms: List[str] = [],
                     context: str = "None",
                     save_results: int = 0,
                     analysis_name: str = '',
                     holdout: str = None,
                     use_barcodes: bool = True):
        """Run the complete analysis workflow.
        Args:
            genes: List of gene symbols to analyze
            email: Email address for NCBI's reference
            background_genes: List of background genes to use for the enrichment analysis
            ranked: Whether the genes are ranked by differential expression
            search_terms: List of search terms to use in the literature search
            context: Context of where the genes came from and what you're studying
            save_results: 0 saves no results, 1 saves an excel file, 2 saves json files of all intermediate results
            analysis_name: Name of the analysis
            use_barcodes: If True, LLM identifies terms by barcode; if False, LLM reproduces
                term fields and hallucination rates are attached to the return value

        Returns:
            Themed results: A dictionary of themed results with the following keys:
                * themes: A list of themes with the following keys:
                    * theme: The name of the theme
                    * description: A brief description of the function of the theme and why you identified it
                    * confidence: A confidence score for the theme, between 0 and 1
                    * terms: A list of terms from the various tools that were used to identify the theme
                * summary: A summary of the results
                * hallucination_metrics: Present when use_barcodes is False
        """
        # Run enrichment analyses
        start = time.time()
        print("Running enrichment analyses...")
        enrichr_results = self.enrichr.analyze(genes, background_genes)
        toppfun_results = self.toppfun.analyze(genes)
        gprofiler_results = self.gprofiler.analyze(genes, background_genes)

        # Search literature
        print("Searching literature...")
        literature_results = self.literature.search_literature(genes, email, search_terms, ranked)

        # Fetch gene summaries
        if ranked or len(genes) <= 15:
            print("Fetching gene summaries...")
            gene_summaries = self.literature.fetch_gene_summaries(genes[:15])
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
                                                               holdout,
                                                               use_barcodes=use_barcodes)

        # Save results
        if save_results > 0:
            print("Saving results...")
            analysis_name = self._sanitize_analysis_name(analysis_name)
            os.makedirs(self.results_dir, exist_ok=True)

            metadata = {
                'genes': genes,
                'email': email,
                'search_terms': search_terms,
                'context': context,
                'date': datetime.now(),
                'runtime': time.time() - start,
                'open_ai_model': self.summarize.summarize_model,
                'enrichr_sources': self.enrichr.sources,
                'toppfun_sources': ['GO:BP', 'GO:MF', 'GO:CC', 'PATHWAY', 'PPI'] + self.toppfun.additional_sources,
                'gprofiler_sources': ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP'] + self.gprofiler.additional_sources,
                'terms_per_source': self.terms_per_source,
                'papers_per_gene': self.papers_per_gene,
                'aggregate_papers': self.aggregate_papers,
                'background_genes': background_genes,
                'ranked': ranked,
            }
            self.summarize.synthesize_analysis(themed_results, analysis_name, self.results_dir, metadata)
        
        if save_results > 1:
            self.summarize.save_intermediate_results(metadata,
                                                     enrichr_results,
                                                     toppfun_results,
                                                     gprofiler_results,
                                                     literature_results,
                                                     analysis_name,
                                                     self.results_dir)

        return themed_results