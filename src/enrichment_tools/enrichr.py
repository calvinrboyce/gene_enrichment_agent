"""Enrichr enrichment analysis tool integration."""

from typing import List, Dict, Any
import requests

class EnrichrAnalyzer:
    """Handler for Enrichr enrichment analysis."""
    
    def __init__(self, sources: Dict[str, str] = {}, terms_per_source: int = 10):
        """Initialize the Enrichr analyzer.
        Args:
            sources: Dictionary of Enrichr sources to use
            terms_per_source: Number of terms to retrieve per source
        """
        self.terms_per_source = terms_per_source
        self.base_url = "https://maayanlab.cloud/Enrichr"
        self.background_base_url = "https://maayanlab.cloud/speedrichr/api"
        self.sources = {
            "GO_Biological_Process_2025": "GO:BP",
            "GO_Molecular_Function_2025": "GO:MF",
            "GO_Cellular_Component_2025": "GO:CC",
            "KEGG_2021_Human": "KEGG",
            "MSigDB_Hallmark_2020": "MSIGDB",
            "Reactome_Pathways_2024": "REACTOME",
            "TRANSFAC_and_JASPAR_PWMs": "TRANSFAC",
            "Allen_Brain_Atlas_10x_scRNA_2021": "ALLEN",
            "GTEx_Tissue_Expression_Up": "GTEX"
        } if not sources else sources

    def analyze(self, genes: List[str], background_genes: List[str] = []) -> Dict[str, Any]:
        """Run Enrichr enrichment analysis and organize results by source.
        
        Args:
            genes: List of gene symbols to analyze
            background_genes (optional): List of background genes to use for the enrichment analysis
        Returns:
            Dict containing:
                - results: Organized results by source (GO:BP, GO:MF, GO:CC, KEGG)
                - summary_stats: Summary statistics of the analysis
                - raw_results: Raw Enrichr API response
                
        Raises:
            ValueError: If the gene list is empty or API call fails
            requests.RequestException: If there is an error communicating with the API
        """
        # If background genes are provided, use the background base URL
        with_background = len(background_genes) > 0
        if not genes:
            raise ValueError("Gene list cannot be empty")
            
        # Add gene list to Enrichr
        user_list_id = self._upload_gene_list(genes, with_background)
        if with_background:
            background_id = self._upload_background_gene_list(background_genes)
        
        # Run enrichment analysis for each source
        organized_results = {}
        
        for source_name, source_key in self.sources.items():
            try:
                if with_background:
                    # print(f'Running enrichment analysis for {source_name} with background')
                    source_results = self._run_background_enrichment(user_list_id, background_id, source_name)
                else:
                    source_results = self._run_enrichment(user_list_id, source_name)
                organized_results[source_key] = self._process_results(source_results)
            except (ValueError, requests.RequestException) as e:
                print(f"Warning: Error processing {source_key} results: {e}")
                organized_results[source_key] = []
        
        # Generate shortlist
        shortlist = self._generate_shortlist(organized_results)
        
        return {
            "results": organized_results,
            "shortlist": shortlist
        }

    def _upload_gene_list(self, genes: List[str], with_background: bool = False) -> str:
        """Upload gene list to Enrichr.
        
        Args:
            genes: List of gene symbols
            with_background: Whether to use the background base URL
        Returns:
            User list ID from Enrichr
            
        Raises:
            ValueError: If the upload fails
            requests.RequestException: If there is an error communicating with the API
        """
        if with_background:
            base_url = self.background_base_url
        else:
            base_url = self.base_url
        try:
            response = requests.post(
                f"{base_url}/addList",
                files={'list': (None, '\n'.join(genes)), 'description': (None, 'Gene list')},
                timeout=30
            )
            response.raise_for_status()
            result = response.json()
            
            if 'userListId' not in result:
                raise ValueError("Invalid response from Enrichr upload API")
                
            return result['userListId']
            
        except requests.RequestException as e:
            raise ValueError(f"Error uploading gene list to Enrichr: {str(e)}")
        
    def _upload_background_gene_list(self, genes: List[str]) -> str:
        """Upload background gene list to Enrichr.
        
        Args:
            genes: List of gene symbols
        Returns:
            Background ID from Enrichr
            
        Raises:
            ValueError: If the upload fails
            requests.RequestException: If there is an error communicating with the API
        """
        try:
            response = requests.post(
                f"{self.background_base_url}/addbackground",
                data={'background': '\n'.join(genes)},
                timeout=30
            )
            response.raise_for_status()
            result = response.json()
            
            if 'backgroundid' not in result:
                raise ValueError("Invalid response from Enrichr background upload API")
                
            return result['backgroundid']
            
        except requests.RequestException as e:
            raise ValueError(f"Error uploading background gene list to Enrichr: {str(e)}")

    def _run_enrichment(self, user_list_id: str, source: str) -> List[List]:
        """Run enrichment analysis for a specific source.
        
        Args:
            user_list_id: Enrichr user list ID
            source: Name of the gene set library to query
            
        Returns:
            List of enrichment results
            
        Raises:
            ValueError: If the API returns invalid data
            requests.RequestException: If there is an error communicating with the API
        """
        try:
            response = requests.get(
                f"{self.base_url}/enrich",
                params={'userListId': user_list_id, 'backgroundType': source},
                timeout=30
            )
            response.raise_for_status()
            result = response.json()
            
            if source not in result:
                raise ValueError(f"Invalid response from Enrichr API for source: {source}")
                
            return result[source]
            
        except requests.RequestException as e:
            raise ValueError(f"Error running Enrichr analysis: {str(e)}")
        
    def _run_background_enrichment(self, user_list_id: str, background_id: str, source: str) -> List[List]:
        """Run enrichment analysis for a specific source with background genes.
        
        Args:
            user_list_id: Enrichr user list ID
            background_id: Enrichr background ID
            source: Name of the gene set library to query
            
        Returns:
            List of enrichment results
            
        Raises:
            ValueError: If the API returns invalid data
            requests.RequestException: If there is an error communicating with the API
        """
        try:
            response = requests.post(
                f"{self.background_base_url}/backgroundenrich",
                data={'userListId': user_list_id, 'backgroundid': background_id, 'backgroundType': source},
                timeout=30
            )
            response.raise_for_status()
            result = response.json()
            
            if source not in result:
                raise ValueError(f"Invalid response from Enrichr API for source: {source}")
                
            return result[source]
            
        except requests.RequestException as e:
            raise ValueError(f"Error running Enrichr analysis: {str(e)}")

    def _process_results(self, raw_results: List[List]) -> List[Dict[str, Any]]:
        """Process and clean Enrichr results.
        
        Args:
            raw_results: Raw results from Enrichr API
            
        Returns:
            List of processed and cleaned results
        """
        processed_results = []
        
        for result in raw_results:
            try:
                # Enrichr results format:
                # [Rank, Term name, P-value, Odds ratio, Combined score, 
                #  Overlapping genes, Adjusted p-value, Old p-value, Old adjusted p-value]
                cleaned_result = {
                    'rank': int(result[0]),
                    'name': result[1],
                    'p_value': float(result[2]),
                    'odds_ratio': float(result[3]),
                    'combined_score': float(result[4]),
                    'genes': result[5],
                    'q_value': float(result[6]),
                    'old_p_value': float(result[7]),
                    'old_q_value': float(result[8])
                }
                
                # Update with cleaned result
                processed_results.append(cleaned_result)
                
            except (IndexError, ValueError, TypeError) as e:
                print(f"Warning: Error processing result: {e}")
                continue
                
        return processed_results
    
    def _generate_shortlist(self, organized_results: Dict[str, List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
        """Generate a shortlist of the names, sources, and p-values of the top 10 results for each category."""
        shortlist = []
        for category in organized_results:
            if not organized_results[category]:
                continue
            for result in organized_results[category][:self.terms_per_source]:
                shortlist.append(
                    {
                        'name': result['name'],
                        'genes': result['genes'],
                        'p_value': result['p_value'],
                        'source': category,
                        'tool': 'enrichr'
                    }
                )
        return shortlist