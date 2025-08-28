"""Literature search functionality using PubMed."""

import time
import xml.etree.ElementTree as ET
from typing import List, Dict, Any
from Bio import Entrez
import re
import requests
from http.client import IncompleteRead
import json

class LiteratureAnalyzer:
    """Class for searching scientific literature related to genes."""

    def __init__(self, num_papers: int = 20):
        """Initialize the literature analyzer."""
        self.num_papers = num_papers
        self.api_key = "c8712aecaa1a952c3d72968f66d3a8853e08"

    def _extract_text_from_element(self, element) -> str:
        """Extract clean text from an XML element, preserving important formatting.
        
        Args:
            element: XML element to process
            
        Returns:
            Clean text with minimal formatting
        """
        text_parts = []
        
        # Add the element's direct text if it exists
        if element.text:
            text_parts.append(element.text.strip())
            
        # Process child elements
        for child in element:
            # Handle italic text
            if child.tag == 'italic':
                if child.text:
                    text_parts.append(child.text.strip())
            # Handle other formatting we want to preserve
            elif child.tag in ['sup', 'sub']:
                if child.text:
                    text_parts.append(child.text.strip())
            # For other elements, recursively extract text
            else:
                text_parts.append(self._extract_text_from_element(child))
            
            # Add any tail text
            if child.tail:
                text_parts.append(child.tail.strip())
                
        return ' '.join(filter(None, text_parts))

    def _highlight_genes_in_text(self, text: str, genes: List[str]) -> tuple[str, List[str]]:
        """Highlight genes in text and track which ones were found.
        
        Args:
            text: Text to process
            genes: List of genes to look for
            
        Returns:
            Tuple of (processed text, list of found genes)
        """
        found_genes = []
        if not text:
            return text, found_genes
            
        processed_text = f" {text} "  # Add spaces for better matching
        
        pattern = re.compile(rf"(?<!\w)({'|'.join(map(re.escape, genes))})(?!\w)", re.IGNORECASE)

        def repl(match):
            found = match.group(0)
            found_genes.append(found.upper())  # or normalized casing
            return f"**{found}**"

        processed_text = pattern.sub(repl, processed_text)
        
        return processed_text.strip(), found_genes

    def search_literature(self,
                          genes: List[str],
                          email: str,
                          search_terms: List[str]) -> List[Dict[str, Any]]:
        """Search for scientific articles related to the gene list using PubMed."""

        Entrez.email = email
        Entrez.tool = 'gene_enrichment_agent'
        Entrez.api_key = self.api_key

        articles = []
        terms = 'AND (' + '[MeSH Terms] OR '.join(search_terms) + '[MeSH Terms]) ' if search_terms else ''
        search_query = f"({'[tw] OR '.join(genes)}[tw]) {terms}AND 2015:3000[PDAT]"

        # Get PubMed IDs
        search_handle = Entrez.esearch(
            db="pubmed",
            term=search_query,
            retmax=self.num_papers,
            sort="relevance",
            timeout=30
        )
        done = False
        while not done:
            try:
                search_results = Entrez.read(search_handle)
                search_handle.close()
                done = True
            except IncompleteRead as e:
                print(f"Incomplete read: {str(e)}")
                time.sleep(1)
                continue

        pubmed_ids = search_results['IdList']
        if not pubmed_ids:
            return []

        # Get PubMed articles
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=pubmed_ids,
            rettype="xml",
            retmode="xml",
            timeout=30
        )

        done = False
        while not done:
            try:
                records = Entrez.read(fetch_handle)
                fetch_handle.close()
                done = True
            except IncompleteRead as e:
                print(f"Incomplete read: {str(e)}")
                time.sleep(1)
                continue

        # Batch elink to PMC to find which articles have full text
        link_handle = Entrez.elink(
            dbfrom="pubmed",
            db="pmc",
            id=pubmed_ids,
            retmode="json",
            timeout=30
        )
        link_results = json.loads(link_handle.read().decode('utf-8'))

        # Build a mapping: PubMed ID -> PMC ID
        pubmed_to_pmc = {}
        for result in link_results['linksets']:
            pmid = result["ids"][0]
            linkset = result.get('linksetdbs', [])
            if linkset and 'links' in linkset[0]:
                pmcid = linkset[0]['links'][0]
                pubmed_to_pmc[pmid] = pmcid

        # Process each article
        for record in records['PubmedArticle']:
            try:
                article_data = record['MedlineCitation']
                pmid = article_data['PMID']

                # Process title
                title = article_data['Article']['ArticleTitle']
                processed_title, title_genes = self._highlight_genes_in_text(title, genes)

                # Process abstract
                abstract = ''.join(article_data['Article'].get('Abstract', {}).get('AbstractText', []))
                processed_abstract, abstract_genes = self._highlight_genes_in_text(abstract, genes)

                article = {
                    'name': processed_title,
                    'year': article_data['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                    'abstract': processed_abstract,
                    'id': pmid,
                    'genes': list(set(title_genes + abstract_genes))
                }

                # If there's a linked PMC ID, fetch full text
                if pmid in pubmed_to_pmc:
                    pmcid = pubmed_to_pmc[pmid]
                    try:
                        time.sleep(.1)
                        pmc_handle = Entrez.efetch(
                            db="pmc",
                            id=pmcid,
                            rettype="xml",
                            retmode="xml",
                            timeout=30
                        )
                        pmc_xml = pmc_handle.read().decode('utf-8')
                        pmc_handle.close()

                        root = ET.fromstring(pmc_xml)
                        gene_mentions = []

                        for p in root.findall('.//p'):
                            text = self._extract_text_from_element(p)
                            processed_text, found_genes = self._highlight_genes_in_text(text, genes)
                            if found_genes:
                                clean_text = ' '.join(processed_text.split())
                                gene_mentions.append(clean_text)
                                for gene in found_genes:
                                    if gene not in article['genes']:
                                        article['genes'].append(gene)

                        if gene_mentions:
                            article['gene_mentions'] = gene_mentions[:3]

                    except Exception as e:
                        print(f"Error fetching PMC content for PMID {pmid}: {str(e)}")

                articles.append(article)
            except Exception as e:
                print(f"Error processing article {pmid}: {str(e)}")
                continue

        return articles
    
    def fetch_gene_summaries(self, genes: List[str], organism: str = "Homo sapiens") -> List[Dict[str, Any]]:
        """Fetch gene summaries from the Entrez database using batch operations to minimize API calls."""
        if not genes:
            return []
            
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        
        # Build a single search query for all genes
        gene_terms = ' OR '.join([f'"{gene}"[Gene Name]' for gene in genes])
        search_query = f"({gene_terms}) AND {organism}[Organism]"
        
        # Single search call for all genes
        search_params = {
            "db": "gene",
            "term": search_query,
            "retmax": len(genes) * 2,  # Allow for potential duplicates
            "retmode": "json",
            "api_key": self.api_key
        }
        
        try:
            search_resp = requests.get(search_url, params=search_params)
            search_resp.raise_for_status()
            search_data = search_resp.json()
            
            gene_ids = search_data["esearchresult"]["idlist"]
            
            if not gene_ids:
                return []
            
            # Single summary call for all gene IDs
            summary_params = {
                "db": "gene",
                "id": ",".join(gene_ids),
                "retmode": "json",
                "api_key": self.api_key
            }
            
            summary_resp = requests.get(summary_url, params=summary_params)
            summary_resp.raise_for_status()
            summary_data = summary_resp.json()
            
            # Process results and match genes to their summaries
            gene_summaries = {}
            
            # First, build a mapping of gene names to IDs from the search results
            for gene_id in gene_ids:
                if gene_id in summary_data["result"]:
                    gene_info = summary_data["result"][gene_id]
                    gene_name = gene_info.get("name", "")
                    if gene_name:
                        if gene_name.upper() in gene_summaries:
                            continue
                        summary = {
                            "name": gene_name,
                            "id": int(gene_id),
                            "description": gene_info.get("description", ""),
                            "summary": gene_info.get("summary", "")
                        }
                        gene_summaries[gene_name.upper()] = summary
            
            return list(gene_summaries.values())
            
        except Exception as e:
            print(f"Error fetching gene summaries: {str(e)}")
            return []