"""Literature search functionality using PubMed."""

import time
import xml.etree.ElementTree as ET
from typing import List, Dict, Any
from Bio import Entrez
import re

class LiteratureAnalyzer:
    """Class for searching scientific literature related to genes."""

    def __init__(self, papers_per_gene: int = 2, max_papers: int = 10):
        """Initialize the literature analyzer."""
        self.papers_per_gene = papers_per_gene
        self.max_papers = max_papers

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
        for gene in genes:
            # Use regex for case-insensitive, whole-word matching
            pattern = re.compile(rf"(?<!\w){re.escape(gene)}(?!\w)", re.IGNORECASE)
            if pattern.search(processed_text):
                # Replace all case-insensitive matches with the original gene name in bold
                def repl(match):
                    matched_gene = match.group(0)
                    # Track the gene as found if not already (case-insensitive)
                    if gene not in found_genes:
                        found_genes.append(gene)
                    return f"**{matched_gene}**"
                processed_text = pattern.sub(repl, processed_text)
        
        return processed_text.strip(), found_genes

    def search_literature(self,
                          queries: List[str],
                          email: str,
                          search_terms: List[str]) -> List[Dict[str, Any]]:
        """Search for scientific articles related to the gene list using PubMed.
        
        Args:
            queries: List of queries to search for
            email: Email address for NCBI's reference
            search_terms: List of search terms to use in the literature search
            papers_per_gene: Number of papers to retrieve per gene
            max_papers: Maximum number of papers to retrieve for full list
        Returns:
            List of relevant articles with metadata including PMC full text excerpts if available
        """
        # Set your email for NCBI's reference
        Entrez.email = email
        Entrez.tool = 'gene_enrichment_agent'
        
        articles = []
        for query in queries:
            try:
                genes_to_search = query if isinstance(query, list) else [query]
                terms = 'AND (' + '[MeSH Terms] OR '.join(search_terms) + '[MeSH Terms]) ' if len(search_terms)>0 else ''
                search_query = f"({'[tw] OR '.join(genes_to_search)}[tw]) {terms}AND 2015:3000[PDAT]"
                retmax = self.max_papers if isinstance(query, list) else self.papers_per_gene
                
                # Search PubMed
                search_handle = Entrez.esearch(
                    db="pubmed",
                    term=search_query,
                    retmax=retmax,
                    sort="relevance"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if not search_results['IdList']:
                    continue
                    
                # Fetch detailed article information
                fetch_handle = Entrez.efetch(
                    db="pubmed",
                    id=search_results['IdList'],
                    rettype="xml",
                    retmode="xml"
                )
                records = Entrez.read(fetch_handle)
                fetch_handle.close()
                
                for record in records['PubmedArticle']:
                    article_data = record['MedlineCitation']
                    
                    # Get and process title
                    title = article_data['Article']['ArticleTitle']
                    processed_title, title_genes = self._highlight_genes_in_text(title, genes_to_search)
                    
                    # Get and process abstract
                    abstract = ''.join(article_data['Article'].get('Abstract', {}).get('AbstractText', []))
                    processed_abstract, abstract_genes = self._highlight_genes_in_text(abstract, genes_to_search)
                    
                    # Initialize article with processed text and found genes
                    article = {
                        'name': processed_title,
                        'source': 'pubmed',
                        'tool': 'literature',
                        'year': article_data['Article']['Journal']['JournalIssue']['PubDate'].get('Year', ''),
                        'abstract': processed_abstract,
                        'pmid': article_data['PMID'],
                        'genes': list(set(title_genes + abstract_genes))  # Initialize with genes from title/abstract
                    }
                    
                    # Check if article is available in PMC
                    link_handle = Entrez.elink(
                        dbfrom="pubmed",
                        db="pmc",
                        id=article['pmid']
                    )
                    link_results = Entrez.read(link_handle)
                    link_handle.close()
                    
                    # If available in PMC, try to get relevant excerpt
                    if (link_results and link_results[0].get('LinkSetDb') and 
                        link_results[0]['LinkSetDb'] and link_results[0]['LinkSetDb'][0].get('Link')):
                        pmcid = link_results[0]['LinkSetDb'][0]['Link'][0]['Id']
                        try:
                            pmc_handle = Entrez.efetch(
                                db="pmc",
                                id=pmcid,
                                rettype="xml",
                                retmode="xml"
                            )
                            pmc_xml = pmc_handle.read().decode('utf-8')
                            pmc_handle.close()
                            
                            # Parse XML
                            root = ET.fromstring(pmc_xml)
                            
                            # Find paragraphs mentioning the gene
                            gene_mentions = []
                            # Look for paragraphs in the body
                            for p in root.findall('.//p'):
                                text = self._extract_text_from_element(p)
                                processed_text, found_genes = self._highlight_genes_in_text(text, genes_to_search)
                                if found_genes:
                                    clean_text = ' '.join(processed_text.split())
                                    if clean_text:
                                        gene_mentions.append(clean_text)
                                    # Add any newly found genes to the article's gene list
                                    for gene in found_genes:
                                        if gene not in article['genes']:
                                            article['genes'].append(gene)
                            
                            if gene_mentions:
                                article['gene_mentions'] = gene_mentions[:3]  # Limit to 3 mentions
                        except Exception as e:
                            print(f"Error fetching PMC content: {str(e)}")
                    
                    articles.append(article)
                    time.sleep(0.5)  # Be nice to NCBI servers
                    
            except Exception as e:
                print(f"Error searching literature for query {query}: {str(e)}")
                continue
                
        return {'results': {'pubmed': articles}} 