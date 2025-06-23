"""Summarization functionality for gene enrichment analysis."""

import json
from typing import Dict, List
from openai import OpenAI
from pydantic import BaseModel
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

class SummarizeAnalyzer:
    """Class for summarizing gene enrichment analysis results."""

    def __init__(self, api_key: str, model: str = "o3"):
        """Initialize the summarize analyzer.
        Args:
            api_key: OpenAI API key
            model: OpenAI model to use
        """
        self.client = OpenAI(api_key=api_key)
        self.summarize_model = model
    
    def _get_full_term(self,  name: str, source: str, results: Dict) -> Dict:
        """Get the full term from the results."""
        try:
            for term in results['results'][source]:
                if name.lower() in term['name'].lower():
                    return term
        except KeyError:
            print(f"Term {name} not found in {source} results")
            return None

    def group_results_by_theme(self,
                               enrichr_results: Dict,
                               toppfun_results: Dict,
                               gprofiler_results: Dict,
                               literature_results: Dict,
                               genes: List[str],
                               ranked: bool,
                               context: str) -> str:
        """Group enrichment results into functional themes across all tools.
        
        Args:
            enrichr_results: Results from Enrichr analysis
            toppfun_results: Results from ToppFun analysis
            gprofiler_results: Results from gProfiler analysis
            genes: List of genes to analyze
            ranked: If the genes are ranked
            context: Additional user context
            
        Returns:
            Themed summary of results with each theme containing related terms
        """

        # Generate barcoded shortlist
        barcode_dict = dict()
        combined_shortlist = []
        for tool in [enrichr_results, toppfun_results, gprofiler_results]:
            combined_shortlist.extend(tool.get('shortlist', []))
        combined_shortlist.extend(literature_results['results']['pubmed'])

        barcode = 100000
        for term in combined_shortlist:
            term['barcode'] = barcode
            barcode_dict[barcode] = (term['name'], term['source'], term.pop('tool'))
            barcode += 1
        self.barcode_dict = barcode_dict

        # Pydantic model for output
        class LLMTheme(BaseModel):
            theme: str
            description: str
            barcodes: List[int]

        class LLMThemedResults(BaseModel):
            themes: List[LLMTheme]
            summary: str

        # Send to LLM
        prompt = f"""You will be given a list of genes that characterize a group of cells{', ranked by differential expression.' if ranked else '.'}
        Your goal is to determine the biological processes and pathways that are enriched in these genes.
        You will also be given enrichment results from several different tools to help you with this task: Enrichr, ToppFun, gProfiler, and PubMed.
        
        Your task is to analyze these enrichment results and arrange them into functional themes.
        {'You should focus on themes that involve genes towards the top of the differential expression list.' if ranked else ''}
        Feel free to delete any terms that don't fit a theme.
        You should include literature terms in themes as they fit, but your final theme should just be a Literature Findings theme.

        You will return a list of themes with the following attributes:
        * theme: The name of the theme
        * description: A brief description of the function of the theme and why you identified it
        * barcodes: A list of barcodes, unique identifiers for the terms that are associated with the theme

        You will also provide a summary of the results, including a high level overview of what these cells are enriched for.

        Guidelines:
        * Focus on biological meaning rather than technical categories
        * Prioritize themes with strong support across multiple tools
        {'* Focus on themes that involve genes towards the top of the differential expression list' if ranked else ''}
        """
        
        response = self.client.responses.parse(
            model=self.summarize_model,
            input=[
                {"role": "system", "content": "You are an expert in bioinformatics, immunology, and oncology specializing in gene enrichment analysis."},
                {"role": "user", "content": prompt},
                {"role": "user", "content": 'Context: ' + context},
                {"role": "user", "content": 'Genes: ' + ';'.join(genes)},
                {"role": "user", "content": 'Enrichment results:\n' + json.dumps(combined_shortlist, indent=1)}
            ],
            text_format=LLMThemedResults
        )

        themed_results = response.output_parsed

        # Get full terms from barcodes
        clean_themed_results = {
            'summary': themed_results.summary,
            'themes': []
        }
        result_dict = {
            'enrichr': enrichr_results,
            'toppfun': toppfun_results,
            'gprofiler': gprofiler_results,
            'literature': literature_results
        }
        for theme in themed_results.themes:
            temp_theme = {
                'theme': theme.theme,
                'description': theme.description,
                'terms': []
            }
            for barcode in theme.barcodes:
                if barcode in barcode_dict:
                    name, source, tool = barcode_dict[barcode]
                    term = self._get_full_term(name, source, result_dict[tool])
                    term['tool'] = tool
                    term['source'] = source
                    temp_theme['terms'].append(term)
            clean_themed_results['themes'].append(temp_theme)

        return clean_themed_results

    def _sanitize_sheet_name(self, name: str, max_length: int = 31) -> str:
        """Sanitize a string to be used as an Excel sheet name.
        
        Excel sheet names cannot contain: : / \ ? * [ ]
        They also cannot start with a number or be longer than 31 characters.
        
        Args:
            name: The original name to sanitize
            max_length: Maximum length for the sheet name (default 31)
            
        Returns:
            Sanitized sheet name
        """
        # Replace invalid characters with spaces or underscores
        invalid_chars = [':', '/', '\\', '?', '*', '[', ']']
        sanitized = name
        for char in invalid_chars:
            sanitized = sanitized.replace(char, '_')
        
        # Remove leading/trailing whitespace
        sanitized = sanitized.strip()
        
        # If the name starts with a number, add a prefix
        if sanitized and sanitized[0].isdigit():
            sanitized = 'Theme_' + sanitized
        
        # Truncate to max length
        if len(sanitized) > max_length:
            sanitized = sanitized[:max_length]
        
        # Ensure it's not empty
        if not sanitized:
            sanitized = 'Theme'
            
        return sanitized

    def _extract_term_id(self, term: Dict, tool: str) -> str:
        """Extract term ID (GO:XXXXXX or PMID) from a term based on the tool.
        
        Args:
            term: Term dictionary containing the ID information
            tool: Name of the tool that generated this term
            
        Returns:
            Term ID if found, empty string otherwise
        """
        if tool == 'literature':
            # Extract PMID from pmid attribute
            if 'pmid' in term:
                return f"PMID:{term['pmid']}"
        elif tool == 'enrichr':
            # Extract GO:XXXXXX from name which has format "Term Name (GO:XXXXXX)"
            if '(' in term['name'] and ')' in term['name']:
                go_term = term['name'].split('(')[-1].strip(')')
                if go_term.startswith('GO:'):
                    return go_term
        elif tool == 'gprofiler':
            # GO term is stored in native attribute
            if 'native' in term and term['native'].startswith('GO:'):
                return term['native']
        elif tool == 'toppfun':
            # GO term is stored in id attribute
            if 'id' in term and term['id'].startswith('GO:'):
                return term['id']
        return ''

    def synthesize_analysis(self,
                            themed_results: Dict,
                            genes: List[str],
                            email: str,
                            search_terms: List[str],
                            context: str,
                            results_dir: str) -> str:
        """Synthesize themed results into an Excel file with multiple sheets.
        
        Args:
            themed_results: Dict with the following keys:
                themes: List of themes, each containing:
                        - theme: Name of the theme
                        - description: Description of the theme
                        - terms: List of terms with tool-specific information
                summary: String with description of enrichment results
            genes: List of genes in analysis
            email: Email provided
            search_terms: List of search terms provided
            context: Context provided
            results_dir: String with path to results
                
        Returns:
            Path to the generated Excel file
        """
        # Create a new workbook
        wb = Workbook()
        
        # Create summary sheet
        summary_sheet = wb.active
        summary_sheet.title = "Summary"
        
        # Add metadata section
        summary_sheet['A1'] = "Analysis Metadata"
        summary_sheet['A1'].font = Font(bold=True, size=14)
        
        summary_sheet['A3'] = "Email:"
        summary_sheet['B3'] = email
        summary_sheet['A3'].font = Font(bold=True)
        
        summary_sheet['A4'] = "Search Terms:"
        summary_sheet['B4'] = ", ".join(search_terms)
        summary_sheet['A4'].font = Font(bold=True)
        
        summary_sheet['A6'] = "Context:"
        summary_sheet['B6'] = context
        summary_sheet['A6'].font = Font(bold=True)
        summary_sheet['B6'].alignment = Alignment(wrap_text=True)
        
        summary_sheet['A8'] = "Genes:"
        summary_sheet['B8'] = ", ".join(genes)
        summary_sheet['A8'].font = Font(bold=True)
        summary_sheet['B8'].alignment = Alignment(wrap_text=True)
        
        # Add summary section
        summary_sheet['A10'] = "Analysis Summary"
        summary_sheet['A10'].font = Font(bold=True, size=14)
        summary_sheet['A11'] = "Description:"
        summary_sheet['A11'].font = Font(bold=True)
        summary_sheet['B11'] = themed_results['summary']
        summary_sheet['B11'].alignment = Alignment(wrap_text=True)
        summary_sheet['A13'] = "Themes:"
        summary_sheet['A13'].font = Font(bold=True)
        for i, theme in enumerate(themed_results['themes']):
            summary_sheet['B'+str(13+i)] = theme['theme']
        
        # Adjust column widths
        summary_sheet.column_dimensions['A'].width = 20
        summary_sheet.column_dimensions['B'].width = 100
        
        # Collect all tools first
        tool_set = set()
        for theme in themed_results['themes']:
            for term in theme['terms']:
                if term['tool'] != 'literature':
                    tool_set.add(term['tool'])
        sorted_tools = sorted(tool_set)
        
        # Create header row
        header = ["Name", "Source", "ID", "Genes"] + [tool for tool in sorted_tools]
        
        # Process each theme
        for theme in themed_results['themes'][:-1]:
            # Create new sheet for theme
            theme_sheet = wb.create_sheet(title=self._sanitize_sheet_name(theme['theme']))
            
            # Add theme name
            theme_sheet['A1'] = theme['theme']
            theme_sheet['A1'].font = Font(bold=True)
            
            # Add theme description
            theme_sheet['A2'] = theme['description']
            theme_sheet['A2'].alignment = Alignment(wrap_text=True)
            
            # Add Terms header
            theme_sheet['A4'] = "Terms:"
            theme_sheet['A4'].font = Font(bold=True)
            
            # Add column headers
            for col, header_text in enumerate(header, 1):
                cell = theme_sheet.cell(row=5, column=col)
                cell.value = header_text
                cell.font = Font(bold=True)
            
            # Group terms by ID
            id_terms = {}  # key: ID -> value: {tool: (term_name, source, pvalue, genes)}
            non_id_terms = []  # list of terms without IDs
            
            for term in theme['terms']:
                term_id = self._extract_term_id(term, term['tool'])
                term_name = term['name'].replace(',', ';')
                source = term['source']
                genes_list = term.get('genes', [])
                genes_str = " ".join(sorted(genes_list)) if genes_list else ""
                
                # Extract p-value
                pvalue = term.get('p_value', term.get('pvalue', term.get('adjPvalue', '')))
                if isinstance(pvalue, float):
                    pvalue = f"{pvalue:.2e}"
                
                if term_id:
                    if term_id not in id_terms:
                        id_terms[term_id] = {tool: ('', '', '', '') for tool in sorted_tools}
                    id_terms[term_id][term['tool']] = (term_name, source, pvalue, genes_str)
                else:
                    non_id_terms.append((term_name, source, genes_str, term['tool'], pvalue))
            
            # Add ID term rows first
            current_row = 6
            for term_id, tool_values in id_terms.items():
                # Find the most informative term name and source
                term_names = [v[0] for v in tool_values.values() if v[0]]
                sources = [v[1] for v in tool_values.values() if v[1]]
                all_genes = set()
                for v in tool_values.values():
                    if v[3]:  # genes string
                        all_genes.update(v[3].split(';'))
                
                if not term_names or not sources:
                    continue
                    
                # Use the longest term name as it's likely most descriptive
                term_name = max(term_names, key=len)
                # Prefer GO source if available, otherwise take the first one
                source = next((s for s in sources if 'GO' in s.upper()), sources[0])
                genes_str = " ".join(sorted(all_genes)) if all_genes else ""
                
                # Add row data
                theme_sheet.cell(row=current_row, column=1, value=term_name)
                theme_sheet.cell(row=current_row, column=2, value=source)
                theme_sheet.cell(row=current_row, column=3, value=term_id)
                theme_sheet.cell(row=current_row, column=4, value=genes_str)
                
                # Add p-values for each tool
                for col, tool in enumerate(sorted_tools, 5):
                    p_value = float(tool_values[tool][2]) if tool_values[tool][2] else ''
                    theme_sheet.cell(row=current_row, column=col, value=p_value)
                
                current_row += 1
            
            # Add non-ID term rows
            for term_name, source, genes_str, tool, pvalue in non_id_terms:
                theme_sheet.cell(row=current_row, column=1, value=term_name)
                theme_sheet.cell(row=current_row, column=2, value=source)
                theme_sheet.cell(row=current_row, column=3, value="")
                theme_sheet.cell(row=current_row, column=4, value=genes_str)
                
                # Add p-values for each tool
                for col, t in enumerate(sorted_tools, 5):
                    theme_sheet.cell(row=current_row, column=col, value=float(pvalue) if t == tool else '')
                
                current_row += 1
            
            # Adjust column widths
            theme_sheet.column_dimensions['A'].width = 80
            for col in range(2, len(header) + 1):
                theme_sheet.column_dimensions[chr(64 + col)].width = 20
        
        # Literature Findings theme
        theme = themed_results['themes'][-1]
        theme_sheet = wb.create_sheet(title=self._sanitize_sheet_name(theme['theme']))
            
        # Add theme name
        theme_sheet['A1'] = theme['theme']
        theme_sheet['A1'].font = Font(bold=True)
        
        # Add theme description
        theme_sheet['A2'] = theme['description']
        theme_sheet['A2'].alignment = Alignment(wrap_text=True)
        # Add headers for literature information
        headers = ['Paper Title', 'Year', 'PMID', 'Abstract', 'Genes', 'Gene Mentions']
        for col, header in enumerate(headers, 1):
            cell = theme_sheet.cell(row=4, column=col)
            cell.value = header
            cell.font = Font(bold=True)
        
        # Add paper information
        current_row = 5
        for term in theme['terms']:
            if term['tool'] == 'literature':
                # Paper title
                theme_sheet.cell(row=current_row, column=1, value=term['name'])
                theme_sheet.cell(row=current_row, column=1).alignment = Alignment(wrap_text=True)
                
                # Year
                theme_sheet.cell(row=current_row, column=2, value=int(term.get('year')))
                
                # PMID
                theme_sheet.cell(row=current_row, column=3, value=int(term.get('pmid')))
                
                # Abstract
                theme_sheet.cell(row=current_row, column=4, value=term.get('abstract', ''))
                theme_sheet.cell(row=current_row, column=4).alignment = Alignment(wrap_text=True)
                
                # Genes
                genes = term.get('genes', [])
                theme_sheet.cell(row=current_row, column=5, value=', '.join(genes))
                
                # Gene mentions
                gene_mentions = term.get('gene_mentions', [])
                for quote in gene_mentions:
                    theme_sheet.cell(row=current_row, column=6, value=quote)
                    theme_sheet.cell(row=current_row, column=6).alignment = Alignment(wrap_text=True)
                    current_row += 1
                
                current_row += 1
        
        # Adjust column widths
        theme_sheet.column_dimensions['A'].width = 80  # Title
        theme_sheet.column_dimensions['B'].width = 10  # Year
        theme_sheet.column_dimensions['C'].width = 15  # PMID
        theme_sheet.column_dimensions['D'].width = 60  # Abstract
        theme_sheet.column_dimensions['E'].width = 20  # Genes
        theme_sheet.column_dimensions['F'].width = 60  # Gene mentions

        
        # Save the workbook
        output_path = f"{results_dir}/enrichment_analysis.xlsx"
        wb.save(output_path)
        return output_path 