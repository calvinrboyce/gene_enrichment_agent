import os
import csv
import dotenv
from gene_enrichment_agent import GeneEnrichmentAgent

def get_user_input():
    """Get user input for email, search terms, and context."""
    print("=== Gene Enrichment Analysis Setup ===")
    
    # Get email
    email = input("Please enter your email address (required for NCBI literature search): \n").strip()
    while not email or '@' not in email:
        print("Please enter a valid email address.")
        email = input("Email: \n").strip()
    print()
    
    # Get search terms
    search_terms_input = input("Enter search terms for literature analysis (comma-separated, or press Enter for none): \n").strip()
    search_terms = [term.strip() for term in search_terms_input.split(',')] if search_terms_input else []
    print()
    
    # Get context
    context = input("Enter context for the analysis (description of your study): \n").strip()
    if not context:
        context = "Gene enrichment analysis"
    print()
    
    return email, search_terms, context

def read_gene_lists_from_csv(csv_file):
    """Read gene lists from CSV file where first row contains analysis names."""
    gene_lists = {}
    
    try:
        with open(csv_file, 'r', encoding='utf-8-sig') as file:  # utf-8-sig handles BOM
            reader = csv.reader(file)
            headers = next(reader)  # First row contains analysis names
            
            # Initialize lists for each column
            for header in headers:
                gene_lists[header] = []
            
            # Read gene data
            for row in reader:
                for i, gene in enumerate(row):
                    if i < len(headers) and gene.strip():  # Only add non-empty genes
                        gene_lists[headers[i]].append(gene.strip())
        
        return gene_lists
    
    except FileNotFoundError:
        print(f"Error: Could not find {csv_file}")
        return {}
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return {}

def read_background_genes(txt_file):
    """Read background genes from a text file."""
    with open(txt_file, 'r') as file:
        return file.read().splitlines()

def main():
    # Load environment variables
    dotenv.load_dotenv()
    openai_api_key = os.getenv("OPENAI_API_KEY")
    
    if not openai_api_key:
        print("Error: OPENAI_API_KEY not found in environment variables.")
        print("Please set your OpenAI API key in a .env file or environment variable.")
        return
    
    # Initialize the agent
    print("Initializing Gene Enrichment Agent...")
    agent = GeneEnrichmentAgent(openai_api_key, papers_per_gene=0)
    
    # Get user input
    email, search_terms, context = get_user_input()

    # Read background genes from text file
    background_genes = read_background_genes("background_genes.txt")
    
    # Read gene lists from CSV
    csv_file = "enrichment_data.csv"
    gene_lists = read_gene_lists_from_csv(csv_file)
    
    if not gene_lists:
        print("No gene lists found. Exiting.")
        return
    
    # Display available analyses
    print(f"\nFound {len(gene_lists)} gene lists to analyze:")
    for analysis_name, genes in gene_lists.items():
        print(f"  - {analysis_name}: {len(genes)} genes")
    
    # Confirm before proceeding
    print(f"\nProceeding with analysis using:")
    print(f"  Email: {email}")
    print(f"  Search terms: {search_terms if search_terms else 'None'}")
    print(f"  Context: {context}")
    
    proceed = input("\nProceed with analysis? (y/n): ").strip().lower()
    if proceed != 'y':
        print("Analysis cancelled.")
        return
    
    # Run analysis for each gene list
    print("\n=== Starting Analysis ===")
    all_results = {}
    
    for analysis_name, genes in gene_lists.items():
        if not genes:  # Skip empty gene lists
            continue
            
        print(f"\n--- Analyzing {analysis_name} ({len(genes)} genes) ---")
        
        try:
            results = agent.run_analysis(
                genes=genes,
                email=email,
                background_genes=background_genes,
                ranked=True,
                search_terms=search_terms,
                context=context,
                analysis_name=analysis_name
            )
            
            all_results[analysis_name] = results
            
            # Display summary
            print(f"\nResults for {analysis_name}:")
            print(f"Summary: {results['summary']}")
            
        except Exception as e:
            print(f"Error analyzing {analysis_name}: {e}")
            continue
    
    print(f"\n=== Analysis Complete ===")
    print(f"Successfully analyzed {len(all_results)} gene lists.")
    print("Check the 'gene_enrichment_analysis_results' directory for detailed results.")

if __name__ == "__main__":
    main()