#!/usr/bin/env python3
"""
Primer Validation Script for Mouse Genes
Checks if primers target the correct genes and generate 70-150bp amplicons

Usage:
    python check_primers.py                          # Interactive file selection
    python check_primers.py input_file.xlsx          # Specify file directly
    python check_primers.py input_file.xlsx -o output.xlsx  # Custom output name
"""

import pandas as pd
import requests
import time
import re
import sys
import os
import argparse

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def fetch_gene_sequence(gene_name, organism="Mus musculus"):
    """
    Fetch gene sequence from NCBI for a given gene name and organism
    Returns genomic_sequence, cds_sequence, header
    """
    print(f"  Fetching sequences for {gene_name}...", end=" ")
    
    # Search for the gene
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    genomic_seq = None
    cds_seq = None
    header = None
    
    try:
        # First, try to get mRNA/CDS sequence (this is what primers are usually designed for)
        mrna_search_params = {
            "db": "nucleotide",
            "term": f"{gene_name}[Gene] AND {organism}[Organism] AND refseq[filter] AND mrna[filter]",
            "retmode": "json",
            "retmax": 5
        }
        
        mrna_search_response = requests.get(search_url, params=mrna_search_params, timeout=10)
        mrna_search_data = mrna_search_response.json()
        
        if mrna_search_data.get("esearchresult", {}).get("idlist"):
            mrna_id = mrna_search_data["esearchresult"]["idlist"][0]
            
            # Fetch mRNA sequence (this is the spliced transcript without introns)
            mrna_fetch_params = {
                "db": "nucleotide",
                "id": mrna_id,
                "rettype": "fasta",
                "retmode": "text"
            }
            
            mrna_fetch_response = requests.get(fetch_url, params=mrna_fetch_params, timeout=10)
            mrna_fasta = mrna_fetch_response.text
            
            if mrna_fasta and '>' in mrna_fasta:
                mrna_lines = mrna_fasta.strip().split('\n')
                header = mrna_lines[0]
                cds_seq = ''.join(mrna_lines[1:]).upper()
        
        # Now try to get genomic DNA sequence (with introns)
        genomic_search_params = {
            "db": "nucleotide",
            "term": f"{gene_name}[Gene] AND {organism}[Organism] AND refseq[filter] AND genomic[filter]",
            "retmode": "json",
            "retmax": 5
        }
        
        genomic_search_response = requests.get(search_url, params=genomic_search_params, timeout=10)
        genomic_search_data = genomic_search_response.json()
        
        if genomic_search_data.get("esearchresult", {}).get("idlist"):
            genomic_id = genomic_search_data["esearchresult"]["idlist"][0]
            
            # Fetch genomic sequence
            genomic_fetch_params = {
                "db": "nucleotide",
                "id": genomic_id,
                "rettype": "fasta",
                "retmode": "text"
            }
            
            genomic_fetch_response = requests.get(fetch_url, params=genomic_fetch_params, timeout=10)
            genomic_fasta = genomic_fetch_response.text
            
            if genomic_fasta and '>' in genomic_fasta:
                genomic_lines = genomic_fasta.strip().split('\n')
                if not header:  # Use genomic header if we didn't get mRNA
                    header = genomic_lines[0]
                genomic_seq = ''.join(genomic_lines[1:]).upper()
        
        # Print results
        if not genomic_seq and not cds_seq:
            print("❌ No sequences found")
            return None, None, None
        
        status_parts = []
        if cds_seq:
            status_parts.append(f"mRNA: {len(cds_seq)} bp")
        if genomic_seq:
            status_parts.append(f"Genomic: {len(genomic_seq)} bp")
        
        print(f"✓ ({', '.join(status_parts)})")
        return genomic_seq, cds_seq, header
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return None, None, None

def find_primer_binding(primer_seq, gene_seq, primer_name="primer"):
    """
    Find where a primer binds in the gene sequence
    Returns list of binding positions
    """
    positions = []
    primer_upper = primer_seq.upper()
    
    # Search for exact matches
    start = 0
    while True:
        pos = gene_seq.find(primer_upper, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    
    return positions

def check_primer_pair(fwd_seq, rev_seq, genomic_seq, cds_seq, gene_name):
    """
    Check if primer pair works and calculate amplicon size
    Checks mRNA first (since primers are usually designed for mRNA), then genomic
    Returns dict with results
    """
    result = {
        'gene': gene_name,
        'fwd_primer': fwd_seq,
        'rev_primer': rev_seq,
        'fwd_binds': False,
        'rev_binds': False,
        'amplicon_size': None,
        'status': 'FAIL',
        'details': '',
        'sequence_type': None
    }
    
    if genomic_seq is None and cds_seq is None:
        result['details'] = 'Gene sequence not found in NCBI'
        return result
    
    # Try mRNA/CDS first (this is what most primers are designed for)
    if cds_seq:
        cds_result = check_single_sequence(fwd_seq, rev_seq, cds_seq, gene_name)
        if cds_result['status'] in ['PASS', 'WARN']:
            cds_result['sequence_type'] = 'mRNA'
            return cds_result
    
    # If mRNA didn't work, try genomic DNA
    if genomic_seq:
        genomic_result = check_single_sequence(fwd_seq, rev_seq, genomic_seq, gene_name)
        if genomic_result['status'] in ['PASS', 'WARN']:
            genomic_result['sequence_type'] = 'Genomic'
            return genomic_result
        
        # If both failed, report genomic failure with note that mRNA was also tried
        if cds_seq:
            genomic_result['sequence_type'] = 'Genomic (mRNA also checked)'
            return genomic_result
        else:
            genomic_result['sequence_type'] = 'Genomic'
            return genomic_result
    
    # Only mRNA available and it failed
    if cds_seq:
        cds_result['sequence_type'] = 'mRNA'
        return cds_result
    
    return result

def check_single_sequence(fwd_seq, rev_seq, gene_seq, gene_name):
    """
    Check if primer pair works on a single sequence
    Returns dict with results
    """
    result = {
        'gene': gene_name,
        'fwd_primer': fwd_seq,
        'rev_primer': rev_seq,
        'fwd_binds': False,
        'rev_binds': False,
        'amplicon_size': None,
        'status': 'FAIL',
        'details': ''
    }
    
    if gene_seq is None:
        result['details'] = 'Sequence not available'
        return result
    
    # Find forward primer binding sites
    fwd_positions = find_primer_binding(fwd_seq, gene_seq)
    
    # Find reverse primer binding sites (need to check reverse complement)
    rev_comp = reverse_complement(rev_seq)
    rev_positions = find_primer_binding(rev_comp, gene_seq)
    
    result['fwd_binds'] = len(fwd_positions) > 0
    result['rev_binds'] = len(rev_positions) > 0
    
    if not result['fwd_binds']:
        result['details'] = 'Forward primer does not bind to gene'
        return result
    
    if not result['rev_binds']:
        result['details'] = 'Reverse primer does not bind to gene'
        return result
    
    # Calculate all possible amplicon sizes
    amplicon_sizes = []
    for fwd_pos in fwd_positions:
        for rev_pos in rev_positions:
            if rev_pos > fwd_pos:  # Reverse primer must be downstream
                # Amplicon size is from start of forward primer to end of reverse primer
                size = (rev_pos + len(rev_comp)) - fwd_pos
                amplicon_sizes.append(size)
    
    if not amplicon_sizes:
        result['details'] = 'Primers bind but in wrong orientation (no amplicon)'
        return result
    
    # Use the smallest amplicon size (most likely product)
    result['amplicon_size'] = min(amplicon_sizes)
    
    # Check if amplicon is in the target range (70-150 bp)
    if 70 <= result['amplicon_size'] <= 150:
        result['status'] = 'PASS'
        result['details'] = f'Amplicon size: {result["amplicon_size"]} bp ✓'
    else:
        result['status'] = 'WARN'
        if result['amplicon_size'] < 70:
            result['details'] = f'Amplicon too short: {result["amplicon_size"]} bp (need 70-150 bp)'
        else:
            result['details'] = f'Amplicon too long: {result["amplicon_size"]} bp (need 70-150 bp)'
    
    return result

def calculate_tm(seq):
    """
    Calculate melting temperature using Wallace rule for short primers
    and more accurate formula for longer primers
    """
    seq = seq.upper()
    length = len(seq)
    
    if length < 14:
        # Wallace rule: Tm = 2(A+T) + 4(G+C)
        tm = 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))
    else:
        # More accurate formula for primers 14+ bases
        # Tm = 64.9 + 41 * (GC_count - 16.4) / length
        gc_count = seq.count('G') + seq.count('C')
        tm = 64.9 + 41 * (gc_count - 16.4) / length
    
    return tm

def calculate_primer_stats(primer_seq):
    """Calculate GC content and Tm for a primer"""
    gc_count = primer_seq.upper().count('G') + primer_seq.upper().count('C')
    gc_content = (gc_count / len(primer_seq)) * 100
    
    try:
        tm = calculate_tm(primer_seq)
    except:
        tm = None
    
    return {
        'length': len(primer_seq),
        'gc_content': round(gc_content, 1),
        'tm': round(tm, 1) if tm else None
    }

def load_primers_from_excel(filepath):
    """
    Load primers from Excel file in the IDT template format
    Returns DataFrame with primer information
    """
    print(f"\nLoading primers from: {filepath}")
    
    # Read the Excel file starting from row 5 (header row)
    df = pd.read_excel(filepath, header=4)
    
    # Clean column names
    df.columns = ['Oligo_Name', 'Mod_5', 'Sequence', 'Mod_3', 'Scale', 'Purification', 
                  'Format', 'Concentration', 'Number_Tubes', 'Notes'] + [f'Extra_{i}' for i in range(len(df.columns)-10)]
    
    # Extract only rows with actual data (have sequences)
    primers_df = df[df['Sequence'].notna()].copy()
    
    # Extract gene names from primer names (removing _Fwd and _Rev)
    primers_df['Gene'] = primers_df['Oligo_Name'].str.replace('_Fwd', '').str.replace('_Rev', '').str.replace('_fwd', '').str.replace('_rev', '')
    
    print(f"✓ Found {len(primers_df)} primers")
    print(f"✓ Found {primers_df['Gene'].nunique()} unique genes")
    
    return primers_df[['Oligo_Name', 'Sequence', 'Gene']]

def validate_all_primers(primers_df, organism="Mus musculus"):
    """
    Validate all primer pairs in the DataFrame
    Returns results DataFrame
    """
    # Organize primers into pairs
    genes = primers_df['Gene'].unique()
    results = []
    
    print(f"\n{'='*70}")
    print(f"Validating {len(genes)} primer pairs against {organism} genes")
    print(f"{'='*70}\n")
    
    # Cache for gene sequences
    gene_sequences = {}
    
    for i, gene in enumerate(genes, 1):
        print(f"[{i}/{len(genes)}] Checking {gene}...")
        
        # Get forward and reverse primers for this gene
        gene_primers = primers_df[primers_df['Gene'] == gene]
        fwd_primers = gene_primers[gene_primers['Oligo_Name'].str.contains('Fwd|fwd', na=False)]
        rev_primers = gene_primers[gene_primers['Oligo_Name'].str.contains('Rev|rev', na=False)]
        
        if len(fwd_primers) == 0 or len(rev_primers) == 0:
            print(f"  ⚠ Warning: Missing forward or reverse primer for {gene}")
            results.append({
                'Gene': gene,
                'Status': 'ERROR',
                'Amplicon_Size_bp': None,
                'Sequence_Type': 'N/A',
                'Details': 'Missing forward or reverse primer',
                'Fwd_Primer': None,
                'Rev_Primer': None,
                'Fwd_Length': None,
                'Rev_Length': None,
                'Fwd_GC%': None,
                'Rev_GC%': None,
                'Fwd_Tm_C': None,
                'Rev_Tm_C': None
            })
            continue
        
        fwd_seq = fwd_primers.iloc[0]['Sequence']
        rev_seq = rev_primers.iloc[0]['Sequence']
        
        # Fetch gene sequence if not in cache
        if gene not in gene_sequences:
            genomic_seq, cds_seq, gene_header = fetch_gene_sequence(gene, organism)
            gene_sequences[gene] = (genomic_seq, cds_seq, gene_header)
            time.sleep(0.35)  # Be nice to NCBI servers
        else:
            genomic_seq, cds_seq, gene_header = gene_sequences[gene]
            print(f"  Using cached sequences for {gene}")
        
        # Check the primer pair
        result = check_primer_pair(fwd_seq, rev_seq, genomic_seq, cds_seq, gene)
        
        # Calculate primer stats
        fwd_stats = calculate_primer_stats(fwd_seq)
        rev_stats = calculate_primer_stats(rev_seq)
        
        # Print result
        status_symbol = '✓' if result['status'] == 'PASS' else ('⚠' if result['status'] == 'WARN' else '✗')
        
        # Add sequence type to print message
        seq_type_msg = f" [{result.get('sequence_type', 'Unknown')}]" if result.get('sequence_type') else ""
        print(f"  {status_symbol} {result['details']}{seq_type_msg}\n")
        
        # Prepare details with sequence type
        details_with_note = result['details']
        if result.get('sequence_type'):
            details_with_note = f"[{result['sequence_type']}] {details_with_note}"
        
        # Add to results
        results.append({
            'Gene': gene,
            'Status': result['status'],
            'Amplicon_Size_bp': result['amplicon_size'],
            'Sequence_Type': result.get('sequence_type', 'N/A'),
            'Details': details_with_note,
            'Fwd_Primer': fwd_seq,
            'Rev_Primer': rev_seq,
            'Fwd_Length': fwd_stats['length'],
            'Rev_Length': rev_stats['length'],
            'Fwd_GC%': fwd_stats['gc_content'],
            'Rev_GC%': rev_stats['gc_content'],
            'Fwd_Tm_C': fwd_stats['tm'],
            'Rev_Tm_C': rev_stats['tm']
        })
    
    return pd.DataFrame(results)

def print_summary(results_df):
    """Print a summary of the validation results"""
    print(f"\n{'='*70}")
    print("VALIDATION SUMMARY")
    print(f"{'='*70}\n")
    
    total = len(results_df)
    passed = len(results_df[results_df['Status'] == 'PASS'])
    warned = len(results_df[results_df['Status'] == 'WARN'])
    failed = len(results_df[results_df['Status'] == 'FAIL'])
    errors = len(results_df[results_df['Status'] == 'ERROR'])
    
    print(f"Total primer pairs checked: {total}")
    print(f"✓ PASS (70-150 bp):         {passed} ({passed/total*100:.1f}%)")
    print(f"⚠ WARN (wrong size):        {warned} ({warned/total*100:.1f}%)")
    print(f"✗ FAIL (no amplicon):       {failed} ({failed/total*100:.1f}%)")
    print(f"⚠ ERROR (missing primers):  {errors} ({errors/total*100:.1f}%)")
    
    if warned > 0:
        print(f"\n{'='*70}")
        print("PRIMERS WITH INCORRECT AMPLICON SIZE:")
        print(f"{'='*70}\n")
        warn_df = results_df[results_df['Status'] == 'WARN'][['Gene', 'Amplicon_Size_bp', 'Details']]
        for _, row in warn_df.iterrows():
            print(f"  {row['Gene']}: {row['Details']}")
    
    if failed > 0:
        print(f"\n{'='*70}")
        print("PRIMERS THAT FAILED:")
        print(f"{'='*70}\n")
        fail_df = results_df[results_df['Status'] == 'FAIL'][['Gene', 'Details']]
        for _, row in fail_df.iterrows():
            print(f"  {row['Gene']}: {row['Details']}")

def select_input_file():
    """Interactive file selection"""
    print("\n" + "="*70)
    print("PRIMER VALIDATION TOOL")
    print("="*70)
    print("\nPlease enter the path to your Excel file:")
    print("(You can drag and drop the file here, or type the path)")
    filepath = input("\nFile path: ").strip().strip('"').strip("'")
    
    if not os.path.exists(filepath):
        print(f"\n❌ Error: File not found: {filepath}")
        sys.exit(1)
    
    return filepath

def main():
    parser = argparse.ArgumentParser(
        description='Validate primer pairs against mouse genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python check_primers.py                              # Interactive mode
  python check_primers.py primers.xlsx                 # Specify input file
  python check_primers.py primers.xlsx -o results.xlsx # Custom output name
        """
    )
    parser.add_argument('input_file', nargs='?', help='Input Excel file with primers')
    parser.add_argument('-o', '--output', help='Output Excel file name (default: primer_validation_results.xlsx)')
    parser.add_argument('--organism', default='Mus musculus', help='Organism name (default: Mus musculus)')
    
    args = parser.parse_args()
    
    # Get input file
    if args.input_file:
        input_file = args.input_file
        if not os.path.exists(input_file):
            print(f"\n❌ Error: File not found: {input_file}")
            sys.exit(1)
    else:
        input_file = select_input_file()
    
    # Set output file
    if args.output:
        output_file = args.output
    else:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = f"{base_name}_validation_results.xlsx"
    
    # Load primers
    try:
        primers_df = load_primers_from_excel(input_file)
    except Exception as e:
        print(f"\n❌ Error loading Excel file: {str(e)}")
        print("\nMake sure the file is in the correct format (IDT template)")
        sys.exit(1)
    
    # Validate primers
    results_df = validate_all_primers(primers_df, args.organism)
    
    # Print summary
    print_summary(results_df)
    
    # Save results
    print(f"\n{'='*70}")
    print(f"Saving results to: {output_file}")
    
    # Sort by status (FAIL first, then WARN, then PASS)
    status_order = {'FAIL': 0, 'ERROR': 1, 'WARN': 2, 'PASS': 3}
    results_df['_sort'] = results_df['Status'].map(status_order)
    results_df = results_df.sort_values('_sort').drop('_sort', axis=1)
    
    results_df.to_excel(output_file, index=False, sheet_name='Validation Results')
    
    print(f"✓ Results saved successfully!")
    print(f"\n{'='*70}\n")

if __name__ == "__main__":
    main()
