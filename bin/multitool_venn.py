#!/usr/bin/env python3

import argparse
import pandas as pd
import subprocess
from pathlib import Path


TOOL_SCORE_COLUMN = {
    'vs2': 'vs2_score',
    'vb': 'vb_score',
    'dvf': 'dvf_score',
    'gn': 'gn_score'
}

def main():
    parser = argparse.ArgumentParser(
        description='Generate Venn diagram for viral contigs identified by multiple tools'
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='all_viral_ctgs.score.tsv'
    )
    parser.add_argument(
        '-t', '--tools', required=True,
        help='Comma-separated tools: vs2,vibrant,dvf,gn (>=2 tools)'
    )
    parser.add_argument(
        '-o', '--outdir', required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--venn_dir', default=str(Path(__file__).parent.resolve()),
        help='Directory containing venn2.R / venn3.R / venn4.R'
    )
    args = parser.parse_args()
    tools = [t.strip() for t in args.tools.split('-')]
    if not (2 <= len(tools) <= 4):
        raise ValueError('Number of tools must be between 2 and 4')
    for t in tools:
        if t not in TOOL_SCORE_COLUMN:
            raise ValueError(f'Unsupported tool: {t}')
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Load table
    df = pd.read_csv(args.input, sep='\t')
    contig_col = 'Contig'
    if contig_col not in df.columns:
        raise ValueError("Column 'Contig' not found in input file")
    list_files = []

    # Generate contig lists for each tool
    for tool in tools:
        score_col = TOOL_SCORE_COLUMN[tool]
        if score_col not in df.columns:
            raise ValueError(f"Column '{score_col}' not found in input file")
        out_list = outdir / f'{tool}.contigs'
        df.loc[df[score_col] > 0, contig_col].to_csv(
            out_list, index=False, header=False
        )
        list_files.append(str(out_list))

    # Select Venn script
    venn_script = Path(args.venn_dir) / f'venn{len(tools)}.R'
    if not venn_script.exists():
        raise FileNotFoundError(f'{venn_script} not found')
    pdf_out = outdir / f'viral_contigs_{args.tools}_comparison_venn.pdf'
    cmd = ['Rscript', str(venn_script)] + list_files + [str(pdf_out)]
    subprocess.run(cmd, check=True)
    print('Venn diagram generated:', pdf_out)

if __name__ == '__main__':
    main()
