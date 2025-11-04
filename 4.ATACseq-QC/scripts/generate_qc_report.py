#!/usr/bin/env python3
"""
generate_qc_report.py
Combine QC metrics and plots into a summarized ATAC-seq QC PDF report.
"""
import sys, os
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate ATAC-seq QC summary PDF report.")
    parser.add_argument('--input', required=True, help='Path to qc_results directory')
    parser.add_argument('--output', required=True, help='Output PDF path')
    args = parser.parse_args()

    pdf = SimpleDocTemplate(args.output)
    styles = getSampleStyleSheet()
    story = [Paragraph("ATAC-seq Quality Control Report", styles['Title']), Spacer(1, 12)]

    metrics = [
        ['Metric', 'Value', 'Threshold'],
        ['TSS Enrichment', '-', '>5'],
        ['Mitochondrial %', '-', '<50%'],
        ['NRF', '-', '>0.5'],
        ['PBC', '-', '>0.5'],
        ['FRiP', '-', '>0.1'],
        ['IDR', '-', '<0.2'],
        ['Correlation', '-', '>0.7']
    ]

    table = Table(metrics, hAlign='LEFT')
    table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.grey),
                               ('TEXTCOLOR',(0,0),(-1,0),colors.whitesmoke),
                               ('ALIGN',(0,0),(-1,-1),'CENTER'),
                               ('GRID',(0,0),(-1,-1),1,colors.black)]))
    story.append(table)
    story.append(Spacer(1, 12))
    story.append(Paragraph("Visual QC plots and interpretation will be inserted here.", styles['Normal']))

    pdf.build(story)
    print(f"QC report saved to {args.output}")

if __name__ == '__main__':
    main()
