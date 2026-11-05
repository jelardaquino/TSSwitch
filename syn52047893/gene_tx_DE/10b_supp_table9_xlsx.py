"""Write Supplementary Table 9 (.xlsx, one sheet) from the TSV made by 10_supp_table9_gene_tx_DE.R.

Run with openpyxl on PYTHONPATH (not with -I, which ignores it):  PYTHONPATH=<dir with openpyxl> python3 10b_supp_table9_xlsx.py
"""
import csv
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

TABLES = Path("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/emtab14167_replication/tables")
SRC = TABLES / "SuppTable9_gene_tx_DE_DIU.tsv"
OUT = TABLES / "Supplementary_Table_9.xlsx"

MODELS = {"noCov": "~ condition"}
BASE = {
    "gene_id": "Ensembl gene ID",
    "gene_name": "Gene symbol",
    "DIU": "DIU gene",
    "switch_ids": "Switching isoforms (dIF in AD)",
    "promoter_switching_gene": "Promoter-switching gene",
}
PER_MODEL = {
    "log2FC": "Gene log2FC (shrunken)",
    "p": "Gene p",
    "padj": "Gene FDR",
    "geneDE": "Gene DEG",
    "n_txDE": "DE transcripts (n)",
    "txDE_ids": "DE transcripts (direction in AD)",
    "region": "Venn region",
}
NUMERIC = {"log2FC", "p", "padj", "n_txDE"}

TITLE = ("Supplementary Table 9. Genes in Supplementary Figure 5: gene-level and transcript-level differential "
         "expression and differential isoform usage in the long-read data (6 AD vs 6 CT).")
NOTES = [
    "One row per gene that is a gene DEG, has at least one differentially expressed transcript, or is one of the 86 DIU genes. "
    "DESeq2 (~ condition) on Bambu gene and transcript counts; DE = FDR < 0.05 and |apeglm-shrunken log2FC| > 1. "
    "DIU: IsoformSwitchAnalyzeR, FDR < 0.05, |dIF| > 0.1. Positive log2FC and dIF = higher in AD.",
    "'not tested' = removed by the expression filter (≥ 10 counts in ≥ 3 samples of either group) or p-value set to NA by "
    "DESeq2. Filter on 'Venn region' to reproduce each region of Supplementary Figure 5.",
]
HEADER_ROW = len(NOTES) + 3

HEADER_FILL = PatternFill("solid", fgColor="DDE5F0")
MODEL_FILL = {"noCov": PatternFill("solid", fgColor="EEF2F7"), "withCov": PatternFill("solid", fgColor="F5EFE6")}

with open(SRC, newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

cols = [(k, h, None, k) for k, h in BASE.items()]
for m, design in MODELS.items():
    cols += [(f"{m}_{k}", h, m, k) for k, h in PER_MODEL.items()]

wb = Workbook()
ws = wb.active
ws.title = "Supplementary Table 9"

ws.cell(row=1, column=1, value=TITLE).font = Font(bold=True)
for i, note in enumerate(NOTES, start=2):
    ws.cell(row=i, column=1, value=note)

for j, (_, header, m, _) in enumerate(cols, start=1):
    c = ws.cell(row=HEADER_ROW, column=j, value=header)
    c.font = Font(bold=True)
    c.fill = MODEL_FILL[m] if m else HEADER_FILL
    c.alignment = Alignment(wrap_text=True, vertical="top")
ws.row_dimensions[HEADER_ROW].height = 60

for i, r in enumerate(rows, start=HEADER_ROW + 1):
    for j, (key, _, _, short) in enumerate(cols, start=1):
        raw = r[key]
        if raw == "":
            continue
        if short in NUMERIC:
            c = ws.cell(row=i, column=j, value=int(float(raw)) if short == "n_txDE" else float(raw))
            if short in {"p", "padj"}:
                c.number_format = "0.00E+00"
            elif short == "log2FC":
                c.number_format = "0.000"
        else:
            ws.cell(row=i, column=j, value=raw)

for j, (key, _, _, short) in enumerate(cols, start=1):
    ws.column_dimensions[get_column_letter(j)].width = 40 if short in {"switch_ids", "txDE_ids"} else 15
ws.freeze_panes = ws.cell(row=HEADER_ROW + 1, column=3)
ws.auto_filter.ref = f"A{HEADER_ROW}:{get_column_letter(len(cols))}{HEADER_ROW + len(rows)}"

wb.save(OUT)
print("wrote", OUT, "|", len(rows), "genes |", len(wb.sheetnames), "sheet")
