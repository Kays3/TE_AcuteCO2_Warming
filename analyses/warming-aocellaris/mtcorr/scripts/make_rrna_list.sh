#!/usr/bin/env bash
# Amendment 5, A5.1: every gene with gene_biotype "rRNA" in the RefSeq GTF used for counting (GCF_022539595.1), one gene_id per line.
# Usage: make_rrna_list.sh GTF > rrna_genes.txt
awk -F'\t' '$3 == "gene" && $9 ~ /gene_biotype "rRNA"/' "$1" | sed -E 's/.*gene_id "([^"]+)".*/\1/' | sort -u
