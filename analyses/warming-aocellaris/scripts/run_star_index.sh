#!/usr/bin/env bash
cd ~/jim_te_warm
F=~/KD/workspace/envs/te-ocean-acidification/bin
R=ref/GCF_022539595.1_ASM2253959v1
exec /usr/bin/time -v -o logs/star_index.time taskset -c 0-20 $F/STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir ref/star_149 --genomeFastaFiles ${R}_genomic.fna --sjdbGTFfile ${R}_genomic.gtf --sjdbOverhang 149 --genomeSAindexNbases 13 \
  --outTmpDir ~/jim_te_warm/tmp_star_index > logs/star_index.log 2>&1 < /dev/null
