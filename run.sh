 bowtie2-build input/mpa_vJan21_TOY_CHOCOPhlAnSGB_202103.fna index/mpa_vJan21_TOY_CHOCOPhlAnSGB_202103.fna

bowtie2 --seed 1992 \
  --quiet --no-unal \
  --very-sensitive \
  -S output/SRS014476-Supragingival_plaque.sam \
  -x ./index/mpa_vJan21_TOY_CHOCOPhlAnSGB_202103.fna \
  -p 4 -U input/SRS014476-Supragingival_plaque.fasta -f
