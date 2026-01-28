#!/bin/bash
# denovo.sh
# Parameters:
# -M : maximum distance allowed between stacks (ustacks)
# -n : number of mismatches allowed between sample loci (cstacks)
# -T : number of threads

denovo_map.pl -M 3 -n 3 -T 8 \
  --samples /home/redempta/Redempta/HaploNet/denovo/test \
  --popmap ./populations_tilapia.txt \
  -o ./stacks_output \
  -X "gstacks: --var-alpha 0.001 --gt-alpha 0.001 --min-mapq 40"

