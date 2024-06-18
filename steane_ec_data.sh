#!/bin/bash

l=48
dlist=(2 4 6)
qlist=(2 4 6)
r=3
for i in 0 1 2
do
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.01 -n 4000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.015 -n 4000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.0175 -n 4000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.02 -n 13000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.021 -n 13000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.022 -n 13000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.023 -n 13000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.024 -n 13000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.025 -n 13000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.0275 -n 4000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.03 -n 4000 -x 0 -z 1 -h"
  sbatch -n 1 --time=72:00:00 --mem=8gb --wrap="srun -n 1 julia ./steane_ec_data.jl -l $l -d ${dlist[i]} -q ${qlist[i]} -r $r -p 0.035 -n 4000 -x 0 -z 1 -h"
done
