#!/bin/bash
l=48
d=6
q=6
r=3
n=1000
f="state_prep_entropy_L_${l}_d_${d}_q_${q}_r_${r}_n_${n}"
echo $f
julia --project=. state_prep_plot.jl -f $f -r 2 4 6
