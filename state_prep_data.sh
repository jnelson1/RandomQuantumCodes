#!/bin/bash

l=48
d=6
q=6
r=3
n=1000
f="state_prep_entropy_L_${l}_d_${d}_q_${q}_r_${r}_n_${n}"
julia ./state_prep_data.jl -l $l -d $d -q $q -r $r -n $n -f $f -p 0.005 -h
for p in 0.0075 0.01 0.0125 0.015 0.0175 0.02 0.0225 0.025 0.0275 0.03 0.0325 0.035;
do
  julia ./state_prep_data.jl -l $l -d $d -q $q -r $r -n $n -f $f -p $p
done
