#!/bin/bash
rm result.txt

dossier=OR-Library
 
printf $dossier " \n" >> result.txt

printf " & J & I & v & K & id & Iter & Var & cols I & cols J & CPU & CPU(Master) & Gap & Dual b. & Primal b. & LR & LR(Cplex) & Opt(Cplex) & CPU(OptCplex) \\\\\\ \n " >> result.txt
for v in {1..10} ; do
  for K in {0..2*$v}
    for id in "cap61" ; do
      ./Instances
      for met in 100 200 3000 3001 ; do
        rm logs/$met.txt
        #rm convergence/${n}_${T}_$id.csv
        ./bin/SCIP_RCFLP_BP.linux.x86_64.gnu.opt.cpx $v $K $dossier $id $met >> logs/$met.txt
        #python3 convergence/plot.py $n $T $id $met
      done
      printf "\\hline \n" >> result.txt	
    done
    printf "\\hline \n" >> result.txt	
done