#!/bin/bash
rm result.txt

dossier=OR-Library
 
printf $dossier " \n" >> result.txt

printf " & J & I & v & K & id & Iter & Var & cols I & cols J & CPU & CPU(Master) & Gap & Dual b. & Primal b. & LR & LR(Cplex) & Opt(Cplex) & CPU(OptCplex) \\\\\\ \n " >> result.txt
for v in 10 ; do
  for K in 8 9 10 11 ; do
    for id in "cap61" ; do
      python3 ./Instances/convertinstance.py $dossier $id $v $K
      for met in 3002 3003 ; do
        rm logs/$met.txt
        #rm convergence/${n}_${T}_$id.csv
        ./bin/SCIP_RCFLP_BP.linux.x86_64.gnu.opt.cpx $v $K $dossier $id $met >> logs/$met.txt
        #python3 convergence/plot.py $n $T $id $met
      done
      rm Instances/Reliable/${v}_${K}_${dossier}_$id.txt
      printf "\\hline \n" >> result.txt	
    done
  done
  printf "\\hline \n" >> result.txt	
done