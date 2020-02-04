#!/bin/bash
#BSUB -q q_residual
#BSUB -oo output
#BSUB -eo error
#BSUB -n 128
module load use.own
module load vasp/5.4.4
#BSUB -R "span[ptile=32]"
#BSUB -m "g3_a g3_b"
#BSUB -J Ir5

################################################################################################
#                                    Gets data from input.bh                                   #
################################################################################################

Simbolo_1=$(grep "cluster_ntyp" $1 | cut -d "[" -f 2 | cut -d ":" -f 1)
Simbolo_2=$(grep "cluster_ntyp" $1 | cut -d "[" -f 3 | cut -d ":" -f 1) 2>/dev/null
N_Simbolo_1=$(grep "cluster_ntyp" $1 | cut -d "[" -f 2 | cut -d ":" -f 2 | cut -d "]" -f 1)
N_Simbolo_2=$(grep "cluster_ntyp" $1 | cut -d "[" -f 3 | cut -d ":" -f 2 | cut -d "]" -f 1) 2>/dev/null
n=$(grep "cluster_ntyp"  $1 | awk '{print $4}' | wc -c )  #Este es un criterio para determinar si es bimetálico o no

continue=$(grep "continue" $1 | cut -d " " -f 3 )

randomness=$(grep "randomness" $1 | awk '{print $3}')
kick=$(grep "kick_type" $1 | awk '{print $3}')
file_name=$(grep "file_name" $1 | awk '{print $3}')

Pseudo_Dir=$( grep "pseudo_dir" $1 | awk '{print $3}' )
pseudotype=$(grep "pseudo_type" $1 | awk '{print $3}' )
step_width=$(grep "step_width" $1 | awk '{print $3}')
Temperature=$(grep "temperature_K" $1 | awk '{ print $3 }')
if [ $n -gt 3 ]
then
Nat=$(($N_Simbolo_1+$N_Simbolo_2)) #Numero de atomos del cluster
else
Nat=$(echo $N_Simbolo_1)
fi
Ncore=ncore=$(grep '#BSUB -n' $0 | head -1 | awk '{print $3}' )
iteraciones=$(grep "iterations" $1 | awk '{ print $3 }' )
path=$(grep "initialization_file" $1 | awk '{ print $3 }' )
Npath=$(echo $path | wc -c )
m=1 #contador de coords rechazadas
Sel=$(grep "Selective"  input/POSCAR | wc -l )   #Determina si hay un selective dynamics
NPOSCAR=$(cat input/POSCAR | grep . | wc -l ) #Numero de lineas del poscar sin cluster


swap_step=10      # SWAP STEP

./basin_hopping $continue $randomness $kick $iteraciones $step_width $Temperature $file_name $swap_step $Sel $Nat $n $Ncore $Simbolo_1 $Simbolo_2 $N_Simbolo_1 $N_Simbolo_2
 2>> Ir5.err >> Ir5.out
