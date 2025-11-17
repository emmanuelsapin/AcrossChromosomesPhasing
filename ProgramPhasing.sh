#!/bin/sh
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --qos=blanca-ibg
#SBATCH --mem=120G
 

echo Welcome to the script 

NbIndiv=$1
PathInput=$2  
PathOutput=$3


#  sbatch /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/ProgramPhasing.sh 1000 /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/data/chr /pl/active/KellerLab/Emmanuel/gameticphasing/Fileforgithub/
#./GenerateRandomData data/chr HAP
mkdir "$PathOutput"
 
./ProgramPhasing -NbIndiv "$NbIndiv" -PathInput "$PathInput" -PathOutput "$PathOutput"
 