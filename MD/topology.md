;generate topology by the .pdb file of interest, first download your protein from RCSB website in .pdb format, if your structure contains water first remove the water 

````bash
grep -v HOH name.pdb > name_clean.pdb

; then you need to download your forcefiles and keep the forcefield file in the directory that you intend to run MD

```bash
export GMXLIB=$pwd/(YOUR FORCEFIELD)

;to generate topology file from a .pdb file: 

```bash
gmx pdb2gmx -f name_clean.pdb -o name_processed.gro -water tip3p

 
