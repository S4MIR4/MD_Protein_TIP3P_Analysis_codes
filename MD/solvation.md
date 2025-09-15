;define the box
````bash
gmx editconf -f name_processed.gro -o name_newbox.gro -c -d 1.2 -bt cubic

;solvate the box with water
````bash
gmx solvate -cp name_newbox.gro -cs spc216.gro -o name_solv.gro -p topol.top
