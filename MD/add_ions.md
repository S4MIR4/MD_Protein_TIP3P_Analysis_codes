```
;ions.mdp
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching 
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions

```


;assemble your .tpr  (binary file includes the atomic_level description of our system)
````bash
gmx grompp -f ions.mdp -c name_solv.gro -p topol.top -o ions.tpr

;to generate the ions: add negative ions if your total charge of your system is posotive and vice versa
````bash
gmx genion -s ions.tpr -o name_solv_ions.gro -p topol.top -pname NA -nname CL -neutral


; check the topol.top [ molecules ] section to make sure the solvation and ion section went through.
