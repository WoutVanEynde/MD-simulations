#MD simulation script for protein with metal ions
#Make sure your protein.pdb file is prepared and energy-minimized

"
Protonated/uncharged Asp			ASH
Protonated/uncharged Glu			GLH
Deprotonated/uncharged Lys			LYN
His protonated at epsilon position		HIE
His protonated at delta position		HID
Charged His (protonated at both positions)	HIP
Deprotonated Cys or Cys bound to a metal	CYM
Cys involved in disulfide bridge		CYX 
"

#STEP 1: create a topology
gmx pdb2gmx -ff amber99sb -water tip3p -ignh -merge all -f *.pdb -o protein.gro -p topol.top -quiet

#STEP 2: Use topol.top and pymol to apply restraints and change 3-letter codon for residues to force deprotonate them
# Open the protein.gro in pymol and look for the interacting residues with the metal ions
# Change these residues their 3 letter codon, f.e. CYS to CYM
 
#STEP 3: add following line into mdp files
disre = Simple   ; distance restraint
nstdisreout = 0  ; MPI 

#STEP 4: create box 
gmx editconf -f protein.gro -o protein_newbox.gro -c -d 1.0 -bt cubic -quiet

#STEP 5: solvate system
gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top -quiet

#STEP 6: add ions
gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr -maxwarn 1 -quiet

gmx genion -s ions.tpr -o protein_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -quiet << EOL
SOL
EOL

#STEP 7: energy minimization and generate potential file
gmx grompp -f minim.mdp -c protein_solv_ions.gro -r protein_solv_ions.gro -p topol.top -o em.tpr -quiet
gmx mdrun -v -deffnm em -quiet
gmx energy -f em.edr -o potential.xvg -quiet << EOL
Potential
EOL

#STEP 8: nvt and generate temperature file
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -quiet
gmx mdrun -v -deffnm nvt -quiet
gmx energy -f nvt.edr -o temperature.xvg -quiet << EOL
Temperature
EOL

#STEP 9: npt and generate pressure and density files
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1 -quiet
gmx mdrun -v -deffnm npt -quiet
gmx energy -f npt.edr -o pressure.xvg -quiet << EOL
Pressure
EOL
gmx energy -f npt.edr -o density.xvg -quiet << EOL
Density
EOL

#STEP 10: generate tpr file on computer
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr -quiet

#STEP 11: run MD
gmx mdrun -v -deffnm md -quiet -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4
