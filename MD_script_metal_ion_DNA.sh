#MD simulation script for protein with DNA and metal ions.
#Make sure your protein.pdb file is prepared and energy-minimized.
#The .gro and .g96 file formats do not support chain identifiers. Therefore it is useful to enter a .pdb file name at the -o option when you want to convert a multi-chain .pdb file. 

#STEP 1: create a topology
gmx pdb2gmx -ff amber99sb -water tip3p -ignh -merge all -f *.pdb -o protein.gro -p topol.top -quiet

#STEP 2: Use topol.top and pymol to apply restraints and change 3-letter codon for residues to force deprotonate them
# Open the protein.gro in pymol and look for the interacting residues with the metal ions
# Change these residues their 3 letter codon, f.e. CYS to CYM
# Create a table like this in the topology file before the line "; Include Position restraint file" (more info at https://manual.gromacs.org/2021-current/reference-manual/functions/restraints.html?highlight=restraints):

"
#For our androgen receptor in particular:
#DC1 H41	DG21 06
#DC1 N3		DG21 H1
#DC1 02 	DG21 H21

[ distance_restraints ]
; ai	aj	type	index	type'	low	up1	up2	fac
2196	3624    1	4	2	0.0	0.3	0.4	1.0
2198	3626    1	4	2	0.0	0.3	0.4	1.0
2200	3629    1	4	2	0.0	0.3	0.4	1.0
2925	2892    1	5	2	0.0	0.3	0.4	1.0
2927	2894    1	5	2	0.0	0.3	0.4	1.0
2929	2897	1	5	2	0.0	0.3	0.4	1.0
" 

#STEP 3: add following line into minim.mdp, nvt.mdp, npt.mdp and md.mdp files
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

#STEP 7: energy minimization and generate potential file (CHARMM-GUI START)
gmx grompp -f minim.mdp -c protein_solv_ions.gro -r protein_solv_ions.gro -p topol.top -o em.tpr -quiet
gmx mdrun -v -deffnm em -quiet -nb gpu
gmx energy -f em.edr -o potential.xvg -quiet << EOL
Potential
EOL

#STEP 8: create group
gmx make_ndx -f em.gro -o index.ndx -quiet << EOL
1|12|13
14|15|16
quit
EOL

#STEP 9: check group
gmx make_ndx -f em.gro -o index2.ndx -n index.ndx -quiet << EOL
1|13
quit
EOL

#STEP 10: modify tc-grps in nvt.mdp, npt.mdp, md.mdp
#tc-grps = Protein_DNA_ZN Na+_Cl-_OPC

#STEP 11: nvt and generate temperature file
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx -quiet
gmx mdrun -v -deffnm nvt -quiet -nb gpu
gmx energy -f nvt.edr -o temperature.xvg -quiet << EOL
Temperature
EOL

#STEP 12: npt and generate pressure and density files
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -n index.ndx -quiet
gmx mdrun -v -deffnm npt -quiet -nb gpu
gmx energy -f npt.edr -o pressure.xvg -quiet << EOL
Pressure
EOL
gmx energy -f npt.edr -o density.xvg -quiet << EOL
Density
EOL

#STEP 13: generate tpr file on computer
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr -n index.ndx -quiet

#STEP 14: run MD; adjust -ntmpi and gputasks if given error
gmx mdrun -v -deffnm md -quiet -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4

#STEP 15: continue MD
gmx mdrun -deffnm md -quiet -v -cpi md.cpt -append -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4

#STEP 16: append MD (picoseconds)
gmx convert-tpr -s md.tpr -extend 300000 -o md2.tpr
gmx mdrun -v -quiet -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4 -s md2.tpr -cpi md.cpt -deffnm md
