#MD simulation script for protein with metal ions bound to ligand.
#acpype (https://github.com/alanwilter/acpype).
#Make sure your protein.pdb file is prepared and energy-minimized.
#The .gro and .g96 file formats do not support chain identifiers. Therefore it is useful to enter a .pdb file name at the -o option when you want to convert a multi-chain .pdb file.
#L1.mol2 and L2.mol2 files are in the working directory and are named properly in second line and in second last column!

#STEP 1: prepare ligand file with -n defining the charge
acpype -i L1.mol2 -n 0
acpype -i L2.mol2 -n 0

#STEP 2: use topol.top to change 3-letter codon for residues to force deprotonate them

#STEP 3: create a topology
gmx pdb2gmx -ff amber99sb -water tip3p -ignh -merge all -f *.pdb -o protein.pdb -p topol.top -quiet

#STEP 4: copy ligand2_NEW.pdb and ligand1_GMX.itp from the folder created by acpype
cp ../L1.acpype/L1_NEW.pdb .
cp ../L1.acpype/L1_GMX.itp .
mv L1_GMX.itp L1.itp
cp ../L2.acpype/L2_NEW.pdb .
cp ../L2.acpype/L2_GMX.itp .
mv L2_GMX.itp L2.itp

#STEP 5: add ligand coordinate file to protein.gro
grep -h ATOM protein.pdb L1_NEW.pdb >| complex.pdb
grep -h ATOM complex.pdb L2_NEW.pdb >| complex2.pdb
rm complex.pdb
mv complex2.pdb complex.pdb

#STEP 6: Distance restraints on H-bonds on end of DNA termini

"
#For our androgen receptor in particular:
#DC1 H41	DG21 06
#DC1 N3		DG21 H1
#DC1 02 	DG21 H21

[ distance_restraints ]
; ai	aj	type	index	type'	low	up1	up2	fac
6295	7722    1	4	2	0.0	0.3	0.4	1.0
6297	7724    1	4	2	0.0	0.3	0.4	1.0
6299	7727    1	4	2	0.0	0.3	0.4	1.0
7028	6995    1	5	2	0.0	0.3	0.4	1.0
7030	6997    1	5	2	0.0	0.3	0.4	1.0
7032	7000	1	5	2	0.0	0.3	0.4	1.0 
" 

#STEP 7: include ligand topology to topol.top!
cat topol.top | sed '/water\ topology/i\#include "L1.itp"' >| topol_first_complex.top
cat topol_first_complex.top | sed '/water\ topology/i\#include "L2.itp"' >| topol_second_complex.top
rm topol_first_complex.top
mv topol_second_complex.top topol_first_complex.top

#step 8: add ligand reference to topol.top, should be same name as .itp file!
(cat topol_first_complex.top ; echo "L1          1") > topol_second_complex.top
(cat topol_second_complex.top ; echo "L2          1") > topol_third_complex.top
rm topol_second_complex.top
mv topol_third_complex.top topol_second_complex.top

#STEP 9: create gaff file, this is forcefield for the ligand; if multiple, different ligands, merge the gaff files manually! 
grep -Eh 'atomtypes|sigma|0.00000  0.00000' L1.itp >| gaff.ff
grep -Eh 'atomtypes|sigma|0.00000  0.00000' L2.itp >| gaff2.ff

#STEP 10: include ligand gaff information reference to topol.top
cat topol_second_complex.top | sed '/forcefield\.itp\"/a\#include "gaff.ff"' >| topol_complex.top

#STEP 11: delete ligand gaff information from ligand itp file
grep -Ev 'atomtypes|sigma|0.00000  0.00000' L1.itp >| L1_2.itp
grep -Ev 'atomtypes|sigma|0.00000  0.00000' L2.itp >| L2_2.itp

#STEP 12: clean up and change name
rm topol.top topol_first_complex.top topol_second_complex.top 
rm L1.itp
mv L1_2.itp L1.itp
rm L2.itp
mv L2_2.itp L2.itp

#STEP 13: Generate a file analogous to posre.itp for our ligands
gmx genrestr -f L1_NEW.pdb -o posre_L1.itp -fc 1000 1000 1000 << EOL
L1
EOL
gmx genrestr -f L2_NEW.pdb -o posre_L2.itp -fc 1000 1000 1000 << EOL
L2
EOL

#STEP 14: Include position restraint reference of proteins and ligands to topol.top
cat topol_complex.top | sed '/L1\.itp\"/a\#ifdef POSRES\n#include "posre_L1.itp"\n#endif' >| topol_complex2.top
mv topol_complex2.top topol_complex.top
cat topol_complex.top | sed '/L2\.itp\"/a\#ifdef POSRES\n#include "posre_L2.itp"\n#endif' >| topol_complex2.top
mv topol_complex2.top topol_complex.top

#STEP 15: create box 
gmx editconf -f complex.pdb -o complex_newbox.gro -c -d 1.0 -bt cubic -quiet

#STEP 16: solvate system
gmx solvate -cp complex_newbox.gro -cs spc216.gro -o complex_solv.gro -p topol_complex.top -quiet

#STEP 17: create a index group that merges the protein and compound
gmx make_ndx -f complex_solv.gro -o index.ndx -quiet << EOL
1|12|14
quit
EOL

#STEP 18: check group
gmx make_ndx -f complex_solv.gro -o index2.ndx -n index.ndx -quiet << EOL
quit
EOL

#STEP 19: run centralization
gmx grompp -f ions.mdp -c complex_solv.gro -p topol_complex.top -o topol_complex_centered.tpr -quiet -maxwarn 1
gmx trjconv -f complex_solv.gro -o complex_solv_centered.gro -pbc mol -ur compact -center -s topol_complex_centered.tpr -n index.ndx -quiet << EOL
20
System
EOL

#STEP 20: add ions
gmx grompp -f ions.mdp -c complex_solv_centered.gro -p topol_complex.top -o ions.tpr -maxwarn 1 -quiet
gmx genion -s ions.tpr -o complex_solv_centered_ions.gro -p topol_complex.top -pname NA -nname CL -neutral -quiet << EOL
SOL
EOL

#STEP 21: energy minimization and generate potential file (CHARMM-GUI START)
gmx grompp -f minim.mdp -c complex_solv_centered_ions.gro -r complex_solv_centered_ions.gro -p topol_complex.top -o em.tpr -quiet
gmx mdrun -nb gpu -v -deffnm em -quiet -gpu_id 0
gmx energy -f em.edr -o potential.xvg -quiet << EOL
Potential
EOL

#STEP 21: new index file
gmx make_ndx -f em.gro -o index.ndx -quiet << EOL
1|12|13|14
15|16|17
q
EOL

#STEP 22: check group
gmx make_ndx -f complex_solv.gro -o index2.ndx -n index.ndx -quiet << EOL
quit
EOL

#STEP 23: modify tc-grps in nvt.mdp, npt.mdp, md.mdp
#tc-grps = Protein_DNA_DHT_ZN Na+_Cl-_OPC

#STEP 24: nvt and generate temperature file
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol_complex.top -o nvt.tpr -n index.ndx -quiet
gmx mdrun -nb gpu -v -deffnm nvt -quiet -gpu_id 0
gmx energy -f nvt.edr -o temperature.xvg -quiet << EOL
Temperature
EOL

#STEP 25: npt and generate pressure and density files
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol_complex.top -o npt.tpr -n index.ndx -quiet
gmx mdrun -nb gpu -v -deffnm npt -quiet -gpu_id 0
gmx energy -f npt.edr -o pressure.xvg -quiet << EOL
Pressure
EOL
gmx energy -f npt.edr -o density.xvg -quiet << EOL
Density
EOL

#STEP 26: generate tpr file on computer
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol_complex.top -o md.tpr -n index.ndx -quiet

#STEP 27: run MD; adjust -ntmpi and gputasks if given error
gmx mdrun -v -deffnm md -quiet -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4

#STEP 28: continue MD
gmx mdrun -deffnm md -quiet -v -cpi md.cpt -append -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4

#STEP 29: append MD (picoseconds)
gmx convert-tpr -s md.tpr -extend 300000 -o md2.tpr
gmx mdrun -v -quiet -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4 -s md2.tpr -cpi md.cpt -deffnm md
