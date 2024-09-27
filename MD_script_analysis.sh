#Recentering and rewrapping coordinates of the trajectory files

#STEP 1: make index file
gmx make_ndx -f md.gro -o index.ndx -quiet << EOL
1|12|13
q
EOL

#STEP 2: makes broken protein whole
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC_whole.xtc -pbc whole -quiet -n index.ndx << EOL
17
EOL

#STEP 3: eliminates jumps from one side to the other side of the box
gmx trjconv -s md.tpr -f  md_noPBC_whole.xtc -pbc nojump -o md_noPBC_whole_nojump.xtc -quiet -n index.ndx << EOL
17
EOL

#STEP 4: center protein in box
gmx trjconv -s md.tpr -f md_noPBC_whole_nojump.xtc -o md_noPBC_whole_nojump_center.xtc -center -quiet -n index.ndx << EOL
Protein
17
EOL

#STEP 5: center all molecules in box and put all atoms at the closest distance from center of the box
gmx trjconv -s md.tpr -f md_noPBC_whole_nojump_center.xtc -o md_noPBC_whole_nojump_center_mol_com.xtc -pbc mol -ur compact -quiet -n index.ndx << EOL
17
EOL

#STEP 6: fit the system to reference in structure file
gmx trjconv -s md.tpr -f md_noPBC_whole_nojump_center_mol_com.xtc -o md_noPBC_whole_nojump_center_mol_com_fit.xtc -fit rot+trans -quiet -n index.ndx << EOL
Backbone
17
EOL
rm md_noPBC_whole.xtc md_noPBC_whole_nojump.xtc md_noPBC_whole_nojump_center.xtc md_noPBC_whole_nojump_center_mol_com.xtc

#STEP 6: output a .pdb file that contains the exact particles as the trajectory at t=0, so other programs know how to interpret the numbers in the trajectory
gmx trjconv -s md.tpr -f md.xtc -o md.pdb -pbc mol -ur compact -dump 0 -quiet -n index.ndx << EOL
17
EOL

#STEP 8: make indeces of chains
gmx make_ndx -f md.gro -o index2.ndx -quiet << EOL
splitch 1
17|18
q
EOL

#STEP 9: visualise using VMD or pymol
pymol md.gro md_noPBC_whole_nojump_center_mol_com_fit.xtc

#STEP 10: calculate RMSD
gmx rms -s md.tpr -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -o rmsd.xvg -tu ns -quiet << EOL
Backbone
Backbone
EOL

#STEP 11: visualise structural stability
gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set output 'rmsd.png'
set title 'RMSD'
set ylabel 'RMSD (nm)'
set xlabel 'Time (ns)'
unset key  
set grid
plot 'rmsd.xvg' u 1:2 w lines lw 2 lc rgb '#27ad81'
quit
EOL
okular rmsd.png

gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set output 'RMSD DR3.png'
set title 'RMSD DR3'
set ylabel 'RMSD (nm)'
set xlabel 'Time (ns)'
set xrang [0:200]
set grid
set key top right
plot './2_first_run/rmsd.xvg' using 1:2 smooth sbezier w lines lw 2 lc rgb 'red' title 'First simulation', \
     './2_second_run/rmsd.xvg' using 1:2 smooth sbezier w lines lw 2 lc rgb 'orange' title 'Second simulation', \
     './2_third_run/rmsd.xvg' using 1:2 smooth sbezier w lines lw 2 lc rgb 'yellow' title 'Third simulation'
quit
EOL

#STEP 12: visualise hydrogen bonds
gmx hbond -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -s md.tpr -n index.ndx << EOL
1
12
EOL

gnuplot << EOL 
set terminal pngcairo enhanced font 'Verdana,10'
set output "hbnum.png"
set xlabel "Time (ps)"
set ylabel "Number of Hydrogen Bonds"
set title "Hydrogen Bond Formation"
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5
set yrange [0:*]     # Set y-axis minimum to 0 and maximum to auto
datafile = "hbnum.xvg"
stats datafile using 2 nooutput
avg = STATS_mean
plot datafile using 1:2 with lines ls 1 title "Number of Hydrogen Bonds", \
     avg with lines lt -1 lw 2 title sprintf("Average = %.2f", avg)
EOL

#STEP 13: B-factor visualisation
gmx rmsf -oq B_factor.pdb -ox B_factor_average_trajectory.pdb -f md_mmpbsa.xtc -s md.tpr -n index.ndx << EOL
20
EOL
rm B_factor_average_trajectory.pdb

gmx rmsf -f md_mmpbsa.xtc -s md.tpr -n index.ndx -o rmsf_backbone.xvg<< EOL
3
EOL

gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set output 'rmsf.png'
set title 'RMSF, backbone'
set ylabel 'RMSD (nm)'
set xlabel 'Residue number'
unset key  
set grid
plot 'rmsf_backbone.xvg' u 1:2 w lines lw 2 lc rgb '#27ad81'
quit
EOL

#SAUSAGE PLOT: Adjust in MOE if file messed up
show cartoon
cartoon putty
unset cartoon_smooth_loops
unset cartoon_flat_sheets
set ray_shadows,0
set ray_trace_mode,1
set ray_trace_color, black
spectrum b, rainbow, minimum=0, maximum=100
ramp_new color_bar, none, [0, 100], rainbow
sele resn ZN
cmd.color("gray", "sele")
sele resn 0YB
cmd.hide("everything","sele")
sele resn 0YA
cmd.hide("everything","sele")
sele resn ROH
cmd.hide("everything","sele")
set ray_shadows,0
set ray_trace_mode,1
set ray_trace_color, black
set ray_trace_gain, 0.005

# extra
for i in range(1, 80): cmd.bond(f"/B_factor//C/{i}/O3'", f"/B_factor//C/{i+1}/P")
set label_color, black, B_factor

# B_factor_averaging:
awk '/^ATOM/ && ($3=="CA" || $3=="C4\x27") {chain=substr($0,22,1); resnum=substr($0,23,4); print chain resnum}' B_factor.pdb > extracted_chain_and_residue.txt
awk '/^ATOM/ && ($3=="CA" || $3=="C4\x27") {print substr($0,61,7)}' B_factor.pdb > extracted_B_factor.txt
paste extracted_B_factor*.txt > merged_B_factor.txt
awk '{printf "%s %.2f\n", $0, ($1+$2+$3)/3}' merged_B_factor.txt > averaged_B_factor.txt
awk '{print $4}' averaged_B_factor.txt > averaged_B_factor2.txt
sed -i 's/\,/./g' averaged_B_factor2.txt
paste extracted_chain_and_residue.txt averaged_B_factor2.txt > input_insertion_py.txt
vim -c "%s/NGLNA/NGL A/g | %s/NGLNB/NGL B/g | %s/CLEUA/CLE A/g | %s/CLEUB/CLE B/g" -c "wq" B_factor.pdb
python /home/wout/Scripts/1_molecular_dynamics/insert-bfactor.py -p B_factor.pdb -b input_insertion_py.txt -o B_factor_averaged.pdb -d 0.0
vim -c "%s/NGL A/NGLNA/g | %s/NGL B/NGLNB/g | %s/CLE A/CLEUA/g | %s/CLE B/CLEUB/g" -c "wq" B_factor_averaged.pdb

rm B_factor.pdb
rm extracted_B_factor*.txt extracted_chain_and_residue.txt
rm averaged_B_factor.txt merged_B_factor.txt averaged_B_factor2.txt

#STEP 14: angle visulisation
gmx make_ndx -f md.pdb -o index_angle3.ndx -quiet << EOL
ri 26 & 3
ri 43 & 3
ri 97 & 3
ri 114 & 3
15 | 16 | 17
17 | 18 | 15
name 19 A564_A581_B564
name 20 B564_B581_A564
ri 66 & 3
ri 137 & 3
ri 148 & a P
ri 161 & a P
16 | 15 | 18
15 | 16 | 17
name 19 B604_A604_C19
name 20 A564_B604_C6
ri 47 & 3
ri 177 & a P
ri 118 & 3
ri 154 & a P
15 | 18 | 17
15 | 16 | 17
name 19 A585T_C12_B585T
name 20 A585T_D12_B585T
quit
EOL

gmx angle -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -n index_angle3.ndx -od bend5av.xvg -ov bend5.xvg -quiet << EOL
19
EOL
gmx angle -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -n index_angle3.ndx -od bend6av.xvg -ov bend6.xvg -quiet << EOL
20
EOL

gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set datafile separator ','
set output 'averaged_bending_over_time.png'
set title 'Bending'
set ylabel 'Angle (degree)'
set xlabel 'Frames'  
set grid
plot 'angle1_averaged.csv' u 2 with lines lw 2 lc rgb 'red' title 'A564 A581 B564', \
     'angle2_averaged.csv' u 2 with lines lw 2 lc rgb 'green' title 'B564 B581 A564'
set key
quit
EOL

#STEP 15: distance calculations (Use MOE for residue indeces)
gmx make_ndx -f md.pdb -o index_distance.ndx -quiet << EOL 
ri 29 & a NZgmx make_ndx -f md.tpr -o index_mmpbsa << EOL
1 | 12 | 13
14 | 15| 16
1 | 13
ri 144-164
ri 167-187
21 | 22
name 23 stable_DNA
20 | 23
name 24 Protein_stable_DNA_ZN
ri 147-161
ri 170-184
25 | 26
name 27 Protein_truncated_DNA_ZN
q
EOL
ri 33 & a CD 
ri 36 & a NE2
ri 38 & a OH  
ri 148 & a OP2
ri 149 & a OP2
ri 178 & a OP1
ri 179 & a OP1
15 | 19
15 | 20
15 | 16
17 | 21
18 | 22 
name 23 A567K_C4G
name 24 A567K_C5T
name 25 A567K_A571E
name 26 A574Q_D11A
name 27 A576Y_D12C
ri 41 & a O
ri 42 & a OG
ri 47 & a OG1
ri 112 & a O
ri 113 & a OG
ri 118 & a OG1
28 | 33
29 | 32
30 | 31
name 34 A579A_B585T
name 35 A580S_B580S
name 36 A585T_B579A
ri 96 & a NZ
ri 101 & a CZ
ri 157 & a O6
ri 171 & a O6
ri 25 & a NZ
ri 30 & a CZ
ri 148 & a O6
ri 180 & a O6
37 | 40
38 | 39
41 | 43
42 | 44
name 45 CAN_B563K_D6G
name 46 CAN_B568K_C15G
name 47 NON_CAN_A563K_C6G
name 48 NON_CAN_A568R_D15G    
q
EOL

gmx distance -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -s md.pdb -n index_distance.ndx -select 23 24 25 26 27 34 35 36 45 46 47 48 -oallstat -oall -quiet

gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set datafile separator ","
set output 'A585T B579A.png'
set title 'A585T B579A'
set ylabel 'Distance (nm)'
set xlabel 'Time (ps)'
set grid
set key top right
plot './C3/distance_averaged.csv' using 1:9 smooth sbezier w lines lw 2 lc rgb 'red' title 'C3', \
     './MMTV/distance_averaged.csv' using 1:9 smooth sbezier w lines lw 2 lc rgb 'orange' title 'MMTV', \
     './DR3/distance_averaged.csv' using 1:9 smooth sbezier w lines lw 2 lc rgb 'yellow' title 'DR3'
quit
EOL

gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set output 'MMTV NON_CAN_A568R_D15G.png'
set title 'MMTV NON_CAN_A568R_D15G'
set ylabel 'Distance (nm)'
set xlabel 'Time (ps)'
set grid
set key top right
plot './2_first_run/dist.xvg' using 1:13 smooth sbezier w lines lw 2 lc rgb 'red' title 'First run', \
     './2_second_run/dist.xvg' using 1:13 smooth sbezier w lines lw 2 lc rgb 'orange' title 'Second run', \
     './2_third_run/dist.xvg' using 1:13 smooth sbezier w lines lw 2 lc rgb 'yellow' title 'Third run'
quit
EOL

#violinplot
seaborn script in jupyter notebook

#STEP 16: SASA
gmx make_ndx -f md.gro -o index_sasa.ndx -quiet << EOL
1|13
1|12|13
q
EOL
gmx sasa -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -s md.pdb -n index_sasa.ndx -quiet -o sasa_complex.xvg << EOL
19
EOL
gmx sasa -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -s md.pdb -n index_sasa.ndx -quiet -o sasa_protein.xvg << EOL
18
EOL
gmx sasa -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -s md.pdb -n index_sasa.ndx -quiet -o sasa_dna.xvg << EOL
12
EOL

#MASS RENAMING: for f in *.png; do mv "$f" "C3_first_$f"; done

#STEP 17: delete frames when structure is unstable using RMSD

#STEP 18: save as .dcd file and .pdb file in VMD

#STEP 19: continue analysis using MD_script_PCA.sh

#EXTRA LEVER ARM CALCULATIONS:

gmx make_ndx -f md.gro -o index_la.ndx -quiet << EOL
ri 25-42 & 4
ri 96-113 & 4
name 18 non_canonical_la
name 19 canonical_la
quit
EOL

gmx rms -s md.tpr -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -n index_la.ndx -o nc_la.xvg << EOL
18
18
EOL
gmx rms -s md.tpr -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -n index_la.ndx -o can_la.xvg << EOL
19
19
EOL

gmx trjconv -s md.tpr -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -o md_nc_la.xtc -fit rot+trans -quiet -n index_la.ndx << EOL
Backbone
18
EOL
gmx trjconv -s md.tpr -f md.xtc -o md_nc_la.pdb -pbc mol -ur compact -dump 0 -quiet -n index_la.ndx << EOL
18
EOL
gmx trjconv -s md.tpr -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -o md_can_la.xtc -fit rot+trans -quiet -n index_la.ndx << EOL
Backbone
19
EOL
gmx trjconv -s md.tpr -f md.xtc -o md_can_la.pdb -pbc mol -ur compact -dump 0 -quiet -n index_la.ndx << EOL
19
EOL

gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set output 'C3_nc_la_rmsd.png'
set title 'C3 non canonical lever arm RMSD'
set ylabel 'RMSD (nm)'
set xlabel 'Time (ps)'
set grid
set key top right
plot './2_first_run/nc_la.xvg' using 1:2 smooth sbezier w lines lw 2 lc rgb 'red' title 'First run', \
     './2_second_run/nc_la.xvg' using 1:2 smooth sbezier w lines lw 2 lc rgb 'orange' title 'Second run', \
     './2_third_run/nc_la.xvg' using 1:2 smooth sbezier w lines lw 2 lc rgb 'yellow' title 'Third run'
quit
EOL

gmx rms -s '/home/wout/0_AR/7_ARE/5_isolated_leverarm/C3_can_LA.pdb'  -f '/home/wout/0_AR/7_ARE/1_molecular_dynamics/C3/2_first_run/md_la.xtc' -fit rot+trans

#GBSA STABLE DNA SELECTION:

gmx make_ndx -f md.tpr -o index_mmpbsa << EOL
1 | 12 | 13
14 | 15| 16
1 | 13
ri 144-164
ri 167-187
21 | 22
name 23 stable_DNA
20 | 23
name 24 Protein_stable_DNA_ZN
ri 147-161
ri 170-184
25 | 26
name 27 Truncated_DNA
20 | 27
name 28 Protein_truncated_DNA_ZN
q
EOL

