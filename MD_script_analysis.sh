#Recentering and rewrapping coordinates of the trajectory files

#STEP 1: make index file
gmx make_ndx -f md.gro -o index.ndx -quiet << EOL
1|12|13
q
EOL

#STEP 2: makes broken protein whole
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC_whole.xtc -pbc whole -quiet -n index.ndx << EOL
18
EOL

#STEP 3: eliminates jumps from one side to the other side of the box
gmx trjconv -s md.tpr -f  md_noPBC_whole.xtc -pbc nojump -o md_noPBC_whole_nojump.xtc -quiet -n index.ndx << EOL
18
EOL

#STEP 4: center protein in box
gmx trjconv -s md.tpr -f md_noPBC_whole_nojump.xtc -o md_noPBC_whole_nojump_center.xtc -center -quiet -n index.ndx << EOL
Protein
18
EOL

#STEP 5: center all molecules in box and put all atoms at the closest distance from center of the box
gmx trjconv -s md.tpr -f md_noPBC_whole_nojump_center.xtc -o md_noPBC_whole_nojump_center_mol_com.xtc -pbc mol -ur compact -quiet -n index.ndx << EOL
18
EOL

#STEP 6: fit the system to reference in structure file
gmx trjconv -s md.tpr -f md_noPBC_whole_nojump_center_mol_com.xtc -o md_noPBC_whole_nojump_center_mol_com_fit.xtc -fit rot+trans -quiet -n index.ndx << EOL
Backbone
18
EOL
rm md_noPBC_whole.xtc md_noPBC_whole_nojump.xtc md_noPBC_whole_nojump_center.xtc md_noPBC_whole_nojump_center_mol_com.xtc

#STEP 6: output a .pdb file that contains the exact particles as the trajectory at t=0, so other programs know how to interpret the numbers in the trajectory
gmx trjconv -s md.tpr -f md.xtc -o md.pdb -pbc mol -ur compact -dump 0 -quiet -n index.ndx << EOL
18
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
gmx rms -s md.tpr -f md_noPBC_whole_nojump_center_mol_com.xtc -o rmsd.xvg -tu ns -quiet << EOL
Backbone
Backbone
EOL

#STEP 11: visualise structural stability
gnuplot << EOL 
set terminal png size 1000, 800 enhanced font "Helvetica,14"
set output 'rmsd.png'
set title 'RMSD, backbone'
set ylabel 'RMSD (nm)'
set xlabel 'Time (ns)'
#set xrang [0:2000]
unset key  
set grid
plot 'rmsd.xvg' u 1:2 w lines lw 2 lc rgb '#27ad81'
quit
EOL
okular rmsd.png

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
gmx rmsf -oq B_factor.pdb -ox B_factor_average_trajectory.pdb -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -s md.tpr -n index.ndx << EOL
20
EOL
#Adjust in MOE if file messed up
show cartoon
cartoon putty
unset cartoon_smooth_loops
unset cartoon_flat_sheets
set ray_shadows,0
set ray_trace_mode,1
spectrum b, rainbow, minimum=0, maximum=100
ramp_new color_bar, B_factor, [0, 100], rainbow
# for i in range(1, 80): cmd.bond(f"/B_factor//C/{i}/O3'", f"/B_factor//C/{i+1}/P")
# set label_color, black, B_factor

#STEP 14: delete frames when structure is unstable using RMSD

#STEP 15: save as .dcd file and .pdb file in VMD

#STEP 16: continue analysis using MD_script_PCA.sh
