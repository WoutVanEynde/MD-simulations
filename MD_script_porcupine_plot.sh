#STEP 1: build and diagonalize the covariance matrix
gmx covar -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -s md.tpr << EOL
C-alpha
C-alpha
EOL

#STEP 2: project an MD trajectory along a 1st eigenvector
gmx anaeig -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -v eigenvec.trr -eig eigenval.xvg -s md.tpr -proj proj_1_2.xvg -b 0 -e 200000 -first 1 -last 1 << EOL
C-alpha
C-alpha
EOL

#STEP 3: extract time scale from proj_1_2.xvg
awk '{print $1}' proj_1_2.xvg > timescale.csv
vim timescale.csv
:set number
:1,24d

#STEP 4: project an MD trajectory in two dimensions of selected eigenvectors (essential subspace)
gmx anaeig -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -v eigenvec.trr -eig eigenval.xvg -s md.tpr -2d 2dproj_1_2.xvg  -b 0 -e 200000 -first 1 -last 2 <<EOL
C-alpha
C-alpha
EOL
vim 2dproj_1_2.xvg
set number
:1,17d

#STEP 5: merge timescale into 2dproj_1_2.xvg
paste 2dproj_1_2.xvg timescale.csv  > 2dproj_timescale.csv
vim 2dproj_timescale.csv

#STEP 6: plot PCA 2D projection
gnuplot << EOL
set terminal png size 600, 400
set output '2D_projection_1_2.png'
set title '2D-projection' 
set ylabel 'projection on eigenvector 2 (nm)'
set xlabel 'projection on eigenvector 1 (nm)' 
#set xrang [-8:6]
#set yrang [-5:4]
set grid
unset key
set palette defined ( 0 "#000090",\
        1 "#000fff",\
        2 "#0090ff",\
        3 "#0fffee",\
        4 "#90ff70",\
        5 "#ffee00",\
        6 "#ff7000",\
        7 "#ee0000",\
        8 "#7f0000")
plot '2dproj_timescale.csv' u 1:2:3 with points palette
quit
EOL
okular 2D_projection_1_2.png 

#STEP 7: porcupine
gmx anaeig -v eigenvec.trr -s md.tpr -f md_noPBC_whole_nojump_center_mol_com_fit.xtc -b 0 -e 200000 -first 1 -last 1 -extr extreme1.gro -nframes 2000 <<EOL
C-alpha
C-alpha
EOL

#STEP 8: visualise using modevectors.py
sele resn CLEU
remove sele
sele resn NGLN
remove sele
split_states extreme1, 1, 2
split_states extreme1, 1999, 2000
set_name extreme1_0001, 0001
set_name extreme1_2000, 2000
delete extreme1_0002
delete extreme1_1999
run modevectors.py
modevectors 0001, 2000, cutoff=0, head_length=1, head=0.4, headrgb=(1,.2,.1), tailrgb=(1,.2,.1), notail=0
bg_color black
# Dont do preset.publication(selection='all'), will alter modevectors
set cartoon_ring_mode, 3
set ray_shadows,0
set ray_trace_mode,1
set ray_trace_color, black
bg_color black

#STEP 9: clean backupfiles
rm -f \#*\#
