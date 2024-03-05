# Script to automate gmx

for dir in */; do
    cd $dir
    echo "$dir"
    cp -r 1_startfiles 2_run
    cp -r 1_startfiles 2_run_third
    cp -r 1_startfiles 2_run_second    
    for run_dir in 2_run 2_run_second 2_run_third; do
	cd "$run_dir"
	echo "$run_dir"
rm step4.0_minimization.mdp step4.1_equilibration.mdp step5_production.mdp step3_input.pdb step3_input.psf README index.ndx
mv step3_input.gro protein_solv_ions.gro
/usr/local/gromacs/bin/gmx grompp -f minim.mdp -c protein_solv_ions.gro -r protein_solv_ions.gro -p topol.top -o em.tpr -quiet -maxwarn 1
/usr/local/gromacs/bin/gmx mdrun -v -deffnm em -quiet -nb gpu
/usr/local/gromacs/bin/gmx energy -f em.edr -o potential.xvg -quiet << EOL
Potential
EOL
/usr/local/gromacs/bin/gmx make_ndx -f em.gro -o index.ndx -quiet << EOL
1|13|14|15|16
17|18|19
quit
EOL
/usr/local/gromacs/bin/gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx -quiet -maxwarn 1
/usr/local/gromacs/bin/gmx mdrun -v -deffnm nvt -quiet -nb gpu
/usr/local/gromacs/bin/gmx energy -f nvt.edr -o temperature.xvg -quiet << EOL
Temperature
EOL
/usr/local/gromacs/bin/gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -n index.ndx -quiet -maxwarn 1
/usr/local/gromacs/bin/gmx mdrun -v -deffnm npt -quiet -nb gpu
/usr/local/gromacs/bin/gmx energy -f npt.edr -o pressure.xvg -quiet << EOL
Pressure
EOL
/usr/local/gromacs/bin/gmx energy -f npt.edr -o density.xvg -quiet << EOL
Density
EOL
/usr/local/gromacs/bin/gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr -n index.ndx -quiet -maxwarn 1
/usr/local/gromacs/bin/gmx mdrun -v -deffnm md -quiet -gputasks 0001 -nb gpu -pme gpu -npme 1 -ntmpi 4 &
wait $!
    	cd ..
    done
    cd ..
done

echo "The automated calculation has finished"
