# Script to automate MMPBSA calculation

mkdir MMPBSA_results

for dir in */; do
    echo "$dir"
    cd $dir
    mpirun -np 120 gmx_MMPBSA MPI -i mmpbsa.in -cs md.tpr -ci index2.ndx -cg 1 12 -ct md_noPBC_whole_nojump_center_mol_com_fit.xtc -cp topol.top -O &
    wait $!
    echo "The calculation has finished."
    dirname=$(basename `pwd`)
    cp FINAL_RESULTS_MMPBSA.dat MMPBSA_$dirname.dat
    cp MMPBSA_*.dat ../MMPBSA_results
    cd ..
done

echo "The automated calculation has finished"
