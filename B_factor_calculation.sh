for dir in */; do
    cd $dir
    echo "$dir"   
    for run_dir in 2_run 2_run_second 2_run_third; do
	cd "$run_dir"
	echo "$run_dir"
/opt/gromacs/gromacs_2023/bin/gmx rmsf -oq B_factor.pdb -ox B_factor_average_trajectory.pdb -f md_mmpbsa.xtc -s md.tpr -n index.ndx << EOL
21
EOL
	rm B_factor_average_trajectory.pdb
	rm rmsf.xvg
	dirname=$(basename `pwd`)
	awk '/^ATOM/ && ($3=="CA" || $3=="C4\x27") {chain=substr($0,22,1); resnum=substr($0,23,4); print chain resnum}' B_factor.pdb > extracted_chain_and_residue_$dirname.txt
	awk '/^ATOM/ && ($3=="CA" || $3=="C4\x27") {print substr($0,61,7)}' B_factor.pdb > extracted_B_factor_$dirname.txt
	cp extracted_chain_*.txt ..
	cp extracted_B_factor_*.txt ..
	cp B_factor.pdb ..
	rm -f \#*\#
	cd ..
    done
    paste extracted_B_factor*.txt > merged_B_factor.txt
    mv extracted_chain_and_residue_2_run.txt extracted_chain_and_residue.txt
    #rm extracted_chain_and_residue_2_run_second.txt extracted_chain_and_residue_2_run_third.txt
    awk '{printf "%s %.2f\n", $0, ($1+$2+$3)/3}' merged_B_factor.txt > averaged_B_factor.txt
    awk '{print $4}' averaged_B_factor.txt > averaged_B_factor2.txt
    sed -i 's/\,/./g' averaged_B_factor2.txt
    paste extracted_chain_and_residue.txt averaged_B_factor2.txt > input_insertion_py.txt
    vim -c "%s/NGLYA/NGL A/g | %s/NGLYB/NGL B/g | %s/CTHRA/CTH A/g | %s/CTHRB/CTH B/g" -c "wq" B_factor.pdb
    python /home/wout/Scripts/1_molecular_dynamics/B_factor/insert-bfactor.py -p B_factor.pdb -b input_insertion_py.txt -o B_factor_averaged.pdb -d 0.0
    vim -c "%s/NGL A/NGLYA/g | %s/NGL B/NGLYB/g | %s/CTH A/CTHRA/g | %s/CTH B/CTHRB/g" -c "wq" B_factor_averaged.pdb
    dirname2=$(basename `pwd`)
    mv B_factor_averaged.pdb B_factor_averaged_$dirname2.pdb
    rm extracted*.txt input_insertion_py.txt merged_B_factor.txt B_factor.pdb averaged_B_factor.txt averaged_B_factor2.txt 
    cd ..
done
