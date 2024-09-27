for dir in */; do
    cd $dir
    echo "$dir"   
    for run_dir in 2_run 2_run_second 2_run_third; do
	cd "$run_dir"
	echo "$run_dir"
/opt/gromacs/gromacs_2023/bin/gmx rmsf -oq B_factor.pdb -ox B_factor_average_trajectory.pdb -f md_mmpbsa.xtc -s md.tpr -n index.ndx << EOL
18
EOL
	rm B_factor_average_trajectory.pdb
	rm rmsf.xvg
	dirname=$(basename `pwd`)
	awk '/^ATOM/ {chain=substr($0,22,1); atomnum=substr($0,8,4); print chain atomnum}' B_factor.pdb > extracted_chain_and_atom_full_$dirname.txt
	awk '/^ATOM/ {print substr($0,61,7)}' B_factor.pdb > extracted_B_factor_full_$dirname.txt
	cp extracted_chain_*.txt ..
	cp extracted_B_factor_*.txt ..
	cp B_factor.pdb ..
	rm -f \#*\#
	cd ..
    done
    paste extracted_B_factor_full_*.txt > merged_B_factor_full.txt
    mv extracted_chain_and_atom_full_2_run.txt extracted_chain_and_atom_full.txt
    rm extracted_chain_and_atom_full_2_run_second.txt extracted_chain_and_atom_full_2_run_third.txt
    awk '{printf "%s %.2f\n", $0, ($1+$2+$3)/3}' merged_B_factor_full.txt > averaged_B_factor_full.txt
    awk '{print $4}' averaged_B_factor_full.txt > averaged_B_factor2_full.txt
    sed -i 's/\,/./g' averaged_B_factor2_full.txt
    paste extracted_chain_and_atom_full.txt averaged_B_factor2_full.txt > input_insertion_py_full.txt
    vim -c "%s/A/A /g | %s/C/C /g | %s/B/B /g" -c "wq" input_insertion_py_full.txt
    #vim -c ":%s,  ,A  ,g" -c "wq" input_insertion_py.txt
    #vim -c ":%s,   ,  ,g" -c "wq" input_insertion_py.txt
    #python3 /home/wout/Scripts/1_molecular_dynamics/B_factor/chain_replacer.py #replace all chain tags to A
    cp B_factor.pdb modified_B_factor_full.pdb
    awk '/^ATOM/ { $0 = substr($0, 1, 60) " " substr($0, 61) }1' modified_B_factor_full.pdb > modified_B_factor_full2.pdb
    vim -c "%s/NILEA/NIL A/g | %s/CGLNA/CGL A/g" -c "wq" modified_B_factor_full2.pdb
    python3 /home/wout/Scripts/1_molecular_dynamics/B_factor/full_replacements_B_factor.py
    vim -c "%s/NIL A/NILEA/g | %s/CGL A/CGLNA/g" -c "wq" B_factor_averaged_full.pdb
    dirname2=$(basename `pwd`)
    mv B_factor_averaged_full.pdb B_factor_averaged_full_$dirname2.pdb
    rm extracted*.txt merged_B_factor_full.txt B_factor.pdb averaged_B_factor_full.txt averaged_B_factor2_full.txt input_insertion_py_full.txt modified_B_factor_full.pdb modified_B_factor_full2.pdb
    cd ..
done
