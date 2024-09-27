# Script to automate gmx

for dir in */; do
    cd $dir
    echo "$dir"
    for run_dir in 2_run 2_run_second 2_run_third; do
	cd "$run_dir"
	echo "$run_dir"
	
cd toppar
cp PROA.itp PROA_adjusted.itp
cp PROB.itp PROB_adjusted.itp
vim PROA_adjusted.itp << EOL
:%s,0YB,0YA,g
:%s,1      0YA,2      0YB,gc
nnnnnnnnnnnnnnnnnnnnnnnnnnnnyyyyyyyyyyyyyyyyyyyyyyyyyyyl
:%s,1      0YA,3      0YC,gc
nnnnnnnnnnnnnnnnnnnnnnnnnnnnyyyyyyyyyyyyyyyyyyyyyyyyyyyl
:wq
EOL
vim PROB_adjusted.itp << EOL
:%s,0YB,0YA,g
:%s,1      0YA,1      0YD,gc
yyyyyyyyyyyyyyyyyyyyyyyyyyyl
:%s,1      0YA,2      0YE,gc
yyyyyyyyyyyyyyyyyyyyyyyyyyyl
:%s,1      0YA,3      0YF,gc
yyyyyyyyyyyyyyyyyyyyyyyyyyyl
:wq
EOL
..

cp topol.top topol_adjusted.top
vim topol_adjusted.top << EOL
:%s,toppar/PROA.itp,toppar/PROA_adjusted.itp,g
:%s,toppar/PROB.itp,toppar/PROB_adjusted.itp,g
:wq
EOL

cp md_mmpbsa.pdb md_adjusted.pdb

vim md_adjusted.pdb << EOL
:%s,0YB,0YA,g
:%s,0YA A   1,0YB A   2,gc
nnnnnnnnnnnnnnnnnnnnnnnnnnnnyyyyyyyyyyyyyyyyyyyyyyyyyyyl
:%s,0YA A   1,0YC A   3,gc
nnnnnnnnnnnnnnnnnnnnnnnnnnnnyyyyyyyyyyyyyyyyyyyyyyyyyyyl
:%s,0YA B   1,0YD B   1,gc
yyyyyyyyyyyyyyyyyyyyyyyyyyyl
:%s,0YA B   1,0YE B   2,gc
yyyyyyyyyyyyyyyyyyyyyyyyyyyl
:%s,0YA B   1,0YF B   3,gc
yyyyyyyyyyyyyyyyyyyyyyyyyyyl
:wq
EOL

rm -f \#*\#
    	cd ..
    done
    cd ..
done

echo "The automated calculation has finished"
