#!/bin/bash
# run in folder like N=0.7dumps

for fol in */; do
	cd $fol
	for i in *.dump; do
		pref=$(echo "$i" | sed 's/.dump//g')
		echo 'starting on' $i
		python /home/kamron/mdext_KF/examples/RDF_1D.py ${i} ${pref} 10 
		# python /home/kamron/mdext_KF/examples/RDF_1D.py 1D_T1.6_P1.0.dump test 10 
		# xx is a placeholder in case triclini output of positions are used 
		# exit 1
	done
	cd ..
done
