#!/bin/bash
# run in data folder

module purge
module load venv/mdextKF

# for PTfol in */; do
# 	cd $PTfol
#     PTfolName=$(echo "$PTfol" | sed 's/\///g')  # remove the / from the folder name
    
for seedfol in */; do
# for seedfol in seed0000; do
    cd $seedfol
    seedfolName=$(echo "$seedfol" | sed 's/\///g')  # remove the / from the folder name
    
    for replicafol in */; do
        # for seedfol in seed0000; do
        cd $replicafol
        replicafolName=$(echo "$replicafol" | sed 's/\///g')  # remove the / from the folder name
        
        # build H5
        # h5Name=${seedfolName}${replicafolName}Grouped.h5  
        h5Name=${seedfolName}Grouped.h5 # HACK - no replica
        echo $h5Name
        python /home/kamron/mdext_KF/examples/BuildH5forTraining.py ${h5Name}
        python /home/kamron/mdext_KF/examples/plotDensityTI.py ${h5Name}
        cp $h5Name ../.
        cd ..
    done

    # plot all replicas after averaging
    python /home/kamron/mdext_KF/examples/plotDensityTIReplicas.py
    cd ..
done
# cd ..


# done
