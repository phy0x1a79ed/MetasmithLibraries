#!/bin/bash -ue
array=( /scratch/st-shallam-1/pwy_group/metasmith/runs/62Tq8B53/nxf_work/51/d7eb050607a129bab15c5b17002a6f /scratch/st-shallam-1/pwy_group/metasmith/runs/62Tq8B53/nxf_work/91/34572e74bf8f1531dbbec6d1817979 )
export nxf_array_task_dir=${array[SLURM_ARRAY_TASK_ID]}
bash $nxf_array_task_dir/.command.run 2>&1 > $nxf_array_task_dir/.command.log
