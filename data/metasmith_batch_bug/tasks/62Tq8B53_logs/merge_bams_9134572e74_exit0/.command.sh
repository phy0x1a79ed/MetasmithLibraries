#!/bin/bash -ue
echo "step 2, sample [[Rsoo8XE7:[408064623690479551], k0BWqyxf:[156226118216467402, 633523990048773908, 566498055328118649], AmC7DKKQ:[506249955841459167], FILES:[[porphyridium_experiment.txt], [1-1-1.g3ah0QAiGjmgmOQv-9mrjFffM.bam, 1-1-1.iCaeK0EUwqlMC7ZR-9mrjFffM.bam, 1-1-1.SYqoyBxlvuPEILkU-9mrjFffM.bam], [samtools.oci]]]]"
echo "merge_bams"
echo "res 8/16 GB/1" >>.command.metadata
echo "lin [{\"Rsoo8XE7\":[408064623690479551],\"k0BWqyxf\":[156226118216467402,633523990048773908,566498055328118649],\"AmC7DKKQ\":[506249955841459167],\"FILES\":[[\"porphyridium_experiment.txt\"],[\"1-1-1.g3ah0QAiGjmgmOQv-9mrjFffM.bam\",\"1-1-1.iCaeK0EUwqlMC7ZR-9mrjFffM.bam\",\"1-1-1.SYqoyBxlvuPEILkU-9mrjFffM.bam\"],[\"samtools.oci\"]]}]" >>.command.metadata
echo "inp Rsoo8XE7,k0BWqyxf,AmC7DKKQ" >>.command.metadata
echo "out iNlpm1XR" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/star_bams"
echo "--bind $b1:$b1" >.command.binds

CONTAINER=/msm_home
DIRECT=/scratch/st-shallam-1/pwy_group/metasmith
function bootstrap {
	if [ -e $CONTAINER ]; then
		$CONTAINER/lib/msm_bootstrap $@
	elif [ -e $DIRECT ]; then
		$DIRECT/lib/msm_bootstrap $@
	else
		echo "critical error: could not find metasmith bootstrap script"
	fi
}

bootstrap /scratch/st-shallam-1/pwy_group/metasmith/runs/62Tq8B53 "2" login01
[ -e .command.success ] && exit 0 || exit 1
