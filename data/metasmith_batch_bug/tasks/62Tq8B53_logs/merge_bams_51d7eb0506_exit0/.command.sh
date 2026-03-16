#!/bin/bash -ue
echo "step 2, sample [[Rsoo8XE7:[408064623690479551], k0BWqyxf:[310195413067453707, 811104384982471761, 305980042534366556, 1073592773866965764, 893397822587109542, 822627013523906702], AmC7DKKQ:[506249955841459167], FILES:[[porphyridium_experiment.txt], [1-1-1.Grq6HkY9LZzPhIwC-9mrjFffM.bam, 1-1-1.UMK0uBW7MDfv9lR2-9mrjFffM.bam, 1-1-1.L8rJHsIV6bp7TLk8-9mrjFffM.bam, 1-1-1.4Ih5b6oUlLPBmjQj-9mrjFffM.bam, 1-1-1.ZjUpCov8SiawjbWC-9mrjFffM.bam, 1-1-1.SZV4oNiEOp2dIAHQ-9mrjFffM.bam], [samtools.oci]]]]"
echo "merge_bams"
echo "res 8/16 GB/1" >>.command.metadata
echo "lin [{\"Rsoo8XE7\":[408064623690479551],\"k0BWqyxf\":[310195413067453707,811104384982471761,305980042534366556,1073592773866965764,893397822587109542,822627013523906702],\"AmC7DKKQ\":[506249955841459167],\"FILES\":[[\"porphyridium_experiment.txt\"],[\"1-1-1.Grq6HkY9LZzPhIwC-9mrjFffM.bam\",\"1-1-1.UMK0uBW7MDfv9lR2-9mrjFffM.bam\",\"1-1-1.L8rJHsIV6bp7TLk8-9mrjFffM.bam\",\"1-1-1.4Ih5b6oUlLPBmjQj-9mrjFffM.bam\",\"1-1-1.ZjUpCov8SiawjbWC-9mrjFffM.bam\",\"1-1-1.SZV4oNiEOp2dIAHQ-9mrjFffM.bam\"],[\"samtools.oci\"]]}]" >>.command.metadata
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
