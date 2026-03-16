#!/bin/bash -ue
echo "step 3, sample [[Rsoo8XE7:[408064623690479551], VwXt5dAr:[178111524324759755], k0BWqyxf:[156226118216467402, 633523990048773908, 566498055328118649, 310195413067453707, 811104384982471761, 305980042534366556, 1073592773866965764, 893397822587109542, 822627013523906702], AmC7DKKQ:[506249955841459167], iNlpm1XR:[955161200972376201, 519030939364663020], ea9yJktW:[1072381707910174208], FILES:[[porphyridium_experiment.txt], [1-1-1.f1CMorcneUoLGMna-O4PhHAkd.fna], [1-1-1.hqXlyWE1DV3e7nDs-iNlpm1XR.bam, 1-1-1.PfamljgePG5VX967-iNlpm1XR.bam], [braker3.oci]]]]"
echo "braker3"
echo "res 16/64 GB/1" >>.command.metadata
echo "lin [{\"Rsoo8XE7\":[408064623690479551],\"VwXt5dAr\":[178111524324759755],\"k0BWqyxf\":[156226118216467402,633523990048773908,566498055328118649,310195413067453707,811104384982471761,305980042534366556,1073592773866965764,893397822587109542,822627013523906702],\"AmC7DKKQ\":[506249955841459167],\"iNlpm1XR\":[955161200972376201,519030939364663020],\"ea9yJktW\":[1072381707910174208],\"FILES\":[[\"porphyridium_experiment.txt\"],[\"1-1-1.f1CMorcneUoLGMna-O4PhHAkd.fna\"],[\"1-1-1.hqXlyWE1DV3e7nDs-iNlpm1XR.bam\",\"1-1-1.PfamljgePG5VX967-iNlpm1XR.bam\"],[\"braker3.oci\"]]}]" >>.command.metadata
echo "inp Rsoo8XE7,VwXt5dAr,iNlpm1XR,ea9yJktW" >>.command.metadata
echo "out 3J2d8bJE,N3PPDtPP" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/eguEpdhP-intermediates/assembly"
echo "--bind $b1:$b1" >.command.binds

CONTAINER=/scratch/st-shallam-1/pwy_group/metasmith
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

bootstrap /scratch/st-shallam-1/pwy_group/metasmith/runs/62Tq8B53 "3" login01
[ -e .command.success ] && exit 0 || exit 1
