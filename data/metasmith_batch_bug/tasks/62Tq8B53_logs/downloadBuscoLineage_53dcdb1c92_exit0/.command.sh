#!/bin/bash -ue
echo "step 1, sample [[H7g3InfC:[54446232961680702], OvFEfok5:[634031039128142113], FILES:[[busco.oci], [busco_source.txt]]]]"
echo "downloadBuscoLineage"
echo "res 1/4 GB/1" >>.command.metadata
echo "lin [{\"H7g3InfC\":[54446232961680702],\"OvFEfok5\":[634031039128142113],\"FILES\":[[\"busco.oci\"],[\"busco_source.txt\"]]}]" >>.command.metadata
echo "inp H7g3InfC,OvFEfok5" >>.command.metadata
echo "out b83cdTFq" >>.command.metadata
echo "" >.command.binds

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

bootstrap /scratch/st-shallam-1/pwy_group/metasmith/runs/62Tq8B53 "1" login01
[ -e .command.success ] && exit 0 || exit 1
