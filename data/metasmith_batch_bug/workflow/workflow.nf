params.testSpread=1
params.home = '/scratch/st-shallam-1/pwy_group/metasmith'
params.workspace = "${params.home}/runs/62Tq8B53"
params.bootstrap_def = '''
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
'''

import groovy.json.JsonSlurper
def in(f, l) {
    def rows = Channel.fromPath(f).splitCsv(header: false)
    if (f in l) {
        rows = Channel.fromList(l[f]).merge(rows)
    }
    return rows.map((row) -> {
        if (row.size()>1) {
            def (ri, rx) = row
            return tuple(ri, file(rx))
        } else {
            def i = [:]
            return tuple(i, file(row[0]))
        }
    })
}


process p01__downloadBuscoLineage {
	label 'xXzlTFjWrx'
	label 'xlocalx'
input:
	tuple val(index),path(_01),path(_02)
output:
	tuple val(index),path("*-1.*-b83cdTFq")
script:
"""
echo "step 1, sample $index"
echo "downloadBuscoLineage"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp H7g3InfC,OvFEfok5" >>.command.metadata
echo "out b83cdTFq" >>.command.metadata
echo "" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "1" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-b83cdTFq"
"""
}

process p02__merge_bams {
	label 'xh2y89BxAx'
input:
	tuple val(index),path(_01),path(_02),path(_03)
output:
	tuple val(index),path("*-1.*-iNlpm1XR.bam")
script:
"""
echo "step 2, sample $index"
echo "merge_bams"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp Rsoo8XE7,k0BWqyxf,AmC7DKKQ" >>.command.metadata
echo "out iNlpm1XR" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/star_bams"
echo "--bind \$b1:\$b1" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "2" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-iNlpm1XR.bam"
"""
}

process p03__braker3 {
	label 'x5c7Jq5xkx'
input:
	tuple val(index),path(_01),path(_02),path(_03),path(_04)
output:
	tuple val(index),path("*-1.*-3J2d8bJE.gff3")
	tuple val(index),path("*-1.*-N3PPDtPP.faa")
script:
"""
echo "step 3, sample $index"
echo "braker3"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp Rsoo8XE7,VwXt5dAr,iNlpm1XR,ea9yJktW" >>.command.metadata
echo "out 3J2d8bJE,N3PPDtPP" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/eguEpdhP-intermediates/assembly"
echo "--bind \$b1:\$b1" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "3" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-3J2d8bJE.gff3" "1-1-1.test$hash-N3PPDtPP.faa"
"""
}

process p04__stringtie_merge {
	label 'xHQOhnhjTx'
input:
	tuple val(index),path(_01),path(_02),path(_03),path(_04)
output:
	tuple val(index),path("*-1.*-Myb7Pxjl.gtf")
script:
"""
echo "step 4, sample $index"
echo "stringtie_merge"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp Rsoo8XE7,X81sU15q,3J2d8bJE,rRC4FxHn" >>.command.metadata
echo "out Myb7Pxjl" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/eguEpdhP-intermediates/stringtie_gtfs"
echo "--bind \$b1:\$b1" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "4" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-Myb7Pxjl.gtf"
"""
}

process p05__gffread_proteins {
	label 'xYObZ0Pxzx'
input:
	tuple val(index),path(_01),path(_02),path(_03),path(_04)
output:
	tuple val(index),path("*-1.*-LmS8pDph.faa")
script:
"""
echo "step 5, sample $index"
echo "gffread_proteins"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp Rsoo8XE7,3J2d8bJE,VwXt5dAr,56zi3dIJ" >>.command.metadata
echo "out LmS8pDph" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/eguEpdhP-intermediates/assembly"
echo "--bind \$b1:\$b1" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "5" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-LmS8pDph.faa"
"""
}

process p06__stringtie_quant {
	label 'xs8VnlfHUx'
input:
	tuple val(index),path(_01),path(_02),path(_03)
output:
	tuple val(index),path("*-1.*-pDpp0Ttl.gtf")
script:
"""
echo "step 6, sample $index"
echo "stringtie_quant"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp k0BWqyxf,Myb7Pxjl,rRC4FxHn" >>.command.metadata
echo "out pDpp0Ttl" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/star_bams"
echo "--bind \$b1:\$b1" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "6" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-pDpp0Ttl.gtf"
"""
}

process p07__busco {
	label 'xK8Rict2dx'
input:
	tuple val(index),path(_01),path(_02),path(_03)
output:
	tuple val(index),path("*-1.*-yTyedsZ0.json")
script:
"""
echo "step 7, sample $index"
echo "busco"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp H7g3InfC,LmS8pDph,b83cdTFq" >>.command.metadata
echo "out yTyedsZ0" >>.command.metadata
echo "" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "7" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-yTyedsZ0.json"
"""
}

process p08__eggnog_mapper {
	label 'xFd8XQrQPx'
input:
	tuple val(index),path(_01),path(_02),path(_03)
output:
	tuple val(index),path("*-1.*-F0MqIpV6.tsv")
script:
"""
echo "step 8, sample $index"
echo "eggnog_mapper"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp 0uxln7h5,LmS8pDph,MUbLkijl" >>.command.metadata
echo "out F0MqIpV6" >>.command.metadata
b1="/arc/project/st-shallam-1/pwy_group/lib/annotation-dbs"
echo "--bind \$b1:\$b1" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "8" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-F0MqIpV6.tsv"
"""
}

process p09__pydeseq2 {
	label 'xEb5DENOEx'
input:
	tuple val(index),path(_01),path(_02),path(_03)
output:
	tuple val(index),path("*-1.*-OJTSKg7W.tsv")
script:
"""
echo "step 9, sample $index"
echo "pydeseq2"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp Rsoo8XE7,pDpp0Ttl,Z9xw077K" >>.command.metadata
echo "out OJTSKg7W" >>.command.metadata
echo "" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "9" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-OJTSKg7W.tsv"
"""
}

process p10__stringtie_count_matrix {
	label 'xOpDNsOczx'
input:
	tuple val(index),path(_01),path(_02),path(_03)
output:
	tuple val(index),path("*-1.*-5TR9TDph.csv")
script:
"""
echo "step 10, sample $index"
echo "stringtie_count_matrix"
echo "res $task.cpus/$task.memory/$task.attempt" >>.command.metadata
echo "lin ${Orchestrator.JsonforEcho(index)}" >>.command.metadata
echo "inp Rsoo8XE7,pDpp0Ttl,dB3q1YB1" >>.command.metadata
echo "out 5TR9TDph" >>.command.metadata
echo "" >.command.binds
${params.bootstrap_def}
bootstrap ${params.workspace} "10" ${params.hostName}
[ -e .command.success ] && exit 0 || exit 1
"""
stub:
def dt = new Random().nextFloat()*params.testSpread
def hash = "${index[0].sort().collectEntries((k, v) -> [k, v.sort()])}".md5()[0..11]
"""
sleep $dt
touch "1-1-1.test$hash-5TR9TDph.csv"
"""
}

workflow {
main:
o = new Orchestrator(Channel.fromList([null])) // cant create channels in groovy
l = new JsonSlurper().parseText(file("workflow.lineage_of_given.json").text)
o.child2parent["VwXt5dAr"] = (["Rsoo8XE7"] as Set)
o.child2parent["k0BWqyxf"] = (["Rsoo8XE7"] as Set)
o.child2parent["X81sU15q"] = (["Rsoo8XE7"] as Set)
(_OvFEfok5) = o.postIn([in("inputs/OvFEfok5", l)], ["OvFEfok5"]) // annotation::busco_source
(_MUbLkijl) = o.postIn([in("inputs/MUbLkijl", l)], ["MUbLkijl"]) // annotation::eggnog_data
(_ea9yJktW) = o.postIn([in("inputs/ea9yJktW", l)], ["ea9yJktW"]) // containers::braker3.oci
(_H7g3InfC) = o.postIn([in("inputs/H7g3InfC", l)], ["H7g3InfC"]) // containers::busco.oci
(_0uxln7h5) = o.postIn([in("inputs/0uxln7h5", l)], ["0uxln7h5"]) // containers::eggnog-mapper.oci
(_56zi3dIJ) = o.postIn([in("inputs/56zi3dIJ", l)], ["56zi3dIJ"]) // containers::gffread.oci
(_Z9xw077K) = o.postIn([in("inputs/Z9xw077K", l)], ["Z9xw077K"]) // containers::pydeseq2.oci
(_dB3q1YB1) = o.postIn([in("inputs/dB3q1YB1", l)], ["dB3q1YB1"]) // containers::python_for_data_science.oci
(_AmC7DKKQ) = o.postIn([in("inputs/AmC7DKKQ", l)], ["AmC7DKKQ"]) // containers::samtools.oci
(_rRC4FxHn) = o.postIn([in("inputs/rRC4FxHn", l)], ["rRC4FxHn"]) // containers::stringtie.oci
(_Rsoo8XE7) = o.postIn([in("inputs/Rsoo8XE7", l)], ["Rsoo8XE7"]) // transcriptomics::experiment
(_VwXt5dAr) = o.postIn([in("inputs/VwXt5dAr", l)], ["VwXt5dAr"]) // sequences::assembly
(_k0BWqyxf) = o.postIn([in("inputs/k0BWqyxf", l)], ["k0BWqyxf"]) // transcriptomics::star_bam
(_X81sU15q) = o.postIn([in("inputs/X81sU15q", l)], ["X81sU15q"]) // transcriptomics::stringtie_gtf
k = ['b83cdTFq']
(_b83cdTFq) = o.post([*p01__downloadBuscoLineage(o.group('OvFEfok5', [_H7g3InfC, _OvFEfok5], k, 1))], k)
k = ['iNlpm1XR']
(_iNlpm1XR) = o.post([*p02__merge_bams(o.group('Rsoo8XE7', [_Rsoo8XE7, _k0BWqyxf, _AmC7DKKQ], k, 1))], k)
k = ['3J2d8bJE', 'N3PPDtPP']
(_3J2d8bJE, _N3PPDtPP) = o.post([*p03__braker3(o.group('Rsoo8XE7', [_Rsoo8XE7, _VwXt5dAr, _iNlpm1XR, _ea9yJktW], k, 1))], k)
k = ['Myb7Pxjl']
(_Myb7Pxjl) = o.post([*p04__stringtie_merge(o.group('Rsoo8XE7', [_Rsoo8XE7, _X81sU15q, _3J2d8bJE, _rRC4FxHn], k, 1))], k)
k = ['LmS8pDph']
(_LmS8pDph) = o.post([*p05__gffread_proteins(o.group('Rsoo8XE7', [_Rsoo8XE7, _3J2d8bJE, _VwXt5dAr, _56zi3dIJ], k, 1))], k)
k = ['pDpp0Ttl']
(_pDpp0Ttl) = o.post([*p06__stringtie_quant(o.group('k0BWqyxf', [_k0BWqyxf, _Myb7Pxjl, _rRC4FxHn], k, 1))], k)
k = ['yTyedsZ0']
(_yTyedsZ0) = o.post([*p07__busco(o.group('LmS8pDph', [_H7g3InfC, _LmS8pDph, _b83cdTFq], k, 1))], k)
k = ['F0MqIpV6']
(_F0MqIpV6) = o.post([*p08__eggnog_mapper(o.group('LmS8pDph', [_0uxln7h5, _LmS8pDph, _MUbLkijl], k, 1))], k)
k = ['OJTSKg7W']
(_OJTSKg7W) = o.post([*p09__pydeseq2(o.group('Rsoo8XE7', [_Rsoo8XE7, _pDpp0Ttl, _Z9xw077K], k, 1))], k)
k = ['5TR9TDph']
(_5TR9TDph) = o.post([*p10__stringtie_count_matrix(o.group('Rsoo8XE7', [_Rsoo8XE7, _pDpp0Ttl, _dB3q1YB1], k, 1))], k)

publish:
_OJTSKg7W = o.publish(_OJTSKg7W)
_5TR9TDph = o.publish(_5TR9TDph)
_F0MqIpV6 = o.publish(_F0MqIpV6)
_3J2d8bJE = o.publish(_3J2d8bJE)
_yTyedsZ0 = o.publish(_yTyedsZ0)
}

output {
	_3J2d8bJE{
		path 'transcriptomics-braker3_gff'
		index { path '_manifests/transcriptomics-braker3_gff.3J2d8bJE.json' }
	}
	_yTyedsZ0{
		path 'annotation-busco_results'
		index { path '_manifests/annotation-busco_results.yTyedsZ0.json' }
	}
	_F0MqIpV6{
		path 'annotation-eggnog_results'
		index { path '_manifests/annotation-eggnog_results.F0MqIpV6.json' }
	}
	_OJTSKg7W{
		path 'transcriptomics-diff_count_table'
		index { path '_manifests/transcriptomics-diff_count_table.OJTSKg7W.json' }
	}
	_5TR9TDph{
		path 'transcriptomics-gene_count_table'
		index { path '_manifests/transcriptomics-gene_count_table.5TR9TDph.json' }
	}
}