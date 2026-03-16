#!/bin/bash
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# >>> agent setup commands
module load gcc/9.4.0
module load apptainer/1.3.1
# <<<
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_DIR="./_metasmith/logs.$TIMESTAMP"
LOG_LATEST="./_metasmith/logs.latest"
mkdir -p $LOG_DIR
[ -e $LOG_LATEST ] && rm "$LOG_LATEST"; ln -s "./logs.$TIMESTAMP" "$LOG_LATEST"
[ -e workflow.params.yml ] || touch workflow.params.yml
[ -e workflow.config.nf ] || touch workflow.config.nf
echo "start time was [$TIMESTAMP]"
export BINDS="--bind /arc/project/st-shallam-1/pwy_group:/arc/project/st-shallam-1/pwy_group"
nohup ../../msm api run_workflow -a key=62Tq8B53 host=$(hostname) log_dir=$LOG_DIR stub_delay=${1:-0} >$LOG_DIR/agent.log 2>&1 &