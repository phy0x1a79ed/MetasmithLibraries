#!/bin/bash
# dev script version 1.1
HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# this file contains a list of commands useful for dev,
# providing automation for some build tasks
#
# example workflow 1, pip:
# dev.sh --idev # create a local conda dev env
# # add pypi api token as file to ./secrets [https://pypi.org/help/#apitoken]
# # make some changes to source
# # bump up ./src/*/version.txt
# dev.sh -bp # build the pip package
# dev.sh -up # test upload to testpypi
# dev.sh -upload-pypi # release to pypi index for pip install
#
# example workflow 2, conda:
# dev.sh --idev # create a local conda dev env
# dev.sh -bp # build the pip package
# dev.sh -bc # build conda package from pip package
# dev.sh -uc # publish to conda index
#
# example workflow 3, containerization:
# dev.sh --idev # create a local conda dev env
# dev.sh -bd # build docker image
# dev.sh -ud # publish to quay.io
# dev.sh -bs # build apptainer image from local docker image

case $1 in
    ###################################################
    # environments

    --ibase) # base only
        cd $HERE/envs
        echo "creating new conda env: $NAME"
        sleep 2
        $CONDA env create --no-default-packages -n $NAME -f ./base.yml
    ;;
    --git-prune-local) # remove local branches not on remote
        git fetch -p
        git branch -r \
            | awk '{print $1}' \
            | egrep -v -f /dev/fd/0 <(git branch -vv \
            | grep origin) \
            | awk '{print $1}' \
            | xargs git branch -d
    ;;

    ###################################################
    # build
    -b) # update std xgdbs
        $(which msm) && msm=msm || msm="$HERE/../Metasmith/dev.sh -r"
        echo $msm
        $msm build \
        --types $HERE/data_types \
        --uniques $HERE/resources/* \
        --transforms $HERE/transforms/*
    ;;
    ###################################################
    # test
    --test-binning)
        pytest tests/test_*.py -v --ignore=tests/cache
    ;;
    --test-comebin)
        pytest tests/test_binning_workflow.py::TestBinningWorkflowExecution::test_comebin_e2e \
            -v --ignore=tests/cache \
            -s --log-cli-level=INFO
    ;;
    --test-semibin2)
        pytest tests/test_binning_workflow.py::TestBinningWorkflowExecution::test_semibin2_e2e \
            -v --ignore=tests/cache \
            -s --log-cli-level=INFO
    ;;
    --test-metabat2)
        pytest tests/test_binning_workflow.py::TestBinningWorkflowExecution::test_metabat2_e2e \
            -v --ignore=tests/cache \
            -s --log-cli-level=INFO
    ;;
    --test-annotation)
        pytest tests/test_annotation_workflow.py
    ;;
    ###################################################
    *)
        echo "bad option"
        echo $1
    ;;
esac
