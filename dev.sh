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
    --create-envs) # create the per-tool conda test envs from envs/tools/*.yml (idempotent)
        # each recipe's `name:` is the env name referenced by resources/env/<tool>.env `conda:`.
        # override the conda frontend with CONDA=... (default: mamba).
        CONDA=${CONDA:-mamba}
        existing="$($CONDA env list | awk '{print $1}')"
        for recipe in "$HERE"/envs/tools/*.yml; do
            name=$(awk -F': *' '/^name:/{print $2; exit}' "$recipe")
            if echo "$existing" | grep -qxF "$name"; then
                echo "[skip] env exists: $name"
                continue
            fi
            echo "[create] $name  <- $(basename "$recipe")"
            $CONDA env create -n "$name" -f "$recipe" \
                || echo "[WARN] failed to solve/create env: $name (recipe $recipe)"
        done
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
        # msm build's STEP positional must precede the flags; --types,
        # --uniques, --transforms are now single-value/repeatable. Build
        # the flag list by repeating each flag once per resolved path.
        if command -v msm >/dev/null 2>&1; then
            msm=msm
        else
            msm="$HERE/../Metasmith/dev.sh -r"
        fi
        echo "$msm"
        args=(build all --types "$HERE/data_types")
        for d in "$HERE"/resources/*/; do args+=(--uniques "${d%/}"); done
        for d in "$HERE"/transforms/*/; do args+=(--transforms "${d%/}"); done
        $msm "${args[@]}"
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
