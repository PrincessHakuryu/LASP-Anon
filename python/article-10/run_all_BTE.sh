#!/usr/bin/env bash
set -u
shopt -s nullglob

BTE_EXE="/home15/liangqm/ChihayaAnon-main/ShengBTE"

export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=close

dirs=(./*-BTE)

if [[ ${#dirs[@]} -eq 0 ]]; then
    echo "No *-BTE directories found in $(pwd)"
    exit 1
fi

echo "Found ${#dirs[@]} *-BTE directories."
echo

for d in "${dirs[@]}"; do
    echo "============================================================"
    echo "Running ShengBTE in: $d"
    echo "Start time: $(date)"
    echo "============================================================"

    if [[ ! -f "$d/CONTROL" ]]; then
        echo "WARNING: $d has no CONTROL, skip."
        continue
    fi

    if [[ ! -f "$d/FORCE_CONSTANTS_2ND" ]]; then
        echo "WARNING: $d has no FORCE_CONSTANTS_2ND, skip."
        continue
    fi

    if [[ ! -f "$d/FORCE_CONSTANTS_3RD" ]]; then
        echo "WARNING: $d has no FORCE_CONSTANTS_3RD, skip."
        continue
    fi

    (
        cd "$d" || exit 1

        mpirun -np 2 \
            --map-by numa:PE=8 \
            --bind-to core \
            --report-bindings \
            -x OMP_NUM_THREADS=8 \
            -x OMP_PLACES=cores \
            -x OMP_PROC_BIND=close \
            --mca coll_hcoll_enable 0 \
            "$BTE_EXE" > ShengBTE.run.log 2>&1
    )

    status=$?

    if [[ $status -eq 0 ]]; then
        echo "SUCCESS: $d"
        touch "$d/ShengBTE.done"
    else
        echo "FAILED : $d, exit code = $status"
        touch "$d/ShengBTE.failed"
    fi

    echo "End time: $(date)"
    echo
done

echo "All done."
