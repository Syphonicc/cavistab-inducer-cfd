#!/bin/bash
# Usage: scripts/run_case_v1.sh <runName> <sigma> [nProcs=4] [endTime=0.01] [fine=1.0]
# Creates runs/<runName> from case_template, meshes, and runs interPhaseChangeFoam.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
NAME=$1; SIG=$2; NP=${3:-4}; END=${4:-0.01}; FINE=${5:-1.0}
[ -z "$SIG" ] && { echo "usage: $0 runName sigma [nProcs] [endTime] [fine]"; exit 1; }
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
source $WM_PROJECT_DIR/bin/tools/RunFunctions
C=$ROOT/runs/$NAME
[ -d "$C" ] && { echo "$C exists"; exit 1; }
mkdir -p $ROOT/runs && cp -r $ROOT/case_template $C && cd $C
foamDictionary -entry sigma -set $SIG constant/caseParams >/dev/null
foamDictionary -entry endTime -set $END system/controlDict >/dev/null
foamDictionary -entry numberOfSubdomains -set $NP system/decomposeParDict >/dev/null
python3 $ROOT/scripts/make_blockMesh_v2.py . --dat $ROOT/geometry/clarky.dat --fine $FINE | tee log.meshgen
runApplication blockMesh
runApplication checkMesh
restore0Dir
runApplication potentialFoam -pName p_rgh
if [ "$NP" -gt 1 ]; then
    runApplication decomposePar
    runParallel $(getApplication)
else
    runApplication $(getApplication)
fi
