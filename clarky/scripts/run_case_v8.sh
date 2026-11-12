#!/bin/bash
# Usage: [SS_N=.. SS_DNUC=.. SS_CV=.. SS_CC=..] scripts/run_case_v8.sh <runName> <sigma> [nProcs=4] [endTime=0.35] [fine=1.0]
# Creates runs/<runName> from case_template_v5 (Huang 2013 setup, pressure ramp, kOmegaSSTReboud; needs clarky/src/kOmegaSSTReboud compiled), meshes, and runs interPhaseChangeFoam.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
NAME=$1; SIG=$2; NP=${3:-4}; END=${4:-0.35}; FINE=${5:-1.0}
[ -z "$SIG" ] && { echo "usage: $0 runName sigma [nProcs] [endTime] [fine]"; exit 1; }
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
source $WM_PROJECT_DIR/bin/tools/RunFunctions
C=$ROOT/runs/$NAME
[ -d "$C" ] && { echo "$C exists"; exit 1; }
mkdir -p $ROOT/runs && cp -r $ROOT/case_template_v5 $C && cd $C
# sed, not foamDictionary -set: the latter rewrites caseParams with #eval values (pOut) frozen
sed -i -E "s/^sigma( +)[0-9.]+;/sigma\1$SIG;/" constant/caseParams
foamDictionary -entry endTime -set $END system/controlDict >/dev/null
# v8: optional Schnerr-Sauer overrides via env: SS_N, SS_DNUC, SS_CV, SS_CC
for kv in "n:$SS_N" "dNuc:$SS_DNUC" "Cv:$SS_CV" "Cc:$SS_CC"; do
    k=${kv%%:*}; v=${kv#*:}
    [ -n "$v" ] && foamDictionary -entry SchnerrSauerCoeffs/$k -set $v constant/transportProperties >/dev/null
done
grep SchnerrSauerCoeffs constant/transportProperties
foamDictionary -entry numberOfSubdomains -set $NP system/decomposeParDict >/dev/null
python3 $ROOT/scripts/make_blockMesh_v2.py . --dat $ROOT/geometry/clarky.dat --fine $FINE --up 4 --down 5 --H 0.0945 | tee log.meshgen
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
