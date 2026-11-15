#!/bin/bash
# Usage: scripts/run_cascade_v3.sh <runName> <sigma> <alphaDeg> [nProcs=4] [endTime=0.30] [fine=1.0]
# Creates runs/<runName> from case_template_v3 (2D flat-plate cascade, beta 15 deg; inlet sigma via outlet-pressure feedback), meshes, runs
# interPhaseChangeFoam (needs clarky/src/kOmegaSSTReboud compiled). Wrap long runs in systemd-inhibit.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
NAME=$1; SIG=$2; ALPHA=$3; NP=${4:-4}; END=${5:-0.30}; FINE=${6:-1.0}
[ -z "$ALPHA" ] && { echo "usage: $0 runName sigma alphaDeg [nProcs] [endTime] [fine]"; exit 1; }
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
source $WM_PROJECT_DIR/bin/tools/RunFunctions
C=$ROOT/runs/$NAME
[ -d "$C" ] && { echo "$C exists"; exit 1; }
mkdir -p $ROOT/runs && cp -r $ROOT/case_template_v3 $C && cd $C
# sed, not foamDictionary -set (that freezes #eval-derived pOut/Ux/Uy)
sed -i -E "s/^sigma( +)[0-9.]+;/sigma\1$SIG;/; s/^alphaDeg( +)[0-9.]+;/alphaDeg\1$ALPHA;/" constant/caseParams
foamDictionary -entry endTime -set $END system/controlDict >/dev/null
foamDictionary -entry numberOfSubdomains -set $NP system/decomposeParDict >/dev/null
python3 $ROOT/scripts/make_cascade_blockMesh_v1.py . --fine $FINE | tee log.meshgen
echo "params: $(foamDictionary -expand constant/caseParams | grep -E '^(sigma|alphaDeg|pIn|kGain|Ux|Uy) ' | tr -s ' \n' ' ')" | tee -a log.meshgen
# fill controller constants into the coded outlet BC (numbers, since $vars are not expanded in code)
ev(){ foamDictionary -expand constant/caseParams | awk -v k="$1" '$1==k {gsub(";","",$2); print $2}'; }
sed -i "s/@PSTART@/$(ev pStart)/; s/@PIN@/$(ev pIn)/; s/@THOLD@/$(ev tHold)/; s/@TRAMP@/$(ev tRamp)/; s/@KGAIN@/$(ev kGain)/" 0.orig/p_rgh
grep -E "const scalar (pS|pT|tH|tR|k) " 0.orig/p_rgh | tr -s ' ' | tee -a log.meshgen
runApplication blockMesh
runApplication checkMesh
restore0Dir
runApplication potentialFoam -pName p_rgh
if [ "$NP" -gt 1 ]; then runApplication decomposePar; runParallel $(getApplication); else runApplication $(getApplication); fi
