#!/bin/bash
# Log CPU/GPU/NVMe temps + fans + load every INTERVAL s to CSV; exit with ALERT if CPU Tctl
# stays >= LIMIT for NHIT consecutive samples (so the caller gets notified).
# Usage: monitor_temps_v1.sh <csv> [INTERVAL=60] [LIMIT=97] [NHIT=3]
OUT=$1; INT=${2:-60}; LIM=${3:-97}; NHIT=${4:-3}; hits=0
[ -f "$OUT" ] || echo "time,Tctl_C,gpu_edge_C,nvme_C,cpu_fan_rpm,load1" > "$OUT"
while true; do
    s=$(sensors 2>/dev/null)
    t=$(awk '/Tctl/ {gsub(/[+°C]/,"",$2); print $2}' <<<"$s")
    g=$(awk '/edge/ {gsub(/[+°C]/,"",$2); print $2; exit}' <<<"$s")
    n=$(awk '/Composite/ {gsub(/[+°C]/,"",$2); print $2; exit}' <<<"$s")
    f=$(awk '/cpu_fan/ {print $2; exit}' <<<"$s")
    l=$(cut -d' ' -f1 /proc/loadavg)
    echo "$(date +%F_%T),$t,$g,$n,$f,$l" >> "$OUT"
    if awk -v t="$t" -v L="$LIM" 'BEGIN{exit !(t+0>=L)}'; then hits=$((hits+1)); else hits=0; fi
    if [ $hits -ge $NHIT ]; then echo "ALERT: Tctl=${t}C >= ${LIM}C for $NHIT samples at $(date +%T)"; exit 2; fi
    sleep "$INT"
done
