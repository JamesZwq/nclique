#!/bin/bash
# Kill benchmark and ALL descendant processes recursively

kill_tree() {
    local pid=$1
    # Get all children first
    local children=$(pgrep -P "$pid" 2>/dev/null)
    for child in $children; do
        kill_tree "$child"
    done
    kill -9 "$pid" 2>/dev/null && echo "  killed $pid"
}

# Find all bench_v3_all.py processes
pids=$(pgrep -f "bench_v3_all.py" 2>/dev/null)
if [ -n "$pids" ]; then
    echo "Killing benchmark process trees:"
    for pid in $pids; do
        echo "  root: $pid"
        kill_tree "$pid"
    done
fi

# Safety: kill any remaining degeneracy_cliques
remaining=$(pgrep -f degeneracy_cliques 2>/dev/null)
if [ -n "$remaining" ]; then
    echo "Killing remaining degeneracy_cliques:"
    echo "$remaining" | xargs kill -9 2>/dev/null
    echo "  killed $(echo "$remaining" | wc -w) processes"
fi

sleep 1
leftover=$(pgrep -f "bench_v3_all\|degeneracy_cliques" 2>/dev/null | wc -l)
echo "Done. Remaining: $leftover"
