#!/bin/bash
# Benchmark Region Tuple vs ST across multiple graphs, r, s values
# Usage: bash bench_region.sh

BIN=./build/bin/degeneracy_cliques

GRAPHS=(
    "graphs/dblp-core30.edges"
    "graphs/email-Enron.edges"
    "graphs/web-Stanford.edges"
    "graphs/com-dblp.edges"
    "graphs/web-it-2004.edges"
)

# (r, s) pairs to test
RS_PAIRS=(
    "3 4"
    "3 5"
    "4 5"
    "5 6"
    "6 7"
)

echo "================================================================"
echo "Region Tuple vs ST Benchmark"
echo "================================================================"
printf "%-28s %5s %5s | %12s %12s | %8s | %s\n" \
    "Graph" "r" "s" "ST(ms)" "Tuple(ms)" "Speedup" "Match?"
echo "----------------------------------------------------------------"

for graph in "${GRAPHS[@]}"; do
    gname=$(basename "$graph" .edges)
    for rs in "${RS_PAIRS[@]}"; do
        r=$(echo $rs | awk '{print $1}')
        s=$(echo $rs | awk '{print $2}')

        # Run ST with 120s timeout
        st_out=$(timeout 120 bash -c "PIVOTER_RUN_ST=1 $BIN $graph $r $s 2>&1")
        st_exit=$?
        if [ $st_exit -ne 0 ]; then
            st_time="TIMEOUT"
            st_core=""
        else
            st_time=$(echo "$st_out" | grep "NucleusCoreDecomposition took" | grep -oP '[\d.]+')
            st_core=$(echo "$st_out" | grep "core=" | sort -t= -k2 -n)
        fi

        # Run Region Tuple with 300s timeout
        rt_out=$(timeout 300 bash -c "PIVOTER_RUN_REGION=1 $BIN $graph $r $s 2>&1")
        rt_exit=$?
        if [ $rt_exit -ne 0 ]; then
            rt_time="TIMEOUT"
            rt_core=""
            rt_tuples=""
            match="N/A"
        else
            rt_time=$(echo "$rt_out" | grep "NucleusCoreDecomposition took" | grep -oP '[\d.]+')
            rt_core=$(echo "$rt_out" | grep "core=" | sort -t= -k2 -n)
            rt_tuples=$(echo "$rt_out" | grep "Tuples:" | head -1 | grep -oP '\d+')

            # Compare core distributions
            if [ "$st_time" = "TIMEOUT" ]; then
                match="N/A(ST-TO)"
            else
                st_hash=$(echo "$st_core" | md5sum | cut -d' ' -f1)
                rt_hash=$(echo "$rt_core" | md5sum | cut -d' ' -f1)
                if [ "$st_hash" = "$rt_hash" ]; then
                    match="EXACT"
                else
                    # Check max core at least
                    st_max=$(echo "$st_core" | tail -1 | grep -oP 'core=\K\d+')
                    rt_max=$(echo "$rt_core" | tail -1 | grep -oP 'core=\K\d+')
                    if [ "$st_max" = "$rt_max" ]; then
                        match="maxOK($st_max)"
                    else
                        match="DIFF(ST=$st_max,RT=$rt_max)"
                    fi
                fi
            fi
        fi

        # Compute speedup
        if [ "$st_time" != "TIMEOUT" ] && [ "$rt_time" != "TIMEOUT" ] && [ -n "$st_time" ] && [ -n "$rt_time" ]; then
            speedup=$(echo "scale=1; $st_time / $rt_time" | bc 2>/dev/null || echo "?")
            speedup="${speedup}x"
        elif [ "$st_time" = "TIMEOUT" ] && [ "$rt_time" != "TIMEOUT" ]; then
            speedup=">inf"
        else
            speedup="?"
        fi

        printf "%-28s %5s %5s | %12s %12s | %8s | %s\n" \
            "$gname" "$r" "$s" \
            "${st_time:-?}" "${rt_time:-?}" \
            "$speedup" "$match"
    done
    echo ""
done
