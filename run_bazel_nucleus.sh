#!/usr/bin/env bash
# run_bazel_nucleus.sh
# 顺序执行；任何一次 (r, s) 超时（或出错）都停止该 r 的后续 s，继续下一个 r。
# 超时设为 5 小时，/bin/time -v 计时信息写入独立日志。

set -u

export PARLAY_NUM_THREADS=1
TIME_BIN="/bin/time -v"
TIMEOUT_LIMIT="5h"     # 超时时间
BAZEL_TARGET=":NucleusDecomposition_main"
COMMON_FLAGS="-s -rounds 1 -efficient_inline"

# 你的数据集（优先用 .adj；若不存在则尝试 .edges）
DATASETS=(
  "/data1/wenqianz/web-it-2004.adj"
  "/data1/wenqianz/com-dblp.adj"
  "/data1/wenqianz/web-Google.adj"
  "/data1/wenqianz/soc-pokec-relationships.adj"
  "/data1/wenqianz/web-Stanford.adj"
  "/data1/wenqianz/com-youtube.adj"
)

# 若 .adj 不存在，则尝试同名 .edges
resolve_path() {
  local p="$1"
  if [[ -f "$p" ]]; then
    printf "%s" "$p"
    return 0
  fi
  # 尝试把 .adj 换成 .edges
  if [[ "$p" == *.adj ]]; then
    local q="${p%.adj}.edges"
    if [[ -f "$q" ]]; then
      printf "%s" "$q"
      return 0
    fi
  fi
  # 反过来：.edges -> .adj
  if [[ "$p" == *.edges ]]; then
    local q="${p%.edges}.adj"
    if [[ -f "$q" ]]; then
      printf "%s" "$q"
      return 0
    fi
  fi
  return 1
}

LOG_DIR="bazel_runs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"

# 先 build 一次，避免每次 run 都触发编译
echo "[INFO] bazel build ${BAZEL_TARGET}"
bazel build "${BAZEL_TARGET}" || {
  echo "[ERROR] bazel build 失败，退出。"
  exit 2
}

run_one() {
  local ds="$1" r="$2" s="$3" log="$4"
  echo "[INFO] dataset=$(basename "$ds"), r=${r}, s=${s} -> $log"
  # 用 nohup 保持会话外也能继续；不加后台 &，以便拿到退出码做超时判断
  nohup bash -lc "timeout ${TIMEOUT_LIMIT} ${TIME_BIN} bazel run ${BAZEL_TARGET} -- ${COMMON_FLAGS} -r ${r} -ss ${s} -efficient_inline '${ds}'" \
    >"${log}" 2>&1
  local ec=$?
  # GNU timeout 超时返回 124
  if [[ $ec -eq 124 ]]; then
    echo "[WARN] (r=${r}, s=${s}) 超时（>${TIMEOUT_LIMIT) }。"
  elif [[ $ec -ne 0 ]]; then
    echo "[WARN] (r=${r}, s=${s}) 非零退出码：$ec"
  fi
  return $ec
}

# r 的范围与 s 的范围（与你之前要求一致）：
# r=1 -> s=3..30
# r=2 -> s=4..30
# r=3 -> s=4..30
for raw_ds in "${DATASETS[@]}"; do
  ds=$(resolve_path "$raw_ds") || {
    echo "[WARN] 找不到数据集：$raw_ds（尝试 .adj/.edges 都失败），跳过。"
    continue
  }
  base="$(basename "$ds")"

  # r=1
#  r=1
#  for ((s=3; s<=30; ++s)); do
#    log="${LOG_DIR}/${base}.r${r}.s${s}.log"
#    run_one "$ds" "$r" "$s" "$log"
#    ec=$?
#    if [[ $ec -ne 0 ]]; then
#      echo "[INFO] r=${r} 在 s=${s} 处失败/超时，跳过 r=${r} 后续 s。转入 r=2。"
#      break
#    fi
#  done
#
#  # r=2
#  r=2
#  for ((s=4; s<=30; ++s)); do
#    log="${LOG_DIR}/${base}.r${r}.s${s}.log"
#    run_one "$ds" "$r" "$s" "$log"
#    ec=$?
#    if [[ $ec -ne 0 ]]; then
#      echo "[INFO] r=${r} 在 s=${s} 处失败/超时，跳过 r=${r} 后续 s。转入 r=3。"
#      break
#    fi
#  done
#
#  # r=3
#  r=3
#  for ((s=4; s<=30; ++s)); do
#    log="${LOG_DIR}/${base}.r${r}.s${s}.log}"
#    run_one "$ds" "$r" "$s" "$log"
#    ec=$?
#    if [[ $ec -ne 0 ]]; then
#      echo "[INFO] r=${r} 在 s=${s} 处失败/超时，跳过 r=${r} 后续 s。该数据集完成。"
#      break
#    fi
#  done
  r=4
  for ((s=5; s<=30; ++s)); do
    log="${LOG_DIR}/${base}.r${r}.s${s}.log}"
    run_one "$ds" "$r" "$s" "$log"
    ec=$?
    if [[ $ec -ne 0 ]]; then
      echo "[INFO] r=${r} 在 s=${s} 处失败/超时，跳过 r=${r} 后续 s。该数据集完成。"
      break
    fi
  done
done

echo "[DONE] 全部任务已提交并顺序执行完成。日志目录：${LOG_DIR}"