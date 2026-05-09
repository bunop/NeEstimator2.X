#!/usr/bin/env bash
set -euo pipefail

run_case() {
  local bin="$1"
  local tag="$2"
  local out_dir="$3"

  rm -f test_Ne.txt test_NexLD.txt
  "./${bin}" i:test_info.txt o:test_option.txt > "${out_dir}/${tag}.stdout" 2>&1

  if [[ -f test_Ne.txt ]]; then
    cp test_Ne.txt "${out_dir}/${tag}.test_Ne.txt"
  fi
  if [[ -f test_NexLD.txt ]]; then
    cp test_NexLD.txt "${out_dir}/${tag}.test_NexLD.txt"
  fi
}

work_dir=".regression"
baseline_dir="${work_dir}/baseline"
current_dir="${work_dir}/current"

rm -rf "${work_dir}"
mkdir -p "${baseline_dir}" "${current_dir}"

echo "[1/4] Build binaries"
make Ne2x Ne2x_cpp >/dev/null

echo "[2/4] Run baseline (C)"
run_case "Ne2x" "c" "${baseline_dir}"

echo "[3/4] Run wrapper (C++)"
run_case "Ne2x_cpp" "cpp" "${current_dir}"

echo "[4/4] Compare outputs"
cmp -s "${baseline_dir}/c.stdout" "${current_dir}/cpp.stdout"
cmp -s "${baseline_dir}/c.test_Ne.txt" "${current_dir}/cpp.test_Ne.txt"

if [[ -f "${baseline_dir}/c.test_NexLD.txt" || -f "${current_dir}/cpp.test_NexLD.txt" ]]; then
  cmp -s "${baseline_dir}/c.test_NexLD.txt" "${current_dir}/cpp.test_NexLD.txt"
fi

echo "PASS: C and C++ wrapper outputs are identical for test_info/test_option."
