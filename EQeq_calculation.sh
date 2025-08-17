#!/bin/bash
# Env: ADTL/LINUX

FrameworkName=$1

mkdir -p "$FrameworkName"
cp "$FrameworkName.cif" "$FrameworkName/"

cd EQeq
chmod 777 eqeq

# Run EQeq
if ! ./eqeq "../$FrameworkName/$FrameworkName.cif" > "../$FrameworkName/${FrameworkName}_eqeq.log" 2>&1; then
    echo "$FrameworkName" >> ../../error_list_eqeq.txt
    exit 0
fi

# Check if output contains error marker
if grep -q "\-2147483648" "../$FrameworkName/${FrameworkName}_eqeq.log"; then
    echo "$FrameworkName" >> ../../error_list_eqeq.txt
    exit 0
fi

cd ..

# Update progress (with file lock)
(
    flock -n 9 || exit 1
    count=$(cat .charge_progress)
    count=$((count + 1))
    echo $count > .charge_progress
    echo -ne "Progress: $count charges calculated\r"
) 9> .charge_progress.lock
