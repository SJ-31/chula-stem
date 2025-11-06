#!/usr/bin/env bash

backup_root="/data/project/stemcell/shannc/backups"
stem_root="/data/project/stemcell/shannc"
home="/data/home/shannc"

date=$(date -I)
logfile="${backup_root}/backup_${date}.log"

echo "Backup started" >"${logfile}"

backup_fn() {
    archive="${3}/${1}.tar.gz"
    snar="${3}/.${1}.snar"
    source="${2}/${1}"
    tar --verbose \
        --create \
        --file "${archive}" \
        --listed-incremental "${snar}" \
        "${source}"

    echo "Backup of ${source} to ${archive} completed $(date -I)" >>"${logfile}"
}

backup_fn "public_data" "${stem_root}" "${backup_root}"
backup_fn "repos" "${stem_root}" "${backup_root}"
backup_fn "reference" "${stem_root}" "${backup_root}"

# Repos
backup_fn "amr_predict" "${home}" "${backup_root}"
backup_fn "chula-stem" "${home}" "${backup_root}"
backup_fn "too-predict" "${home}" "${backup_root}"

# Outputs
backup_fn "CRC" "${stem_root}/output" "${backup_root}/output"
backup_fn "HCC" "${stem_root}/output" "${backup_root}/output"
backup_fn "PDAC" "${stem_root}/output" "${backup_root}/output"
backup_fn "PHCase" "${stem_root}/output" "${backup_root}/output"
backup_fn "WES_PON" "${stem_root}/output" "${backup_root}/output"
backup_fn "PDAC_TCR_2025-09-24" "${stem_root}/output" "${backup_root}/output"
