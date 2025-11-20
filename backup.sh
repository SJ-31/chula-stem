#!/usr/bin/env bash

backup_root="/data/project/stemcell/shannc/backups"
stem_root="/data/project/stemcell/shannc"
home="/data/home/shannc"

date=$(date -I)
logfile="${backup_root}/backup_${date}.log"

# TODO: add a command to send the entire backup directory to synology with rsync, then
# delete it afterward. You don't have enough space to store everything here

echo "Backup started" >"${logfile}"

backup_fn() {
    archive="${3}/${1}.tar.gz"
    snar="${3}/.${1}.snar"
    source="${2}/${1}"
    tar --verbose \
        --create \
        --file "${archive}" \
        --listed-incremental "${snar}" \
        --directory "${2}" \
        "${1}"

    echo "Backup of ${source} to ${archive} completed $(date -I)" >>"${logfile}"
}

# backup_fn "public_data" "${stem_root}" "${backup_root}"
# backup_fn "reference" "${stem_root}" "${backup_root}"

# Repos
backup_fn "amr_predict" "${home}" "${backup_root}"
backup_fn "chula-stem" "${home}" "${backup_root}"
backup_fn "too-predict" "${home}" "${backup_root}"

# amr_predict data
backup_fn "output" "${stem_root}/repos/amr_predict" \
    "${backup_root}/repo_data/amr_predict"
backup_fn "datasets" "${stem_root}/repos/amr_predict" \
    "${backup_root}/repo_data/amr_predict"
backup_fn "genomes" "${stem_root}/repos/amr_predict" \
    "${backup_root}/repo_data/amr_predict"

# Outputs
backup_fn "CRC" "${stem_root}/output" "${backup_root}/output"
backup_fn "HCC" "${stem_root}/output" "${backup_root}/output"
backup_fn "PDAC" "${stem_root}/output" "${backup_root}/output"
backup_fn "PHCase" "${stem_root}/output" "${backup_root}/output"
backup_fn "WES_PON" "${stem_root}/output" "${backup_root}/output"
backup_fn "PDAC_TCR_2025-09-24" "${stem_root}/output" "${backup_root}/output"
