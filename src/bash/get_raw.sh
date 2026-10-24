#!/usr/bin/env bash

extract() {
    for i in $(tar tf "$1" | grep "[^a-z].fastq.gz"); do
        tar -xvf "$1" "$i"
    done
}

for tar in "${1}"/*.tar; do
    extract "$tar"
done
