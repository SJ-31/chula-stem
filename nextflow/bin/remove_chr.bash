#!/usr/bin/env bash

input="${1}"
output="${2}"

bcftools view "${1}" | \
    sed -e 's/ID=chr/ID=/' -e 's/^chr//' | \
    bcftools view -O z > "${2}"
