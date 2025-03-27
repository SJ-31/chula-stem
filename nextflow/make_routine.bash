#!/usr/bin/env bash

dir="${1}"

if [[ -z "${dir}" ]]; then
    echo "No directory to place routine specified!"
    exit 1
fi

if [[ -e "${dir}"  ]]; then
    echo "Routine dir exists"
    exit 0
fi

mkdir "${dir}"

echo "includeConfig \"\$projectDir/main.config.nf\"" > "${dir}/nextflow.config"
echo "workDir = ???" >> "${dir}/nextflow.config"
ln -sr ./nextflow.config "${dir}/main.config.nf"

ln -sr ./main.nf "${dir}"
ln -sr ./modules "${dir}"
ln -sr ./bin "${dir}"
ln -sr ./lib "${dir}"
ln -sr ./reports "${dir}"
ln -sr ./subworkflows "${dir}"
ln -sr ./config "${dir}"
ln -sr ./input "${dir}"
ln -sr ./workflows "${dir}"
