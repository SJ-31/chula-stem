#!/usr/bin/env bash

conf () {
    yq "$1" meta.yaml
}

get_dbsnp () {
    echo -e "RS,SAO,SSR,VC,NSF,NSM,NSN,SYN,U3,U5,PUB,FREQ,COMMON,CLNSIG,CLNDN,CLNREVSTAT,CLNACC\n" > "$2"
    string="%INFO/RS,%INFO/SAO,%INFO/SSR,%INFO/VC,%INFO/NSF,%INFO/NSM,%INFO/NSN,%INFO/SYN,%INFO/U3,%INFO/U5,%INFO/PUB,%INFO/FREQ,%INFO/COMMON,%INFO/CLNSIG,%INFO/CLNDN,%INFO/CLNREVSTAT,%INFO/CLNACC"
    bcftools query -f "$string" "$1" >> "$2"
}
