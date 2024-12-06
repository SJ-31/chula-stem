#!/usr/bin/env bash

if [[ -z "$2" ]]; then
    cpus=1
else
    cpus="$2"
fi

srun --qos=cpu24h --mem=40G --cpus-per-task="${cpus}" nf-test test "$1"
