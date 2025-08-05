#!/usr/bin/env bash

root="${1}"
prefix="${2}"

mkdir "${root}"
mkdir "${root}"/ignored
mkdir "${root}"/false_b

mkdir "${root}"/"${prefix}"-{1,2,3}
touch "${root}"/"${prefix}"-{1,2,3}/a.txt
mkdir "${root}"/"${prefix}"-{1,2,3}/foobar
touch "${root}"/"${prefix}"-{1,2,3}/foobar/b.txt
