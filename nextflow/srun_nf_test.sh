#!/usr/bin/env bash

srun --qos=cpu24h --mem=40G nf-test test "$1"
