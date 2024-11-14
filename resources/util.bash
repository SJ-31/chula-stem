#!/usr/bin/env bash

conf () {
    yq "$1" meta.yaml
}
