#!/usr/bin/env bash

# If you need to train a model pass t flag
train="false"
while getopts "t" f; do
	case "$f" in
	t)
		train="true"
		;;
	*) ;;
	esac
done

if [[ "${train}" == "true" ]]; then
	srun --qos=gpu40g --partition=gpu --gres=gpu:7g.40gb:1 --mem=30G \
		python ./expression_comparison_external.py
else
	srun --qos=cpu24h --mem=50G \
		python ./expression_comparison_external.py
fi

srun --qos=cpu24h --mem=30G \
	Rscript ./expression_comparison_external.R
