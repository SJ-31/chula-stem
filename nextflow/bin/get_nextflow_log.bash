#!/usr/bin/env bash

filename="$1"
{ echo "== COMMAND.BEGIN =="; \
  cat .command.begin ; \
  echo -e "\n== COMMAND.OUT =="; \
  cat .command.out ; \
  echo -e "\n== COMMAND.LOG =="; \
  cat .command.log; \
  echo -e "\n== COMMAND.ERR ==" ; \
  cat .command.err ; \
  } > "$filename"


