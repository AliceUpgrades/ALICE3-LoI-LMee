#!/usr/bin/env bash
if [ -s *.root ]; then #check that we have a non empty root file.
  exit 0;
else
  exit 1;
fi
