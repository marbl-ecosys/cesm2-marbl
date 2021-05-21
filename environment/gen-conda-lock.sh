#!/usr/bin/env bash

# This script generates a conda lock file, which can be used to create an "exact" replica
# of our analysis kernal environment

conda_env=$(grep name: environment.yml | awk '{print $2}')

source activate ${conda_env}

# generate the lockfiles
conda-lock -f environment.yml -p osx-64 -p linux-64
