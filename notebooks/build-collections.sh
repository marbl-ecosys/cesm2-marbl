#!/bin/bash
#SBATCH -J build-collection
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -p dav
#SBATCH -A NCGD0011
#SBATCH -t 24:00:00
#SBATCH --mem=1GB
#SBATCH -e build-collection.out
#SBATCH -o build-collection.out

# doesn't currently work with batch submission

overwrite=
while [[ $# -gt 0 ]]; do
  key="${1}"
  case ${key} in
    --overwrite-existing)
    overwrite=--overwrite-existing
    shift
    ;;
  *)
    echo ERROR
    exit
    ;;
  esac
done

source activate cesm2-marbl

collections=$(find intake-esm-collection-defs/* -name "*.yaml")

for collection in ${collections[@]}; do
  args="-cdef ${collection} ${overwrite}"
  echo intake-esm-builder ${args}
  intake-esm-builder ${args}
  echo
done
