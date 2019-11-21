#!/bin/bash

source /glade/work/mlevy/miniconda3/etc/profile.d/conda.sh

# Build netcdf file
rm -f ~/.intake_esm/collections/CESM1-CMIP5.nc
conda activate /glade/work/abanihi/softwares/miniconda3/envs/legacy-intake-esm/
intake-esm-builder -cdef glade-cesm1-cmip5-collection.yaml
conda deactivate

# Rebuild collection
conda activate cesm2-marbl
cd ..
jupyter nbconvert --to notebook --inplace \
                  --ExecutePreprocessor.kernel_name=python \
                  --ExecutePreprocessor.timeout=3600 \
                  --execute "build intake collections.ipynb"
