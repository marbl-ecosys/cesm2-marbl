import os
path_to_here = os.path.dirname(os.path.realpath(__file__))
project_tmpdir = '/glade/p/cgd/oce/projects/cesm2-marbl/funnel'


# TODO: need a "set" mechanism
units_info_yaml = f'{path_to_here}/units_info.yml'


# TODO: this should settable and shared across users
cache_database_file = f'{path_to_here}/data/cached_datasets.pkl'
