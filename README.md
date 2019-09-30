# CESM2 Biogeochemistry Component: The Marine Biogeochemistry Library (MARBL)

Includes calculations support manuscript documenting MARBL in CESM2.

## To work with this repo

Use `conda` to create the environment.

```bash
conda env create -f environment.yaml
```

Project directory on glade:

```bash
project_dir=/glade/p/cgd/oce/projects/cesm2-marbl
```

To add a user to this directory
```bash
setfacl --modify default:user:USERNAME:rwx ${project_dir}
```