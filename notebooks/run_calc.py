import os
import subprocess 
import sys
import pathlib

import json 
import yaml

import nbformat
import nbterm 

import traceback
import asyncio

import click

def get_main_notebooks():
    """return a list of notebook files"""
    with open('_config_calc.yml') as fid:
        return [
            f'{f}.ipynb' 
            for f in yaml.safe_load(fid)['notebooks']['main']
        ]


def get_pre_notebooks():
    """return the names of the preliminary computational notebooks in _config_calc.yml"""    
    with open('_config_calc.yml') as fid:
        return [
            f'{f}.ipynb' 
            for f in yaml.safe_load(fid)['notebooks']['pre_notebooks']
        ]

    
def get_project_kernel():
    """return the name of the project kernel stored in _config_calc.yml"""
    with open('_config_calc.yml') as fid:
        return yaml.safe_load(fid)['project_kernel']

def get_conda_kernel_cwd(name: str):
    """get the directory of a conda kernel by name"""
    command = ['conda', 'env', 'list', '--json']
    output = subprocess.check_output(command).decode('ascii')
    envs = json.loads(output)['envs']
    for env in envs:
        env = pathlib.Path(env)
        if name == env.stem:
            return env 

    else:
        return None


def nb_set_kernelname(file_in, kernel_name, file_out=None):
    """set the kernel name to python3"""
    if file_out is None:
        file_out = file_in        
    data = nbformat.read(file_in, as_version=nbformat.NO_CONVERT)        
    data['metadata']['kernelspec']['name'] = kernel_name
    nbformat.write(data, file_out)

    
def nb_get_kernelname(file_in):
    """get the kernel name of a notebook"""
    data = nbformat.read(file_in, as_version=nbformat.NO_CONVERT)
    return data['metadata']['kernelspec']['name']
    

def nb_clear_outputs(file_in, file_out=None):
    """clear output cells"""
    if file_out is None:
        file_out = file_in           
    data = nbformat.read(file_in, as_version=nbformat.NO_CONVERT)
    
    assert isinstance(data['cells'], list), 'cells is not a list'
    
    cells = []
    for cell in data['cells']:
        if cell['cell_type'] == 'code':
            cell['execution_count'] = None
            cell['outputs'] = []
        cells.append(cell)
    data['cells'] = cells
    nbformat.write(data, file_out)

    
def nb_execute_nbterm(notebook_path: str, kernel_cwd=None, output_dir=None):
    try:
        _nb_path = pathlib.Path(notebook_path)
        if not output_dir:
            output_dir = _nb_path.parent

        save_path = pathlib.Path(output_dir) / _nb_path.name
        nb = nbterm.Notebook(nb_path=_nb_path, save_path=save_path) #kernel_cwd=kernel_cwd, 
        asyncio.run(nb.run_all())
        nb.save(save_path)
        print(f"Executed notebook has been saved to: {save_path}")
        return True
    
    except Exception:
        msg = f'Error executing the notebook "{notebook_path}".\n'
        msg += f'See notebook "{notebook_path}" for the traceback.\n'
        print(f'{traceback.format_exc()}\n{msg}')
        return False

    
def nb_execute(notebook_filename, output_dir='.', kernel_name='python3'):
    """
    Execute a notebook.
    see http://nbconvert.readthedocs.io/en/latest/execute_api.html
    """
    import io
    import nbformat
    from nbconvert.preprocessors import ExecutePreprocessor
    from nbconvert.preprocessors import CellExecutionError

    #-- open notebook
    with io.open(notebook_filename, encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=nbformat.NO_CONVERT)

    # config for execution
    ep = ExecutePreprocessor(timeout=None, kernel_name=kernel_name)

    # run with error handling
    try:
        out = ep.preprocess(nb, {'metadata': {'path': './'}})

    except CellExecutionError:
        out = None
        msg = f'Error executing the notebook "{notebook_filename}".\n'
        msg += f'See notebook "{notebook_filename}" for the traceback.\n'
        print(msg)

    finally:
        nb_out = os.path.join(output_dir, os.path.basename(notebook_filename))
        with io.open(nb_out, mode='w', encoding='utf-8') as f:
            nbformat.write(nb, f)
        print(f'wrote: {nb_out}')
        
    return out    


def kernel_munge(kernel_name):
    """return the kernel name as it's rendered in the notebook metadata"""
    return f'conda-env-miniconda3-{kernel_name}-py'


@click.command()
@click.option('--notebook', default=None)
@click.option('--run-pre', is_flag=True)
@click.option('--stop-on-fail', is_flag=True)
def main(run_pre, notebook, stop_on_fail):
    """run notebooks"""
    
    project_kernel = get_project_kernel()
    
    assert os.environ['CONDA_DEFAULT_ENV'] == project_kernel, (
        f'activate "{project_kernel}" conda environment before running'
    )

    if notebook is None:
        notebook_list = get_pre_notebooks() if run_pre else []
        notebook_list = notebook_list + get_main_notebooks()
    else:
        notebook_list = [notebook]    
    
       
    # check kernels
    for nb in notebook_list:
        notebook_kernel = nb_get_kernelname(nb)
        assert notebook_kernel == kernel_munge(project_kernel), (
            f'{nb}: unexpected kernel: {notebook_kernel}'
        )
    
    # run the notebooks
    cwd = os.getcwd()
    failed_list = []
    for nb in notebook_list:
        print('-'*80)
        print(f'executing: {nb}')

        # set the kernel name to fool nbterm into running this
        nb_set_kernelname(nb, kernel_name='python3')

        # clear output
        nb_clear_outputs(nb)

        # run the notebook
        ok = nb_execute(nb, output_dir=cwd)
        
        # set the kernel back
        nb_set_kernelname(nb, kernel_name=kernel_munge(project_kernel))
        
        if not ok:
            print('failed')
            if stop_on_fail:
                sys.exit(1)            
            failed_list.append(nb)
        
        print()
        
    if failed_list:
        print('failed list')  
        print(failed_list)
        sys.exit(1)


if __name__ == '__main__':
    main()
    

    