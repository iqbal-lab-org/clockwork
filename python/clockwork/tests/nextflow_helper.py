import shutil
import sys
import os

import clockwork



config_file = 'tmp.nextflow.config'


def clean_files():
    for filename in ['.nextflow.log', config_file]:
        try:
            os.unlink(filename)
        except:
            pass

    for directory in ['.nextflow']:
        try:
            shutil.rmtree(directory)
        except:
            pass

clean_files()


clockwork_python_modules_dir = os.path.dirname(os.path.abspath(clockwork.__file__))
clockwork_python_root = os.path.abspath(os.path.join(clockwork_python_modules_dir, os.pardir))
clockwork_repo_root = os.path.abspath(os.path.join(clockwork_python_root, os.pardir))
clockwork_scripts_dir = os.path.join(clockwork_python_root, 'scripts')
nextflow_dir = os.path.join(clockwork_repo_root, 'nextflow')
sys.path.append(clockwork_python_root)
data_root_dir = os.path.join(clockwork_python_modules_dir, 'tests', 'data')
os.environ['PATH'] = os.path.join(clockwork_scripts_dir) + ':' + os.environ['PATH']
if 'PYTHONPATH' in os.environ:
    os.environ['PYTHONPATH'] = clockwork_python_root + ':' + os.environ['PYTHONPATH']
else:
    os.environ['PYTHONPATH'] = clockwork_python_root



def write_config_file():
    with open(config_file, 'w') as f:
        print('env.PYTHONPATH = "', clockwork_python_root, ':$PYTHONPATH"', sep='', file=f)
        print('env.PATH = "', os.path.join(clockwork_scripts_dir), ':$PATH"', sep='', file=f)
