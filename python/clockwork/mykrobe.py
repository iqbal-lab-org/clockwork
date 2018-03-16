import os
import shutil
import tempfile
from clockwork import utils

class Error (Exception): pass



def _split_stdout_file_into_log_and_json(infile, out_log, out_json):
    with open(infile) as f_in, open(out_log, 'w') as f_log, open(out_json, 'w') as f_json:
        out_fh = f_log
        for line in f_in:
            if out_fh == f_log and line.startswith('{'):
                out_fh = f_json

            print(line, end='', file=out_fh)


def run_predict(reads, outdir, sample_name, species, panel=None, mykrobe=None):
    if mykrobe is None:
        mykrobe = os.environ.get('CLOCKWORK_MYKROBE', 'mykrobe')

    reads = os.path.abspath(reads)
    cwd = os.getcwd()
    os.mkdir(outdir)
    os.chdir(outdir)
    tmp_out = 'tmp.out'
    command = [mykrobe, 'predict']
    if panel is not None:
        command += ['--panel', panel]

    command += [
        sample_name,
        species,
        '-1', reads,
        '>', tmp_out
    ]
    command = ' '.join(command)
    utils.syscall(command)
    _split_stdout_file_into_log_and_json(tmp_out, 'log.txt', 'out.json')
    os.unlink(tmp_out)
    for d in 'atlas', 'tmp':
        try:
            shutil.rmtree(d)
        except:
            pass
    os.chdir(cwd)

