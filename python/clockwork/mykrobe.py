import json
import os
import shutil
import tempfile
from clockwork import utils

class Error (Exception): pass


def susceptibility_dict_from_json_file(json_file):
    '''Returns dict of susceptibility info from json file made by mykrobe predict'''
    with open(json_file) as f:
        json_data = json.load(f)

    sample_names = list(json_data.keys())
    if len(sample_names) != 1:
        raise Error('Expected one key in json file ' + json_file + ', but got: ' + str(sample_names))

    sample_name = sample_names[0]

    try:
        suscept_data = json_data[sample_name]['susceptibility']
    except:
        raise Error('Error getting susceptibility from file ' + json_file)

    return suscept_data


def run_predict(reads, outdir, sample_name, species, panel=None, custom_probe_and_json=None, mykrobe=None, clean=True):
    '''Runs mykrobe predict.
    reads = list of reads files.
    For a custom panel, use the option
    custom_probe_and_json, which should be a tuple of probe fasta file and
    json variant to resistance file. Using custom_probe_and_json will force
    panel='custom'.'''
    if mykrobe is None:
        mykrobe = os.environ.get('CLOCKWORK_MYKROBE', 'mykrobe')

    reads = [os.path.abspath(x) for x in reads]
    cwd = os.getcwd()
    os.mkdir(outdir)
    os.chdir(outdir)
    json_out = 'out.json'
    command = [mykrobe, 'predict', sample_name, species]

    if custom_probe_and_json is not None:
        assert len(custom_probe_and_json) == 2
        panel = 'custom'
        probe_file = os.path.abspath(custom_probe_and_json[0])
        json_file = os.path.abspath(custom_probe_and_json[1])
        assert os.path.exists(probe_file)
        assert os.path.exists(json_file)
        command += ['--custom_probe_set_path', probe_file,
                    '--custom_variant_to_resistance_json', json_file]

    if panel is not None:
        command += ['--panel', panel]

    command += ['--seq'] +  reads + ['>', json_out]
    command = ' '.join(command)
    completed_process = utils.syscall(command)
    with open('log.txt', 'w') as f:
        print('Command line:', command, file=f)
        print(completed_process.stdout, file=f)

    if clean:
        for d in 'atlas', 'tmp', 'mykrobe':
            try:
                shutil.rmtree(d)
            except:
                pass

    os.chdir(cwd)


class CustomPanel:
    def __init__(self, root_dir):
        self.root_dir = os.path.abspath(root_dir)
        self.probes_fasta = os.path.join(self.root_dir, 'probes.fa')
        self.var_to_res_json = os.path.join(self.root_dir, 'variant_to_resistance.json')
        self.species_file = os.path.join(self.root_dir, 'species.txt')
        self._set_species()


    def setup_files(self, species, probes_fasta, var_to_res_json):
        try:
            os.mkdir(self.root_dir)
        except:
            raise Error('Error mkdir ' + self.root_dir)

        utils.rsync_and_md5(probes_fasta, self.probes_fasta)
        utils.rsync_and_md5(var_to_res_json, self.var_to_res_json)
        with open(self.species_file, 'w') as f:
            print(species, file=f)
        self._set_species()


    def _set_species(self):
        if os.path.exists(self.species_file):
            with open(self.species_file) as f:
                self.species = f.readline().rstrip()
        else:
            self.species = None
