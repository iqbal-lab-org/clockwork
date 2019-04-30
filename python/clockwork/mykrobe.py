import json
import os
import shutil
import tempfile
from clockwork import utils

class Error (Exception): pass


built_in_panels = {
    'tb': {'bradley-2015', 'walker-2015'},
}

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


def run_predict(reads, outdir, sample_name, species, panel=None, custom_probe_and_json=None, mykrobe=None, clean=True, unittest=False, unittest_resistant=False):
    '''Runs mykrobe predict.
    reads = list of reads files.
    For a custom panel, use the option
    custom_probe_and_json, which should be a tuple of probe fasta file and
    json variant to resistance file. Using custom_probe_and_json will force
    panel='custom'.
    mykrobe is too slow to include in unit tests. Setting unittest=True will not run
    mykrobe, instead checking 'mykrobe predict' with no options does not return
    an error, and then writes dummy output files. If unittest_resistant=True, then
    it will make a report with resistance to isoniazid '''
    if mykrobe is None:
        mykrobe = os.environ.get('CLOCKWORK_MYKROBE', 'mykrobe')

    reads = [os.path.abspath(x) for x in reads]
    cwd = os.getcwd()
    os.mkdir(outdir)
    os.chdir(outdir)
    json_out = 'out.json'
    command = [mykrobe, 'predict', '--format', 'json', sample_name, species]

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
    if unittest:
        write_test_json_file(json_out, resistant=unittest_resistant)
        command = mykrobe +  ' predict --help'

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


class Panel:
    def __init__(self, root_dir):
        self.root_dir = os.path.abspath(root_dir)
        self.probes_fasta = os.path.join(self.root_dir, 'probes.fa')
        self.var_to_res_json = os.path.join(self.root_dir, 'variant_to_resistance.json')
        self.json_file = os.path.join(self.root_dir, 'data.json')
        self._load_json_file()


    def setup_files(self, species, panel_name, probes_fasta, var_to_res_json):
        try:
            os.mkdir(self.root_dir)
        except:
            raise Error('Error mkdir ' + self.root_dir)

        self.metadata = {
            'species': species,
            'name': panel_name,
            'is_built_in': species in built_in_panels and panel_name in built_in_panels[species],
        }

        with open(self.json_file, 'w') as f:
            print(json.dumps(self.metadata), file=f)

        if not self.metadata['is_built_in']:
            utils.rsync_and_md5(probes_fasta, self.probes_fasta)
            utils.rsync_and_md5(var_to_res_json, self.var_to_res_json)


    def _load_json_file(self):
        if os.path.exists(self.json_file):
            with open(self.json_file) as f:
                self.metadata = json.load(f)
        else:
            self.metadata = {'species': None, 'is_built_in': False, 'name': None}


def write_test_json_file(outfile, resistant=False):
    with open(outfile, 'w') as f:
        print(r'''{
    "sample_name": {
        "variant_calls": {
            "foo": "bar"
        },
        "susceptibility": {
            "Rifampicin": {
                "predict": "S"
            },
            "Capreomycin": {
                "predict": "S"
            },
            "Isoniazid": {''', file=f)

        if resistant:
            print(r'''                "predict": "R",
                "called_by": {
                    "inhA_I21T-I21T": {
                        "info": {
                            "filter": "PASS",
                            "contamination_depths": [],
                            "conf": 99999997,
                            "expected_depths": [
                                4
                            ],
                            "coverage": {
                                "alternate": {
                                    "percent_coverage": 100.0,
                                    "median_depth": 3.0,
                                    "kmer_count": 81,
                                    "min_non_zero_depth": 3.0
                                },
                                "reference": {
                                    "percent_coverage": 0.0,
                                    "median_depth": 0.0,
                                    "kmer_count": 0,
                                    "min_non_zero_depth": 0.0
                                }
                            }
                        },
                        "genotype": [
                            1,
                            1
                        ],
                        "variant": null,
                        "genotype_likelihoods": [
                            -99999999,
                            -99999999,
                            -1.6995430525350501
                        ],
                        "_cls": "Call.VariantCall"
                    }
                }
            },''', file=f)
        else:
            print(r'''                "predict": "S"
            },''', file=f)

        print(r'''            "Amikacin": {
                "predict": "S"
            },
            "Pyrazinamide": {
                "predict": "S"
            },
            "Ethambutol": {
                "predict": "S"
            },
            "Kanamycin": {
                "predict": "S"
            },
            "Streptomycin": {
                "predict": "S"
            },
            "Quinolones": {
                "predict": "S"
            }
        },
        "version": {
            "mykrobe-atlas": "v0.6.4-0-gb0374e1",
            "mykrobe-predictor": "v0.6.4-0-gb0374e1"
        },
        "kmer": 21,
        "sequence_calls": {},
        "genotype_model": "median_depth",
        "phylogenetics": {
            "lineage": {
                "European_American": {
                    "percent_coverage": 65.833,
                    "median_depth": 3.5
                }
            },
            "sub_complex": {
                "Unknown": {
                    "percent_coverage": -1,
                    "median_depth": -1
                }
            },
            "phylo_group": {
                "Mycobacterium_tuberculosis_complex": {
                    "percent_coverage": 95.435,
                    "median_depth": 4
                }
            },
            "species": {
                "Mycobacterium_tuberculosis": {
                    "percent_coverage": 94.177,
                    "median_depth": 4.0
                }
            }
        },
        "probe_sets": [
            "/usr/local/lib/python2.7/dist-packages/mykrobe/data/panels/tb-species-170421.fasta.gz",
            "/usr/local/lib/python2.7/dist-packages/mykrobe/data/panels/tb-walker-probe-set-feb-09-2017.fasta.gz"
        ],
        "files": [
            "/home/vagrant/clockwork/python/clockwork/tests/data/mykrobe/run_predict.reads.fq.gz"
        ]
    }
}''', file=f)
