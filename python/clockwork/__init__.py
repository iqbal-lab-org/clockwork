from pkg_resources import get_distribution

try:
        __version__ = get_distribution('clockwork').version
except:
    __version__ = 'local'


__all__ = [
    'common_data',
    'contam_remover',
    'cortex',
    'db',
    'db_connection',
    'db_maker',
    'db_schema',
    'isolate_dir',
    'ena',
    'fastqc',
    'fqtools',
    'het_snp_caller',
    'lock_file',
    'picard',
    'read_map',
    'read_pair_importer',
    'read_trim',
    'samtools_qc',
    'spreadsheet_importer',
    'spreadsheet_validator',
    'tasks',
    'utils',
]

from clockwork import *

