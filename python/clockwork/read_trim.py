import os
from clockwork import utils

class Error (Exception): pass

def run_trimmomatic(
        reads1,
        reads2,
        out1,
        out2,
        trimmo_root=None,
        adapters='TruSeq3-PE-2.fa',
        minlen=50,
        verbose=0,
        threads=1,
        qual_trim='',
        adapters_included=True,
        quality_encoding='phred33'):

    if trimmo_root is None:
        trimmo_root = os.environ.get('CLOCKWORK_TRIMMO_DIR', '/bioinf-tools/Trimmomatic-0.36')
    trimmo_root = os.path.abspath(trimmo_root)
    jar_files = [x for x in os.listdir(trimmo_root) if x.endswith('.jar')]
    if len(jar_files) != 1:
        raise Error('Error finding Trimmoatic jar file in directory "' + trimmo_root + '". Found ' + str(len(jar_files)) + ' jar files. Cannot continue')
    jar_file = os.path.join(trimmo_root, jar_files[0])

    if adapters_included:
        adapters = os.path.join(trimmo_root, 'adapters', adapters)

    if not os.path.exists(adapters):
        raise Error('Cannot find adapters file "' + adapters + '".')

    cmd = ' '.join([
        'java -Xmx1000m -jar',
        jar_file,
        'PE',
        '-threads', str(threads),
        reads1,
        reads2,
        out1,
        '/dev/null',
        out2,
        '/dev/null',
        'ILLUMINACLIP:' + os.path.abspath(adapters) + ':2:30:10',
        qual_trim,
        'MINLEN:' + str(minlen),
        '-' + quality_encoding
    ])

    if verbose:
        print('Run trimmomatic:', cmd)
    utils.syscall(cmd)

