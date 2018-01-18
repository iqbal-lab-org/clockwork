import os
import shutil
import tempfile
from clockwork import utils

class Error (Exception): pass

# The picard jar file in the singularity container is /bioinf-tools/picard.jar.
# If we're not using the container, then need to set the env variable
# CLOCKWORK_PICARD_JAR instead
PICARD_JAR = os.environ.get('CLOCKWORK_PICARD_JAR', '/bioinf-tools/picard.jar')

if not os.path.exists(PICARD_JAR):
    raise Error('Picard jar file not found. Please set environment variable CLOCKWORK_PICARD_JAR, or put it here: /bioinf-tools/picard.jar')


def mark_duplicates(sorted_bam_in, bam_out, xmx=2):
    tmpdir = tempfile.mkdtemp(prefix=bam_out + '.tmp.markdups.', dir=os.path.dirname(bam_out))
    m_file = os.path.join(tmpdir, 'covfefe')

    cmd = ' '.join([
        'java',
        '-Xmx' + str(xmx) + 'g',
        '-jar', PICARD_JAR,
        'MarkDuplicates',
        'VALIDATION_STRINGENCY=LENIENT',
        'INPUT=' + sorted_bam_in,
        'OUTPUT=' + bam_out,
        'M=' + m_file
    ])

    try:
        utils.syscall(cmd)
    except:
        shutil.rmtree(tmpdir)
        raise Error('Error runnin mark_duplicates: ' + cmd)

    shutil.rmtree(tmpdir)
