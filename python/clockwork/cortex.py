import glob
import os
import tarfile
import shutil
import pyfastaq
from clockwork import utils

class Error (Exception): pass


def _find_binaries(cortex_root=None, stampy_script=None, vcftools_dir=None):
    if cortex_root is None:
        cortex_root = os.environ.get('CLOCKWORK_CORTEX_DIR', '/bioinf-tools/cortex')
    cortex_root = os.path.abspath(cortex_root)

    if not os.path.isdir(cortex_root):
        raise Error('Cortex directory ' + cortex_root + ' not found. Cannot continue')

    cortex_run_calls = os.path.join(cortex_root, 'scripts', 'calling', 'run_calls.pl')
    if not os.path.exists(cortex_run_calls):
        raise Error('Cortex run_call.pl script ' + cortex_run_calls + ' not found. Cannot continue')

    if stampy_script is None:
        stampy_script = os.environ.get('CLOCKWORK_STAMPY_SCRIPT', '/bioinf-tools/stampy-1.0.31/stampy.py')
    stampy_script = os.path.abspath(stampy_script)

    if not os.path.exists(stampy_script):
        raise Error('Stampy script ' + stampy_script + ' not found. Cannot contine')

    if vcftools_dir is None:
        vcftools_dir = os.environ.get('CLOCKWORK_VCFTOOLS_DIR', '/bioinf-tools/vcftools-0.1.15/')
    vcftools_dir = os.path.abspath(vcftools_dir)

    if not os.path.isdir(vcftools_dir):
        raise Error('vcftools directory ' + vcftools_dir + ' not found. Cannot continue')

    return cortex_root, cortex_run_calls, stampy_script, vcftools_dir


def make_run_calls_index_files(ref_fasta, outprefix, cortex_root=None, stampy_script=None, vcftools_dir=None, mem_height=22):
    cortex_root, cortex_run_calls, stampy_script, vcftools_dir = _find_binaries(cortex_root=cortex_root, stampy_script=stampy_script, vcftools_dir=vcftools_dir)
    cortex_var = os.path.join(cortex_root, 'bin', 'cortex_var_31_c1')
    if not os.path.exists(cortex_var):
        raise Error('Cortex script not found: ' + cortex_var)

    fofn = outprefix + '.fofn'
    with open(fofn, 'w') as f:
        print(os.path.abspath(ref_fasta), file=f)

    utils.syscall(' '.join([
        cortex_var,
        '--kmer_size 31',
        '--mem_height', str(mem_height),
        '--mem_width 100',
        '--se_list', fofn,
        '--max_read_len 10000',
        '--dump_binary', outprefix + '.k31.ctx',
        '--sample_id REF',
    ]))

    os.unlink(fofn)

    utils.syscall(' '.join([
        stampy_script,
        '-G', outprefix + '.stampy',
        ref_fasta,
    ]))


    utils.syscall(' '.join([
        stampy_script,
        '-g', outprefix + '.stampy',
        '-H', outprefix + '.stampy',
    ]))


class CortexRunCalls:
    def __init__(self, ref_dir, reads_infile, outdir, sample_name, cortex_root=None, stampy_script=None, vcftools_dir=None, mem_height=22):
        self.ref_dir = os.path.abspath(ref_dir)
        self.reads_infile = os.path.abspath(reads_infile)
        self.outdir = os.path.abspath(outdir)
        self.sample_name = sample_name
        self.cortex_root = cortex_root
        self.stampy_script = stampy_script
        self.vcftools_dir = vcftools_dir
        self.cortex_log = os.path.join(self.outdir, 'cortex.log')

        if os.path.exists(self.outdir):
            raise Error('Error! Output directory already exists ' + self.outdir)


        self.cortex_outdir = os.path.join(self.outdir, 'cortex.out')
        self.cortex_reads_fofn = os.path.join(self.outdir, 'cortex.in.fofn')
        self.cortex_reads_index = os.path.join(self.outdir, 'cortex.in.index')
        self.cortex_ref_fofn = os.path.join(self.outdir, 'cortex.in.index_ref.fofn')
        self.cortex_root, self.cortex_run_calls, self.stampy_script, self.vcftools_dir = _find_binaries(cortex_root=cortex_root, stampy_script=stampy_script, vcftools_dir=vcftools_dir)
        self.mem_height = mem_height


    def _make_input_files(self):
        try:
            os.mkdir(self.outdir)
        except:
            raise Error('Error making directory ' + self.outdir + '. Cannot continue')

        with open(self.cortex_reads_fofn, 'w') as f:
            print(self.reads_infile, file=f)

        with open(self.cortex_reads_index, 'w') as f:
            print(self.sample_name, self.cortex_reads_fofn, '.', '.', sep='\t', file=f)

        with open(self.cortex_ref_fofn, 'w') as f:
            print(os.path.join(self.ref_dir, 'ref.fa'), file=f)


    def _tidy_files(self):
        shutil.rmtree(os.path.join(self.cortex_outdir, 'tmp_filelists'))

        for filename in glob.glob(os.path.join(self.cortex_outdir, 'binaries', 'uncleaned', '**', '**')):
            if not (filename.endswith('log') or filename.endswith('.covg')):
                os.unlink(filename)

        for filename in glob.glob(os.path.join(self.cortex_outdir, 'calls', '**')):
            if os.path.isdir(filename):
                for filename2 in os.listdir(filename):
                    if not filename2.endswith('log'):
                        os.unlink(os.path.join(filename, filename2))
            elif not (filename.endswith('log') or filename.endswith('callsets.genotyped')):
                os.unlink(filename)

        for filename in glob.glob(os.path.join(self.cortex_outdir, 'vcfs', '**')):
            if not ((filename.endswith('.vcf') and 'FINAL' in filename) or filename.endswith('log') or filename.endswith('aligned_branches')):
                if os.path.isdir(filename):
                    shutil.rmtree(filename)
                else:
                    os.unlink(filename)


    def run(self):
        self._make_input_files()
        ref_fai = os.path.join(self.ref_dir, 'ref.fa.fai')
        genome_size = pyfastaq.tasks.stats_from_fai(ref_fai)['total_length']

        cmd = ' '.join([
            self.cortex_run_calls,
            '--fastaq_index', self.cortex_reads_index,
            '--auto_cleaning yes',
            '--first_kmer 31',
            '--bc yes',
            '--pd no',
            '--outdir', self.cortex_outdir,
            '--outvcf cortex',
            '--ploidy 2',
            '--stampy_hash', os.path.join(self.ref_dir, 'ref.stampy'),
            '--stampy_bin', self.stampy_script,
            '--list_ref_fasta', self.cortex_ref_fofn,
            '--refbindir', self.ref_dir,
            '--genome_size', str(genome_size),
            '--qthresh 5',
            '--mem_height', str(self.mem_height),
            '--mem_width 100',
            '--vcftools_dir', self.vcftools_dir,
            '--do_union yes',
            '--ref CoordinatesAndInCalling',
            '--workflow independent',
            '--logfile', self.cortex_log,
        ])

        utils.syscall(cmd)
        self._tidy_files()
