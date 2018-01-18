import tempfile
import shutil
import os
import pysam
from clockwork import fqtools, picard, utils

class Error (Exception): pass

def map_reads(ref_fasta, reads1, reads2, outfile, rmdup=False, markdup=False, read_group=None):
    '''Maps reads with BWA MEM. By default, outputs SAM file in input read order.
    rmdup=True => remove duplicates using samtools rmdup. Final output is sorted bam
                  Incompatible with markdup=True
    markdup=True => mark duplicates using picard MarkDuplicate. Final output is sorted bam.
                  Incompatible with rmdup=True
    read_group should be a tuple (group_id, group_name). If given, these will be
    put into the BAM'''
    if rmdup and markdup:
        raise Error('Cannot have rmdup and markdup both True.''')

    try:
        expected_read_count = 2 * fqtools.count([reads1, reads2])
    except:
        raise Error('Error counting reads in input files ' + reads1 + ' ' + reads2)

    if rmdup or markdup:
        tmpdir = tempfile.mkdtemp(prefix=outfile + '.tmp.map_reads.', dir=os.path.dirname(outfile))
        sam_file = os.path.join(tmpdir, 'tmp.sam')
    else:
        sam_file = outfile

    # "LB:LIB" is needed, otherwise samtools rmdup segfaults when map_reads_set() is used
    R_option = '' if read_group is None else r'''-R '@RG\tLB:LIB\tID:''' + read_group[0] + r'''\tSM:''' + read_group[1] + "'"

    cmd = ' '.join([
        'bwa mem -M',
        R_option,
        ref_fasta,
        reads1, reads2,
        r''' | awk '/^@/ || !and($2,256)' ''', # remove secondary alignments (but keep header)
        '>', sam_file
    ])

    try:
        utils.syscall(cmd)
    except:
        if rmdup or markdup:
            shutil.rmtree(tmpdir)
        raise Error('Error running BWA MEM: ' + cmd)

    number_in_sam = utils.sam_record_count(sam_file)
    if expected_read_count != number_in_sam:
        if rmdup or markdup:
            shutil.rmtree(tmpdir)
        raise Error('Error! Mismatch in read counts. Expected ' + str(expected_read_count) + ' but got ' + str(number_in_sam))


    if rmdup or markdup:
        sorted_bam = os.path.join(tmpdir, 'tmp.sorted.bam')

        cmd = ' '.join([
            'samtools sort',
            '-o', sorted_bam,
            sam_file
        ])

        try:
            utils.syscall(cmd)
        except:
            shutil.rmtree(tmpdir)
            raise Error('Error running samtools sort: ' + cmd)


        if rmdup:
            cmd = 'samtools rmdup ' + sorted_bam + ' ' + outfile
            try:
                utils.syscall(cmd)
            except:
                shutil.rmtree(tmpdir)
                raise Error('Error running samtools rmdup: ' + cmd)
        else:
            try:
                picard.mark_duplicates(sorted_bam, outfile)
            except:
                shutil.rmtree(tmpdir)
                raise Error('Error running picard mark_duplicates ' + cmd)

        shutil.rmtree(tmpdir)


def map_reads_set(ref_fasta, reads1_list, reads2_list, outfile, rmdup=False, markdup=False, read_group=None):
    '''Same as map reads, but takes a list of file pairs to be mapped'''
    assert len(reads1_list) == len(reads2_list)
    outfile = os.path.abspath(outfile)
    tmpdir = tempfile.mkdtemp(prefix=outfile + '.tmp.map_reads_set.', dir=os.path.dirname(outfile))
    outfiles = []

    for i in range(len(reads1_list)):
        outfiles.append(os.path.join(tmpdir, 'map.' + str(i)))
        map_reads(ref_fasta, reads1_list[i], reads2_list[i], outfiles[-1], rmdup=rmdup, markdup=markdup, read_group=read_group)

    assert len(outfiles) == len(reads1_list)

    if len(outfiles) == 1:
        os.rename(outfiles[-1], outfile)
    else:
        file_string = ', '.join(['outfiles[' + str(i) + ']' for i in range(len(outfiles))])
        exec('pysam.merge("-c", outfile, ' + file_string + ')')

    shutil.rmtree(tmpdir)
