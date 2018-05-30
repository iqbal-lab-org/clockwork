import os
import pathlib

from clockwork import fqtools

class Error (Exception): pass


class FakeContamRemover:
    def __init__(self, reads_in_1, reads_in_2, reads_out_1, reads_out_2, counts_tsv_out):
        self.reads_in_1 = os.path.abspath(reads_in_1)
        self.reads_in_2 = os.path.abspath(reads_in_2)
        self.reads_out_1 = os.path.abspath(reads_out_1)
        self.reads_out_2 = os.path.abspath(reads_out_2)
        self.counts_tsv_out = os.path.abspath(counts_tsv_out)


    @classmethod
    def _write_counts_tsv(cls, reads1, reads2, outfile):
        number_of_read_pairs = fqtools.count([reads1, reads2])
        with open(outfile, 'w') as f:
            print('Name', 'Is_contam', 'Reads', sep='\t', file=f)
            print('Remove_contam_not_run', '0', 2 * number_of_read_pairs, sep='\t', file=f)
            print('Reads_kept_after_remove_contam', '0', 2 * number_of_read_pairs, sep='\t', file=f)


    @classmethod
    def _symlink_reads_file(cls, reads_file, symlink_path):
        '''Makes sylink called symlink_path pointing to reads_file.
        Files must be in the same directory, otherwise raises Error.
        The created symlink is a relative path -- this way the whole
        directory can be moved and the symlink is still valid'''
        if not os.path.exists(reads_file):
            raise Error('Error! file not found: ' + reads_file)

        symlink_path = os.path.abspath(symlink_path)
        symlink_path_dir, symlink_path_name = os.path.split(symlink_path)
        reads_file_dir, reads_file_name = os.path.split(os.path.abspath(reads_file))

        if reads_file_dir != symlink_path_dir:
            raise Error('Can only make symlinks inside same direcory. Filenames were:\n' + reads_file_dir + '\n' + symlink_path_dir)

        assert os.path.exists(symlink_path_dir)
        os.symlink(reads_file_name, symlink_path)
        assert pathlib.Path(symlink_path).exists()


    def run(self):
        FakeContamRemover._write_counts_tsv(self.reads_in_1, self.reads_in_2, self.counts_tsv_out)
        FakeContamRemover._symlink_reads_file(self.reads_in_1, self.reads_out_1)
        FakeContamRemover._symlink_reads_file(self.reads_in_2, self.reads_out_2)

