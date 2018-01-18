import os
from clockwork import db, fqtools, isolate_dir, utils

class Error (Exception): pass

class ReadPairImporter:
    def __init__(self,
        db_ini_file,
        pipeline_root_dir,
        seqrep_id,
        isolate_id,
        sample_id,
        sequence_replicate_number,
        reads_file_1,
        reads_file_2,
        reads_file_md5_1,
        reads_file_md5_2,
    ):
        self.db = db.Db(db_ini_file)
        self.pipeline_root_dir = os.path.abspath(pipeline_root_dir)
        self.seqrep_id = seqrep_id
        self.isolate_id = isolate_id
        self.sample_id = sample_id
        self.sequence_replicate_number = sequence_replicate_number
        self.reads_file_1 = os.path.abspath(reads_file_1)
        self.reads_file_2 = os.path.abspath(reads_file_2)
        self.reads_file_md5_1 = reads_file_md5_1
        self.reads_file_md5_2 = reads_file_md5_2


    @classmethod
    def _copy_reads_file(cls, old_name, new_name, expected_md5):
        old_name_md5 = utils.md5(old_name)
        if old_name_md5 != expected_md5:
            raise Error('MD5 given by submitter ' + expected_md5 + ' does not match calculated MD5 ' + old_name_md5)
        utils.rsync_and_md5(old_name, new_name, expected_md5)


    @classmethod
    def _check_database(cls, database, seqrep_id, isolate_id, sequence_replicate_number):
        rows = database.get_rows_from_table('Seqrep', where='seqrep_id = ' + str(seqrep_id))

        if len(rows) != 1:
            raise Error('Error! Got ' + str(len(rows)) + ' rows from Seqrep table with seqrep_id ' + str(seqrep_id))

        got_isolate = rows[0]['isolate_id']
        if got_isolate != isolate_id:
            raise Error('Error! Expected isolate_id ' + str(isolate_id) + ' but got ' + str(got_isolate) + ' for seqrep_id ' + str(seqrep_id))

        got_sequence_replicate_number = rows[0]['sequence_replicate_number']
        if got_sequence_replicate_number != sequence_replicate_number:
            raise Error('Error! Expected sequence_replicate_number ' + str(sequence_replicate_number) + ' but got ' + str(sequence_replicate_number) + ' for seqrep_id ' + str(seqrep_id))

        if rows[0]['import_status'] == -1:
            raise Error('Error! Import already tried and failed for seqrep ' + str(seqrep_id))
        elif rows[0]['import_status'] == 1:
            raise Error('Error! Already imported seqrep ' + str(seqrep_id))
        elif rows[0]['import_status'] != 0:
            raise Error('Error! import_status should be -1, 0 or 1, but is ' + str(rows[0]['import_status']) + ' for seqrep ' + str(seqrep_id))


    @classmethod
    def _update_database(cls, database, seqrep_id, isolate_id, sequence_replicate_number, import_status=1):
        ReadPairImporter._check_database(database, seqrep_id, isolate_id, sequence_replicate_number)

        try:
            database.update_row('Seqrep', {'seqrep_id': seqrep_id}, {'import_status': import_status})
        except:
            raise Error('Error updating database at end of import for seqrep ' + str(seqrep_id))


    def run(self):
        ReadPairImporter._check_database(self.db, self.seqrep_id, self.isolate_id, self.sequence_replicate_number)
        for filename in self.reads_file_1, self.reads_file_2:
            if not os.path.exists(filename):
                raise Error('Error! Reads file ' + filename + ' not found, Cannot continue.')

        iso_dir = isolate_dir.IsolateDir(self.pipeline_root_dir, self.sample_id, self.isolate_id)
        iso_dir.make_essential_dirs()

        lock_file = os.path.join(iso_dir.reads_dir, 'import_lock.' + str(self.seqrep_id))
        if os.path.exists(lock_file):
            raise Error('Error! Lock file ' + lock_file + ' found. Cannot continue')

        utils.make_empty_file(lock_file)

        try:
            fqtools.validate([self.reads_file_1, self.reads_file_2])
            ReadPairImporter._copy_reads_file(self.reads_file_1, iso_dir.reads_filename('original', self.sequence_replicate_number, 1), self.reads_file_md5_1)
            ReadPairImporter._copy_reads_file(self.reads_file_2, iso_dir.reads_filename('original', self.sequence_replicate_number, 2), self.reads_file_md5_2)
            ReadPairImporter._update_database(self.db, self.seqrep_id, self.isolate_id, self.sequence_replicate_number, import_status=1)
            os.unlink(self.reads_file_1)
            os.unlink(self.reads_file_2)
            for filename in (self.reads_file_1 + '.md5', self.reads_file_2 + '.md5'):
                if os.path.exists(filename):
                    os.unlink(filename)

            self.db.commit_and_close()
            os.unlink(lock_file)
        except:
            ReadPairImporter._update_database(self.db, self.seqrep_id, self.isolate_id, self.sequence_replicate_number, import_status=-1)
            self.db.commit_and_close()
            if os.path.exists(lock_file):
                os.unlink(lock_file)

