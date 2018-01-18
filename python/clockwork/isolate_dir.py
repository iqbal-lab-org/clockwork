import os

class Error (Exception): pass

MAX_SAMPLE_ID = 99999999

class IsolateDir:
    def __init__(self, root_dir, sample_id, isolate_id):
        assert sample_id <= MAX_SAMPLE_ID
        self.root_dir = os.path.abspath(root_dir)
        self.sample_id = sample_id
        self.isolate_id = isolate_id

        self.sample_dir = IsolateDir._make_sample_dir_path(self.root_dir, self.sample_id)
        self.isolate_dir = os.path.join(self.sample_dir, str(isolate_id))
        self.reads_dir = os.path.join(self.isolate_dir, 'Reads')
        self.pipelines_dir = os.path.join(self.isolate_dir, 'Pipelines')


    @classmethod
    def sample_id_to_dir_numbers(cls, sample_id):
        assert sample_id <= MAX_SAMPLE_ID
        id_string = str(sample_id).zfill(len(str(MAX_SAMPLE_ID)))
        return [id_string[i:i+2] for i in range(0, len(id_string), 2)]


    @classmethod
    def _make_sample_dir_path(cls, root_dir, sample_id):
        dir_numbers = IsolateDir.sample_id_to_dir_numbers(sample_id)
        return os.path.join(root_dir, *dir_numbers)


    def make_essential_dirs(self):
        dirs = {
            'sample': self.sample_dir,
            'isolate': self.isolate_dir,
            'reads': self.reads_dir,
            'pipelines': self.pipelines_dir,
        }

        for name, directory in dirs.items():
            try:
                os.makedirs(directory)
            except FileExistsError:
                pass
            except:
                raise Error('Error making ' + name + ' directory: ' + directory)


    def contamination_counts_filename(self, sequence_replicate):
        return os.path.join(self.reads_dir, 'reads.remove_contam.' + str(sequence_replicate) + '.counts.tsv')


    def reads_filename(self, reads_type, sequence_replicate, one_or_two):
        if reads_type not in {'original', 'remove_contam', 'contam'}:
            raise Error('Error! Reads type "' + reads_type + '" not recognised. Cannot continue')

        if one_or_two not in {1, 2}:
            raise Error('Error! Read index must be 1 or 2, but got: ' + str(one_or_two))

        return os.path.join(self.reads_dir, 'reads.' + reads_type + '.' + str(sequence_replicate) + '.' + str(one_or_two) + '.fq.gz')


    def pipeline_dir(self, sequence_replicate, pipeline_name, pipeline_version, reference_id=None):
        final_dir = pipeline_version if reference_id is None else pipeline_version + '.ref.' + str(reference_id)
        return os.path.join(self.pipelines_dir, str(sequence_replicate), pipeline_name, final_dir)


    def xml_submission_file(self, xml_type, sequence_replicate=None):
        if xml_type == 'sample':
            return os.path.join(self.sample_dir, 'ena_sample_submission.xml')
        elif xml_type == 'experiment':
            return os.path.join(self.isolate_dir, 'ena_experiment_submission.xml')
        elif xml_type == 'run':
            if sequence_replicate is None:
                raise Error('Must provide sequence_replicate for xml_type "run"')
            return os.path.join(self.reads_dir, 'reads.remove_contam.' + str(sequence_replicate) + '.ena_run_submission.xml')
        else:
            raise Error('xml_type "' + xml_type + '" not recognised. Cannot continue')

