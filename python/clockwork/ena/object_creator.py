import os
from clockwork.ena import submission_receipt, submit_files, xml_create


class Error(Exception):
    pass


class ObjectCreator:
    def __init__(
        self,
        ini_file,
        object_type,
        object_xml,
        object_alias,
        submit_alias,
        center_name,
        title,
        project_description=None,  # needed for project
        taxon_id=None,  # needed for sample
        study_accession=None,
        sample_accession=None,
        library_name=None,  # needed for experiment
        platform=None,
        instrument=None,  # needed for experiment
        attributes_dict=None,  # optional for sample and experiment
        experiment_accession=None,
        reads_1=None,
        md5_1=None,
        reads_2=None,
        md5_2=None,  # required for run
        use_test_server=False,
        unit_test=None,
        broker_name=None,
    ):
        self.object_type = object_type
        if object_type == "project":
            if project_description is None:
                raise Error("Must provide project_description for object_type=project")
            self.project_description = project_description
        elif object_type == "sample":
            if taxon_id is None:
                raise Error("Must provide taxon_id for object_type=sample")
            self.taxon_id = taxon_id
        elif object_type == "experiment":
            required = [
                study_accession,
                sample_accession,
                library_name,
                platform,
                instrument,
            ]
            if None in required:
                raise Error(
                    "Must provide study_accession, sample_accession, library_name, platform, instrument for object_type=experiment. Got: "
                    + ", ".join([str(x) for x in required])
                )
            self.study_accession = study_accession
            self.sample_accession = sample_accession
            self.library_name = library_name
            self.platform = platform
            self.instrument = instrument
        elif object_type == "run":
            required = [experiment_accession, reads_1, md5_1, reads_2, md5_2]
            if None in required:
                raise Error(
                    "Must provide experiment_accession, reads_1, md5_1, reads_2, md5_2 for object_type=run. Got: "
                    + ", ".join([str(x) for x in required])
                )
            self.experiment_accession = experiment_accession
            self.reads_1 = reads_1
            self.md5_1 = md5_1
            self.reads_2 = reads_2
            self.md5_2 = md5_2
        else:
            raise Error(
                'objct_type "' + object_type + '" not recognised. Cannot continue'
            )

        self.ini_file = os.path.abspath(ini_file)
        self.object_xml_absolute = os.path.abspath(object_xml)
        self.object_xml_dir, self.object_xml = os.path.split(self.object_xml_absolute)
        self.submission_xml = self.object_xml + ".submission.xml"
        self.receipt_xml = self.object_xml + ".submission_receipt.xml"
        self.object_alias = object_alias
        self.submit_alias = submit_alias
        self.center_name = center_name
        self.title = title
        self.attributes_dict = attributes_dict
        self.use_test_server = use_test_server
        self.unit_test = unit_test
        self.broker_name = broker_name

    def _make_xml_files(self):
        if self.object_type == "project":
            xml_create.make_project_xml(
                self.object_xml,
                self.object_alias,
                self.center_name,
                self.title,
                self.project_description,
            )
        elif self.object_type == "sample":
            xml_create.make_sample_xml(
                self.object_xml,
                self.object_alias,
                self.center_name,
                self.title,
                self.taxon_id,
                sample_attributes_dict=self.attributes_dict,
            )
        elif self.object_type == "experiment":
            xml_create.make_experiment_xml(
                self.object_xml,
                self.object_alias,
                self.center_name,
                self.title,
                self.study_accession,
                self.sample_accession,
                self.library_name,
                self.platform,
                self.instrument,
                experiment_attributes_dict=self.attributes_dict,
            )
        elif self.object_type == "run":
            xml_create.make_paired_fastq_run_xml(
                self.object_xml,
                self.object_alias,
                self.center_name,
                self.experiment_accession,
                self.reads_1,
                self.md5_1,
                self.reads_2,
                self.md5_2,
            )
        else:
            raise Error(
                'objct_type "' + self.object_type + '" not recognised. Cannot continue'
            )

        xml_create.make_add_submission_xml(
            self.submission_xml,
            self.submit_alias,
            self.center_name,
            self.object_xml,
            self.object_type,
            broker_name=self.broker_name,
        )

    def run(self):
        original_dir = os.getcwd()
        os.chdir(self.object_xml_dir)
        self._make_xml_files()
        files = {
            "SUBMISSION": self.submission_xml,
            self.object_type.upper(): self.object_xml,
        }
        submit_files.submit_xml_files(
            self.ini_file,
            self.receipt_xml,
            files=files,
            use_test_server=self.use_test_server,
            unit_test=self.unit_test,
            unit_test_obj_type=self.object_type,
        )
        self.submission_receipt = submission_receipt.SubmissionReceipt(self.receipt_xml)
        os.chdir(original_dir)
