import datetime
from xml.dom import minidom
import xml.etree.ElementTree as ET


def element_tree_to_file(tree, outfile):
    with open(outfile, "w") as f:
        print(
            minidom.parseString(ET.tostring(tree)).toprettyxml(indent="   "),
            end="",
            file=f,
        )


def make_add_submission_xml(
    outfile,
    submission_alias,
    center_name,
    source_xml,
    schema,
    broker_name=None,
    hold_for_two_years=True,
):
    if broker_name is None:
        submission = ET.Element(
            "SUBMISSION", {"alias": submission_alias, "center_name": center_name}
        )
    else:
        submission = ET.Element(
            "SUBMISSION",
            {
                "alias": submission_alias,
                "center_name": center_name,
                "broker_name": broker_name,
            },
        )
    actions = ET.SubElement(submission, "ACTIONS")
    action_add = ET.SubElement(actions, "ACTION")
    add = ET.SubElement(action_add, "ADD", {"source": source_xml, "schema": schema})

    if hold_for_two_years:
        today = datetime.datetime.now().date()
        two_years_from_today = datetime.date(today.year + 2, today.month, today.day)
        action_hold = ET.SubElement(actions, "ACTION")
        hold = ET.SubElement(
            action_hold, "HOLD", {"HoldUntilDate": two_years_from_today.isoformat()}
        )

    element_tree_to_file(submission, outfile)


def make_project_xml(outfile, alias, center_name, title_text, description_text):
    project_set = ET.Element("PROJECT_SET")
    project = ET.SubElement(
        project_set, "PROJECT", {"alias": alias, "center_name": center_name}
    )
    title = ET.SubElement(project, "TITLE")
    title.text = title_text
    description = ET.SubElement(project, "DESCRIPTION")
    description.text = description_text
    submission_project = ET.SubElement(project, "SUBMISSION_PROJECT")
    sequencing_project = ET.SubElement(submission_project, "SEQUENCING_PROJECT")
    element_tree_to_file(project_set, outfile)


def make_sample_xml(
    outfile, sample_alias, center_name, title, taxon_id, sample_attributes_dict=None
):
    sample_set = ET.Element("SAMPLE_SET")
    sample = ET.SubElement(
        sample_set, "SAMPLE", {"alias": sample_alias, "center_name": center_name}
    )
    title_element = ET.SubElement(sample, "TITLE")
    title_element.text = title
    sample_name = ET.SubElement(sample, "SAMPLE_NAME")
    taxon_id_element = ET.SubElement(sample_name, "TAXON_ID")
    taxon_id_element.text = str(taxon_id)

    if sample_attributes_dict is not None:
        sample_attributes_element = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")

        for tag, value in sorted(sample_attributes_dict.items()):
            sample_attribute = ET.SubElement(
                sample_attributes_element, "SAMPLE_ATTRIBUTE"
            )
            tag_element = ET.SubElement(sample_attribute, "TAG")
            tag_element.text = tag
            value_element = ET.SubElement(sample_attribute, "VALUE")
            value_element.text = value

    element_tree_to_file(sample_set, outfile)


def make_experiment_xml(
    outfile,
    experiment_alias,
    center_name,
    title,
    study_accession,
    sample_accession,
    library_name,
    platform,
    instrument,
    experiment_attributes_dict=None,
):
    experiment_set = ET.Element("EXPERIMENT_SET")
    experiment = ET.SubElement(
        experiment_set,
        "EXPERIMENT",
        {"alias": experiment_alias, "center_name": center_name},
    )
    title_element = ET.SubElement(experiment, "TITLE")
    title_element.text = title
    study_ref = ET.SubElement(experiment, "STUDY_REF", {"accession": study_accession})
    design = ET.SubElement(experiment, "DESIGN")
    design_description = ET.SubElement(design, "DESIGN_DESCRIPTION")
    sample_descriptor = ET.SubElement(
        design, "SAMPLE_DESCRIPTOR", {"accession": sample_accession}
    )
    library_descriptor = ET.SubElement(design, "LIBRARY_DESCRIPTOR")
    library_name_element = ET.SubElement(library_descriptor, "LIBRARY_NAME")
    library_name_element.text = library_name
    library_strategy = ET.SubElement(library_descriptor, "LIBRARY_STRATEGY")
    library_strategy.text = "WGS"
    library_source = ET.SubElement(library_descriptor, "LIBRARY_SOURCE")
    library_source.text = "GENOMIC"
    library_selection = ET.SubElement(library_descriptor, "LIBRARY_SELECTION")
    library_selection.text = "RANDOM"
    library_layout = ET.SubElement(library_descriptor, "LIBRARY_LAYOUT")
    ET.SubElement(library_layout, "PAIRED")
    platform_element = ET.SubElement(experiment, "PLATFORM")
    manufacturer_element = ET.SubElement(platform_element, platform)
    instrument_element = ET.SubElement(manufacturer_element, "INSTRUMENT_MODEL")
    instrument_element.text = instrument

    if experiment_attributes_dict is not None:
        experiment_attributes_element = ET.SubElement(
            experiment, "EXPERIMENT_ATTRIBUTES"
        )

        for tag, value in sorted(experiment_attributes_dict.items()):
            experiment_attribute = ET.SubElement(
                experiment_attributes_element, "EXPERIMENT_ATTRIBUTE"
            )
            tag_element = ET.SubElement(experiment_attribute, "TAG")
            tag_element.text = tag
            value_element = ET.SubElement(experiment_attribute, "VALUE")
            value_element.text = value

    element_tree_to_file(experiment_set, outfile)


def make_paired_fastq_run_xml(
    outfile, run_alias, center_name, experiment_accession, file_1, md5_1, file_2, md5_2
):
    run_set = ET.Element("RUN_SET")
    run = ET.SubElement(
        run_set, "RUN", {"alias": run_alias, "center_name": center_name}
    )
    experiment_ref = ET.SubElement(
        run, "EXPERIMENT_REF", {"accession": experiment_accession}
    )
    data_block = ET.SubElement(run, "DATA_BLOCK")
    files = ET.SubElement(data_block, "FILES")
    file_1_data = {
        "filename": file_1,
        "filetype": "fastq",
        "checksum_method": "MD5",
        "checksum": md5_1,
    }
    file_2_data = {
        "filename": file_2,
        "filetype": "fastq",
        "checksum_method": "MD5",
        "checksum": md5_2,
    }
    file_1_element = ET.SubElement(files, "FILE", file_1_data)
    file_2_element = ET.SubElement(files, "FILE", file_2_data)
    element_tree_to_file(run_set, outfile)
