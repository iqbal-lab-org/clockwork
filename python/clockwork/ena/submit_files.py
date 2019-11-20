import configparser
import random
import string

from clockwork import utils


class Error(Exception):
    pass


def _make_dummy_success_receipt(outfile, object_type):
    accession = "".join(
        [random.choice(string.ascii_uppercase + string.digits) for _ in range(10)]
    )
    with open(outfile, "w") as f:
        print(
            r"""<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" href="receipt.xsl"?>
<RECEIPT receiptDate="2017-08-31T11:07:50.251+01:00" submissionFile="submission.xml" success="true">
    <"""
            + object_type.upper()
            + r''' accession="'''
            + accession
            + r"""" alias="alias" status="PRIVATE">
        <EXT_ID accession="SAMEA123456789" type="biosample"/>
    </"""
            + object_type.upper()
            + r""">
    <SUBMISSION accession="ERA1234567" alias="alias 42"/>
    <MESSAGES>
        <INFO>Submission has been committed.</INFO>
        <INFO>This submission is a TEST submission and will be discarded within 24 hours</INFO>
    </MESSAGES>
    <ACTIONS>ADD</ACTIONS>
</RECEIPT>
""",
            file=f,
        )


def _make_dummy_fail_receipt(outfile):
    with open(outfile, "w") as f:
        print(
            r"""<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet type="text/xsl" href="receipt.xsl"?>
<RECEIPT receiptDate="2017-09-01T14:13:19.573+01:00" submissionFile="submission.xml" success="false">
    <SUBMISSION alias="Submission alias run 1"/>
    <MESSAGES>
        <ERROR>In submission, alias:"Submission alias run 1", accession:"". The object being added already exists in the submission account.</ERROR>
        <INFO>Submission has been rolled back.</INFO>
        <INFO>This submission is a TEST submission and will be discarded within 24 hours</INFO>
    </MESSAGES>
    <ACTIONS>ADD</ACTIONS>
</RECEIPT>
""",
            file=f,
        )


def parse_config_file(ini_file):
    config = configparser.ConfigParser()
    try:
        config.read(ini_file)
    except:
        raise Error("Error! Parsing config file " + ini_file)

    if "ena_login" not in config:
        raise Error("Error! [ena_login] section not found in config file " + ini_file)

    if not ("user" in config["ena_login"] and "password" in config["ena_login"]):
        raise Error(
            "Error! user and password not found in [ena_login] section of config file "
            + ini_file
        )

    return config["ena_login"]["user"], config["ena_login"]["password"]


def submit_xml_files(
    ini_file,
    outfile,
    files=None,
    use_test_server=False,
    unit_test=None,
    unit_test_obj_type=None,
):
    username, password = parse_config_file(ini_file)
    if files is None:
        files_command = None
    else:
        files_command_list = [
            '-F "' + key + "=@" + value + '"' for key, value in files.items()
        ]
        files_command = " ".join(files_command_list)

    if use_test_server:
        url = "https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20"
    else:
        url = "https://www.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20"

    command_list = [
        "curl -k",
        files_command,
        '"' + url + username + "%20" + password + '"',
        ">",
        outfile,
    ]

    command = " ".join([x for x in command_list if x is not None])
    if unit_test is None:
        utils.syscall(command)
    elif unit_test == "success":
        _make_dummy_success_receipt(outfile, unit_test_obj_type)
    elif unit_test == "fail":
        _make_dummy_fail_receipt(outfile)
    else:
        raise Error("unit_test must be None, success, or fail. Got: " + unit_test)


def upload_file_to_ena_ftp(ini_file, filename, uploaded_name):
    # paranoid about passwords and running ps? Looks like curl is ok:
    # https://unix.stackexchange.com/questions/385339/how-does-curl-protect-a-password-from-appearing-in-ps-output
    # "wipe the next argument out so that the username:password isn't
    # displayed in the system process list"
    username, password = parse_config_file(ini_file)
    cmd = " ".join(
        [
            "curl -T",
            filename,
            "ftp://webin.ebi.ac.uk/" + uploaded_name,
            "--user",
            username + ":" + password,
        ]
    )
    utils.syscall(cmd)
