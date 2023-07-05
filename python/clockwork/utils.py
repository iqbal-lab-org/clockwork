import hashlib
import datetime
import logging
import os
import subprocess
import sys

import pyfastaq


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def syscall(command):
    logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logging.info(f"Return code: {completed_process.returncode}")
    if completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "Output from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise Exception("Error in system call. Cannot continue")

    logging.info(f"stdout:\n{completed_process.stdout.rstrip()}")
    logging.info(f"stderr:\n{completed_process.stderr.rstrip()}")
    return completed_process


def md5(filename):
    """Given a file, returns a string that is the md5 sum of the file"""
    # see https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def load_md5_from_file(filename):
    """Loads md5 from file, where could have been made my 'md5' on a mac,
    or 'md5sum' in linux. Assumes just one md5 in the file - only looks at
    the first line"""
    with open(filename) as f:
        line = f.readline().rstrip()

    # Mac:
    # MD5 (filename) = md5
    # Linux:
    # md5  filename
    try:
        if line.startswith("MD5"):
            md5sum = line.split()[-1]
        else:
            md5sum = line.split()[0]
    except:
        raise Exception("Error getting md5 from file " + filename + ". Unexpected format")

    if len(md5sum) != 32:
        raise Exception(
            "Error getting md5 from file " + filename + ". Expected string of length 32"
        )

    return md5sum


def rsync_and_md5(old_name, new_name, md5sum=None):
    """Copies a file from old_name to new_name using rsync.
    Double-checks the copy was successful using md5. Returns md5.
    If you already know the md5 of the file, then save time
    by providing it with the md5sum option - this will avoid
    calculating it on the old file."""
    if md5sum is None:
        md5sum = md5(old_name)

    syscall("rsync " + old_name + " " + new_name)
    new_md5sum = md5(new_name)

    if new_md5sum != md5sum:
        raise Exception(
            "Error copying file "
            + old_name
            + " -> "
            + new_name
            + "\n. md5s do not match"
        )
    else:
        return md5sum


def date_string_from_file_mtime(filename):
    """Returns a string in the form YYYYMMDD of the last modification
    date of a file"""
    try:
        mtime = os.path.getmtime(filename)
    except:
        raise Exception("Error getting modification time from file " + filename)

    d = datetime.datetime.fromtimestamp(mtime)
    return d.isoformat().split("T")[0].replace("-", "")


def make_empty_file(filename):
    """Makes empty file. Will overwrite if already exists"""
    with open(filename, "w"):
        pass


def sam_record_count(filename):
    """Returns number of sam records in file"""
    count = 0
    with open(filename) as f:
        for line in f:
            if not line.startswith("@"):
                count += 1
    return count


def vcf_has_records(filename):
    """Returns true if there is at least 1 record in VCF file"""
    with open(filename) as f:
        for line in f:
            if not line.startswith("#"):
                return True
    return False


def file_has_at_least_one_line(filename):
    """Returns true if there is at least 1 line in the file. Can be gzipped"""
    f = pyfastaq.utils.open_file_read(filename)
    has_line = any(f)
    f.close()
    return has_line
