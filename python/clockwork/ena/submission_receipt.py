import xml.etree.ElementTree as ET


class Error(Exception):
    pass


class SubmissionReceipt:
    def __init__(self, filename):
        self.filename = filename
        self.tree = ET.parse(self.filename)
        self.root = self.tree.getroot()
        if self.root.tag != "RECEIPT":
            raise Error(
                'Error! Expected first tag to be "RECEIPT" in XML file "'
                + self.filename
                + '"'
            )
        self.successful = (
            "success" in self.root.attrib and self.root.attrib["success"] == "true"
        )
        self.accessions = {
            child.tag: child.attrib["accession"]
            for child in self.root
            if "accession" in child.attrib
        }
