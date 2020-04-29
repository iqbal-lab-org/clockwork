import glob
from setuptools import setup, find_packages


setup(
    name="clockwork",
    version="0.8.3",
    description="Pipeline code for CRyPTIC project",
    packages=find_packages(),
    author="Martin Hunt",
    author_email="mhunt@ebi.ac.uk",
    url="https://github.com/iqbal-lab-org/clockwork",
    scripts=glob.glob("scripts/*"),
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=[
        "cluster_vcf_records >= 0.11.1",
        "python-dateutil >= 2.6.1",
        "openpyxl >= 2.4.7",
        "pyfastaq >= 3.14.0",
        "pymysql >= 0.7.11",
        "pysam >= 0.11.2.1",
        "requests >= 2.9.1",
        "xlsxwriter >= 1.0.0",
    ],
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
