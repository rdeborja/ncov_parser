'''
Setup.py file for the ncov_parser package.
'''
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ncov_parser", # Replace with your own username
    version="0.4.0",
    author="Richard J. de Borja",
    author_email="richard.deborja@oicr.on.ca",
    description="A nCoV package for parsing analysis files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rdeborja/ncov_parser",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    scripts=['bin/get_qc_summary.py',
             'bin/collect_qc_summary.py',
             'bin/create_sample_qc_summary.py']
)
