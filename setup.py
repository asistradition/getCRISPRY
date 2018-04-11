from setuptools import setup, find_packages

setup(
    name="getCRISPRY",
    version="0.1.0",
    license="MIT",
    description="CRISPR Target Design",
    author="Chris Jackson",
    packages=find_packages(exclude=["tests"]),
    install_requires=["biopython"]
)