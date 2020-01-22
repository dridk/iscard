import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="iscard",
    version="0.0.1",
    author="Sacha Schutz",
    author_email="sacha@labsquare.org",
    description="CNV detection from Amplicon sequencing",
    long_description="Iscard use a supervised learning machin to detect CNV from sequencing amplicon strategy",
    long_description_content_type="text/markdown",
    url="https://github.com/dridk/iscard",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: AGPL License",
        "Operating System :: OS Independent",
    ),
)