from setuptools import setup, find_packages

PACKAGE_NAME = "transcriptomics_data_query"

setup(
    name=PACKAGE_NAME,
    version="0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "biopython",
        "GEOparse"
    ],
    package_data={
        "mypackage": ["email_for_ncbi_tracking.txt"],
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            f"configure-ncbi-email={PACKAGE_NAME}.update_ncbi_email:main"],
    }
)
