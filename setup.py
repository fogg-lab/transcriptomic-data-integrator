from setuptools import setup, find_packages

PACKAGE_NAME = "transcriptomic_data_integrator"

setup(
    name=PACKAGE_NAME,
    version="0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "biopython",
        "GEOparse",
        "pandas",
        "requests",
        "mygene"
    ],
    package_data={
        PACKAGE_NAME: ["rscripts/*.R"],
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            f"configure-ncbi-email={PACKAGE_NAME}.update_ncbi_email:main"],
    }
)
