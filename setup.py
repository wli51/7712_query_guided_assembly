from setuptools import setup, find_packages

setup(
    name="query-guided-assembler",
    version="0.1.0",
    description="CPBS7712 course assignment for query-guided assembly and alignment of sequencing reads.",
    author="Weishan Li",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    install_requires=[
        "numpy==2.2.3",
        "pandas==2.2.3",
        "biopython==1.85",
        "tqdm==4.67.1"
    ],
    entry_points={
        "console_scripts": [
            "query-guided-assembly = query_guided_assembler.cli.run:main",
        ]
    },
    python_requires=">=3.11,<3.12",
)