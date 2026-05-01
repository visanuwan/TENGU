from setuptools import setup, find_packages

setup(
    name="TENGU",
    version="0.1.0",
    author="Visanu Wanchai",
    description="A transcriptomic-driven segmentation tool for Visium HD",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/visanuwan/TENGU",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "pandas",
        "numpy",
    ],
    entry_points={
        "console_scripts": [
            "tengu=tengu.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: Creative Commons :: Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)