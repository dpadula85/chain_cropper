#!/usr/bin/env python3
"""
Setup script for ChainCropper package.
"""

from setuptools import setup, find_packages

setup(
    name='chain-cropper',
    version=0.1,
    author="Daniele Padula",
    author_email="daniele.padula@unisi.it",
    description='A tool for cropping molecular side chains from structures and trajectories',
    long_description=open("README.md").read(),
    long_description_content_type='text/markdown',
    url='https://github.com/dpadula85/chain_cropper',
    packages=find_packages(),
    license='GPL-3.0-or-later',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    install_requires=[
        'MDAnalysis>=2.0.0',
        'numpy>=1.19.0',
        'tqdm',
        'joblib',
    ],
    entry_points={
        'console_scripts': [
            'chain-cropper=chain_cropper.cli:main',
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords=[
        "molecular dynamics",
        "computational chemistry",
        "chemistry",
        "materials science",
    ],
)
