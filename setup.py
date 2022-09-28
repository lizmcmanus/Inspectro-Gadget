"""
Setup script for scalp-distance tool.
"""

import os
from setuptools import setup, find_packages

# Get the current version number from inside the module
with open(os.path.join('inspectro_gadget', 'version.py')) as version_file:
    exec(version_file.read())

# Load the long description from the README
with open('README.md') as readme_file:
    long_description = readme_file.read()

# Load the required dependencies from the requirements file
with open("requirements.txt") as requirements_file:
    install_requires = requirements_file.read().splitlines()

setup(
    name='inspectro_gadget',
    version=__version__,
    description='Gives information about neurotransmitter receptor related mRNA expression within MRS voxels.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3.5',
    author='Liz McManus',
    author_email='liz.mcmanus93@googlemail.com ',
    url='https://github.com/lizmcmanus/Inspectro-Gadget',
    packages=find_packages(),
    license='MIT',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    platforms='any',
    keywords=['neuroscience', 'spectroscopy', 'MRI', 'mRNA', 'neurotransmission'],
    install_requires=install_requires,
    include_package_data=True,
    package_data={'': ['data/*.tsv', 'data/mRNA_images/*.nii.gz']},
)
