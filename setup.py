from setuptools import setup, find_packages
'''
setup file
'''
setup(
    name = 'LJengo',
    version = 0.1,
    package = find_packages(),
    install_requires = ['numpy', 'matplotlib', 'tqdm'],
    author = 'Rahul Verma',
    author_email = 'rverma7@ncsu.edu',
    description = 'molecular dynamics simulation package for LJ model',
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
    ],
    python_requires = '>=3.0',
    license = 'MIT',
    keywords = 'LJ MD',
    url = 'https://github.com/umraorahul/Lennard_Jones_MD',
    docs_url='',
    packages = ['LJengo'],
)
