from setuptools import setup, find_packages
#import versioneer
import glob

NAME = 'pointing'
CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
Programming Language :: Python
Natural Language :: English
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Astronomy
Operating System :: MacOS
Operating System :: POSIX
License :: OSI Approved :: MIT License
"""
URL = 'https://github.com/kadrlica/%s'%NAME
DESC = "Plot telescope pointings"
LONG_DESC = "See %s"%URL

setup(
    name=NAME,
    #version=versioneer.get_version(),
    #cmdclass=versioneer.get_cmdclass(),
    version="2.0.0",
    url=URL,
    author='Alex Drlica-Wagner',
    author_email='kadrlica@fnal.gov',
    scripts = ['bin/pointing','bin/pointing.py'],
    install_requires=[
        'numpy >= 1.7',
        'matplotlib >= 1.2.0',
        'basemap >= 1.0.6',
        'setuptools',
    ],
    packages=['pointing'],
    package_data={'pointing':['data/*.txt']},
    description=DESC,
    long_description=LONG_DESC,
    platforms='any',
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f]
)
