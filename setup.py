import fnmatch
import os
import sys
import platform
import argparse

from distutils.core import setup
from distutils.extension import Extension

argparser = argparse.ArgumentParser(add_help=False)
argparser.add_argument('--boost_source', help='boost source directory', required=True)
args, unknown = argparser.parse_known_args()
sys.argv = [sys.argv[0]] + unknown

boost_source = os.path.expanduser(args.boost_source)

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

boost_python_sources = list(find_files(os.path.join(boost_source, 'libs/python/src'), '*.cpp'))
boost_serialization_sources = list(find_files(os.path.join(boost_source, 'libs/serialization/src'), '*.cpp'))

assert(len(boost_python_sources) > 0)
assert(len(boost_serialization_sources) > 0)

extra_compile_args = ['-Wno-unused-variable']
if platform.platform().startswith('Darwin'):
    extra_compile_args.append('-Wno-unneeded-internal-declaration')
    extra_compile_args.append('-Wno-unused-private-field')

extra_compile_args.append('-ftemplate-depth-1024')

extra_link_args = []
if platform.platform().startswith('Darwin'):
    extra_link_args.append('-Wl,-no_compact_unwind')

ext_modules = [
    Extension(
        'pygenes',
        ['src/pygenes.cpp'] + boost_python_sources + boost_serialization_sources,
        language='c++',
        include_dirs=['src', boost_source],
        extra_compile_args=extra_compile_args),
]

setup(
    name='pygenes',
    version='0.0.3',
    description='In memory gene database searchable using an interval tree',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    url='http://bitbucket.org/dranew/pygenes',
    download_url='https://bitbucket.org/dranew/pygenes/get/v0.0.3.tar.gz',
    keywords=['bioinformatics'],
    classifiers=[],
    ext_modules=ext_modules,
)

