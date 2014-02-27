import fnmatch
import os
import sys
import argparse

from distutils.core import setup
from distutils.extension import Extension

argparser = argparse.ArgumentParser(add_help=False)
argparser.add_argument('--boost_source', help='boost source directory', required=True)
args, unknown = argparser.parse_known_args()
sys.argv = [sys.argv[0]] + unknown

def get_filenames(directory, filter):
    directory = os.path.abspath(os.path.expanduser(directory))
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, filter):
            yield os.path.join(root, filename)

boost_python_source = list(get_filenames(args.boost_source + '/libs/python/src', '*.cpp'))
boost_serialization_source = list(get_filenames(args.boost_source + '/libs/serialization/src', '*.cpp'))

ext_modules = [Extension('pygenes',
                         ['src/pygenes.cpp'] + boost_python_source + boost_serialization_source, 
                         language='c++',
                         include_dirs=['src', args.boost_source],
                         extra_compile_args=['-ftemplate-depth-1024'])]

setup(name='pygenes', ext_modules=ext_modules)

