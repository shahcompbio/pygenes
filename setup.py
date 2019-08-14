import os
import sys
import platform
import versioneer

from setuptools import setup, find_packages, Extension


libraries = []
extra_compile_args = ['-g']
extra_link_args = ['-g']
if 'linux' in sys.platform:
    libraries.append('rt')
elif sys.platform == 'darwin':
    extra_compile_args.extend(['-stdlib=libc++'])
    extra_link_args.extend(['-stdlib=libc++', '-mmacosx-version-min=10.9'])

if not os.path.exists('pygenes/pygenes.cpp'):
    from Cython.Build import cythonize

    extensions = [
        Extension(
            'pygenes',
            ['pygenes/pygenes.pyx'],
            language='c++',
            include_dirs=['src'],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
        ),
    ]

    extensions = cythonize(extensions)

else:
    extensions = [
        Extension(
            'pygenes',
            ['pygenes/pygenes.cpp'],
            language='c++',
            include_dirs=['src'],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
        ),
    ]

setup(
    name='pygenes',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='In memory gene database searchable using an interval tree',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    keywords=['bioinformatics'],
    classifiers=[],
    ext_modules=extensions,
)

