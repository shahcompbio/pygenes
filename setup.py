import os
import platform
import versioneer

from setuptools import setup, find_packages, Extension


if not os.path.exists('pygenes/pygenes.cpp'):
    from Cython.Build import cythonize

    extensions = [
        Extension(
            'pygenes',
            ['pygenes/pygenes.pyx'],
            language='c++',
            include_dirs=['src'],
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

