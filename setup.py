import platform
import versioneer

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

extra_compile_args = ['-Wno-unused-variable']
if platform.platform().startswith('Darwin'):
    extra_compile_args.append('-Wno-unneeded-internal-declaration')
    extra_compile_args.append('-Wno-unused-private-field')
    extra_compile_args.append('-stdlib=libc++')
    extra_compile_args.append('-mmacosx-version-min=10.9')

extra_compile_args.append('-ftemplate-depth-1024')

extra_link_args = []
if platform.platform().startswith('Darwin'):
    extra_link_args.append('-Wl,-no_compact_unwind')
    extra_link_args.append('-stdlib=libc++')
    extra_link_args.append('-mmacosx-version-min=10.9')

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

setup(
    name='pygenes',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='In memory gene database searchable using an interval tree',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    keywords=['bioinformatics'],
    classifiers=[],
    ext_modules=cythonize(extensions),
)

