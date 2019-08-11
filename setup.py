import platform
import versioneer


from distutils.core import setup
from distutils.extension import Extension

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
        ['src/pygenes.cpp'],
        language='c++',
        include_dirs=['src'],
        extra_compile_args=extra_compile_args),
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
    ext_modules=ext_modules,
)

