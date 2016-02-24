from distutils.core import setup
from distutils.extension import Extension

ext_modules = [
    Extension(
        'pygenes',
        ['src/pygenes.cpp'],
        language='c++',
        include_dirs=['src'],
        libraries=['boost_python', 'boost_serialization'],
    )
]

setup(
    name='pygenes',
    version='0.0.1',
    description='On-disk gene database searchable using an interval tree',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    url='http://bitbucket.org/dranew/pygenes',
    download_url='https://bitbucket.org/dranew/pygenes/get/v0.0.1.tar.gz',
    keywords=['bioinformatics'],
    classifiers=[],
    ext_modules=ext_modules,
)

