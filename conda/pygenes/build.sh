#!/bin/bash

wget https://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
tar -xzvf boost_1_60_0.tar.gz
$PYTHON setup.py install --boost_source boost_1_60_0
