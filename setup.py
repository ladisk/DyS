# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='DyS',
    version='0.0.6',
    description='Simulation of rigid multibody system dynamics',
    long_description=readme,
    author='Kenneth Reitz',
    author_email='skrinjar.luka@gmail.com',
    url='https://github.com/jankoslavic/DyS',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)