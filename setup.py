# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='DyS',
    version='0.1',
    description='Dynamics simulation software of rigid multibody systems',
    long_description=readme,
    author='Luka Skrinjar, Ales Turel, Janko Slavic',
    author_email='skrinjar.luka@gmail.com',
    url='https://github.com/ladisk/DyS',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)



