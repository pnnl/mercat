import io
from os.path import dirname, join
from setuptools import setup, find_packages


setup(
    name='mercat',
    version="0.1",
    author='Ajay Panyala',
    author_email='ajay@example.com',
    url='https://github.com/pnnl/mercat',
    packages=['mercat'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'mercat=mercat.mercat:mercat_main',
        ],
    }
)


