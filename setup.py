from setuptools import setup, find_packages

setup(
    name='mercat',
    version='0.1',
    author='Ajay Panyala',
    author_email='ajay@example.com',
    url='https://github.com/pnnl/mercat',
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'mercat=mercat.mercat:main',
        ],
    }
)
