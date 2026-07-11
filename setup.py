from pathlib import Path

from setuptools import setup, find_packages

long_description = Path("README.md").read_text(encoding="utf-8")

setup(
    name='xraybinaryorbit',  # Your package name
    version='1.1.0',  # Version of your package
    license='MIT',  # License type
    description='A package for x-ray binary orbital modulations.',  # Short description
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='LAEX',  # Your name
    author_email='graciela.sanjurjo@ua.es',  # Your email address
    url='https://github.com/xragua/xraybinaryorbit/',  # URL of your GitHub repository
    download_url='https://github.com/xragua/xraybinaryorbit/archive/refs/tags/v1.1.0.tar.gz',  # Download URL for the package
    keywords=['astronomy', 'x-ray', 'binary', 'orbit'],  # Keywords describing your package

    # Dependencies
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'pyswarm',
        'matplotlib',
        'astropy',
    ],

    packages=find_packages(),  # Automatically find packages inside your directory

    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Minimum Python version requirement
)
