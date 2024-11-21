from setuptools import setup, find_packages

setup(
    name='xraybinaryorbit',  # Your package name
    version='0.3.1',  # Version of your package
    license='MIT',  # License type
    description='A package for x-ray binary orbital modulations.',  # Short description
    author='LAEX',  # Your name
    author_email='graciela.sanjurjo@ua.es',  # Your email address
    url='https://github.com/xragua/xraybinaryorbit/',  # URL of your GitHub repository
    download_url='https://github.com/xragua/xraybinaryorbit/archive/refs/tags/one.tar.gz',  # Download URL for the package
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
