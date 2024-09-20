from distutils.core import setup

setup(
    name='xraybinaryorbit',  # Your package name
    packages=['xraybinaryorbit'],  # The name of the package directory
    version='0.2.3',  # Initial version
    license='MIT',  # License type
    description='A package for x-ray binary orbital modulations.',  # Short description
    author='LAEX',  # Your name
    author_email='graciela.sanjurjo@ua.es',  # Your email address
    url='https://github.com/xragua/xraybinaryorbit/',  # URL of your GitHub repository
    download_url='https://github.com/xragua/xraybinaryorbit/archive/refs/tags/one.tar.gz',  # Download URL for the package
    keywords=['astronomy', 'x-ray', 'binary','orbit'],  # Keywords describing your package
    install_requires=[  # List of dependencies
        'numpy',
        'pandas',
        'scipy',
        'pyswarm',
        'matplotlib',
        'astropy',
    ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',  # Current development status

    ],
)


