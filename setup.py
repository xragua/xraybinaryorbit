from distutils.core import setup

setup(
    name='xraybinaryorbit',  # Your package name
    packages=['xraybinaryorbit'],  # The name of the package directory
    version='0.1',  # Initial version
    license='MIT',  # License type
    description='A package for x-ray binary orbital modulations.',  # Short description
    author='LAEX',  # Your name
    author_email='graciela.sanjurjo@ua.es',  # Your email address
    url='https://github.com/xragua/xraybinaryorbit/',  # URL of your GitHub repository
    download_url='https://github.com/xragua/xraybinaryorbit/archive/refs/tags/one.tar.gz',  # Download URL for the package
    keywords=['astronomy', 'x-ray, 'binary','orbit'],  # Keywords describing your package
    install_requires=[  # List of dependencies
        'numpy',
        'pandas',
        'scipy',
        'pyswarm',
        'matplotlib',
        'astropy',
    ],
    classifiers=[
        'Development Status :: v_01',  # Current development status
        'Intended Audience :: Astronomers',  # Audience type
        'Topic :: X-ray astronomy orbital modulations:: Analysis Tools',
        'License :: OSI Approved :: MIT License',  # License
        'Programming Language :: Python :: 3',  # Supported Python versions
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)

