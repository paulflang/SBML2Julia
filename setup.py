import os
import setuptools


# extract version
with open(os.path.join(os.path.dirname(__file__),
          "sbml2julia", "_version.py")) as f:
    version = f.read().split('\n')[0].split('=')[-1].strip(' ').strip('"')

name = 'sbml2julia'
with open("README.md", encoding='utf-8') as fh:
    long_description = fh.read()

# install package
setuptools.setup(
    name=name,
    version=version,
    description='Optimization tool based on ODE discretisation.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/paulflang/" + name,
    download_url='https://github.com/paulflang/' + name,
    author="paulflang",
    author_email="paul.lang@wolfson.ox.ac.uk",
    license='MIT',
    keywords='Parameter fitting, ODE discretization, Julia, JuMP',
    packages=setuptools.find_packages(exclude=['tests', 'tests.*']),
    install_requires=['cement >= 3.0.0',
                      'julia',
                      'matplotlib',
                      'numpy',
                      'openpyxl',
                      'pandas',
                      'petab',
                      'python-libsbml',
                      'scipy',
                      'setuptools',
                      ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts': [
            'sbml2julia = sbml2julia.__main__:main',
        ],
    },
)
