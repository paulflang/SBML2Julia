import setuptools
import os

name = 'DisFit'
with open("README.md", "r") as fh:
    long_description = fh.read()
# dirname = os.path.dirname(__file__)
# package_data = {
#     name: [
#         'VERSION',
#     ],
# }

# # get package metadata
# md = pkg_utils.get_package_metadata(dirname, name, package_data_filename_patterns=package_data) # md is a PackageMetadata object that has attributes like version (as extracted from the VERSION file).
# print(20)
# print(dir(md))
# print(md.version)

# install package
setuptools.setup(
    name=name,
    version='0.0.1',
    description='Optimization tool based on ODE discretisation.',
    long_description=long_description,
    url="https://github.com/paulflang/" + name,
    download_url='https://github.com/paulflang/' + name,
    author="paulflang",
    author_email="paul.lang@wolfson.ox.ac.uk",
    license='MIT',
    keywords='Parameter fitting, ODE discretization',
    packages=setuptools.find_packages(exclude=['tests', 'tests.*']), # include all packages (i.e. folders, wit __inti__.py files) except tests
    install_requires=['cement >= 3.0.0', 'importlib', 'python-libsbml', 'matplotlib', 'numpy',
        # 'os',
        'pandas', # 're',
        'scipy', 'setuptools', # 'sys', 'tempfile',
        'julia'], # installs dependencies that are on PyPI
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
            'DisFit = DisFit.__main__:main', 
        ],
    }, # The entry_point says that when I type into the console "DisFit", what will be executed is the function main in calc.__main__.py
)
