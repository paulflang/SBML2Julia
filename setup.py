import setuptools
try:
    import pkg_utils # a Karr Lab package for parsing metadata about the package
except ImportError:
    import pip._internal
    pip._internal.main(['install', '--process-dependency-links', 'git+https://github.com/KarrLab/pkg_utils.git#egg=pkg_utils'])
    import pkg_utils
import os

name = 'calc'
dirname = os.path.dirname(__file__)
package_data = {
    name: [
        'VERSION',
    ],
}

# get package metadata
md = pkg_utils.get_package_metadata(dirname, name, package_data_filename_patterns=package_data) # md is a PackageMetadata object that has attributes like version (as extracted from the VERSION file).
print(20)
print(dir(md))
print(md.version)

# install package
setuptools.setup(
    name=name,
    version=md.version,
    description='A simple calculator app',
    long_description=md.long_description,
    url="https://github.com/paulflang/" + name,
    download_url='https://github.com/paulflang/' + name,
    author="paulflang",
    author_email="paul.lang@wolfson.ox.ac.uk",
    license="MIT",
    keywords='calculator',
    packages=setuptools.find_packages(exclude=['tests', 'tests.*']), # include all packages (i.e. folders, wit __inti__.py files) except tests
    package_data=md.package_data, # includes also data files (as opposed to .py files, I guess). The input has the form {'dir1': [*.pattern1, *.pattern2]}
    install_requires=md.install_requires, # installs dependencies that are on PyPI
    extras_require=md.extras_require, # installs dependencies that are on PyPI
    tests_require=md.tests_require, # installs dependencies that are on PyPI
    dependency_links=md.dependency_links, # installs dependencies that are not on PyPI
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
            'calc = calc.__main__:main', 
        ],
    }, # The entry_point says that when I type into the console "calc", what will be executed is the function main in calc.__main__.py
)
