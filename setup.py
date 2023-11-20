import setuptools
import sys

with open("README.md", "r") as fh:
    at_top = True
    long_description = ''
    for line in fh:
        if at_top and line[:3] == '[![':
            pass                # Skipping the badge-lines in the github README.md
        else:
            at_top = False      # Now starts the "real" README.md
        long_description += line


with open('dnctreek/dnctree/version.py') as fh:
    exec(fh.read())

if sys.version_info.major < 3:
    sys.exit('\n'
             'Sorry, Python 2 is not supported\n'
             'Did you run pip install dnctree?\n'
             'Try \'pip3 install dnctree\'')

elif sys.version_info.minor < 6:
    sys.exit('\nSorry, Python < 3.6 is not supported\n')

requirements = [
    'alv>=1.5',
    'modelmatcher>=1.1.3',
    'biopython>=1.70',
    'ete3>=3.1.1',
    'tree_matching_distance'
]

setuptools.setup(
    name="dnctreek",
    version=__version__,
    author="Lars Arvestad",
    author_email="arve@math.su.se",
    description="Distance-based phylogeny inference using a randomised divide-and-conquer method",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/arvestad/dnctree",
    test_suite = "tests",
    packages=setuptools.find_packages(),
    python_requires='>=3.5',    # I want to merge dictionaries easily. From 3.5, {**x, **y} merges dictionaries x and y
    #entry_points = {
    #    'console_scripts': ['dnctree = dnctree.main:main']
    #    },
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
