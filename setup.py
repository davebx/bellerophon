import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bellerophon",
    author="Dave Bouvier",
    author_email="dave@bx.psu.edu",
    description="Filter reads of a minimum quality that span a junction, retaining the 5´ side of that junction",
    install_requires=['pysam'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/davebx/bellerophon",
    packages=['bellerophon'],
    entry_points={'console_scripts': ['bellerophon=bellerophon.cli:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: Academic Free License (AFL)',
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

