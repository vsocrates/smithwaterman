import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="naive-smith-waterman-vsocrates", 
    version="0.0.1",
    author="Vimig Socrates",
    author_email="vimig.socrates@yale.edu",
    description="A Naive Implementation of Smith-Waterman Local Alignment",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vsocrates/smithwaterman",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'pandas'
    ],    
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
)