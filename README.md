# A Naive Smith-Waterman Implementation

This implementation of Smith-Waterman calculates performs local alignment and outputs a file with the sequences, scoring matrix, alignment score, and a human-readable format of the alignment. Examples of use are below. 

## Dependencies

There are only two dependencies that are required to use this package. They are listed below: 

```
numpy==1.16.6
pandas==0.24.2
```

The system uses **Python 2.7**. I know, I know, we've retired it but bad habits are hard to break :) 

## Installation

This program has been packaged as a python module with `pip` and is really easy to install. Install using the command below. 

```
pip install git+https://github.com/vsocrates/smithwaterman.git
```

It should install all dependencies for you too! I'd create a [virtualenv](https://virtualenv.pypa.io/en/latest/) to keep everything clean. 

## How to Use

There are a few required inputs, and a couple more optional ones described below. In keeping with the skeleton code under which this program was developed, run the program as a module: 

```
python -m naivesw.hw1 -i input.txt -s score.txt -t output.txt
```

If we run the `python -m naivesw.hw1 -h` (help command) we get the following with other optional parameters for penalties.

```
usage: hw1.py [-h] -i INPUT -t OUTPUT -s SCORE [-o OPENGAP] [-e EXTGAP]

Smith-Waterman Algorithm

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file
  -t OUTPUT, --output OUTPUT
                        output file
  -s SCORE, --score SCORE
                        score file
  -o OPENGAP, --opengap OPENGAP
                        open gap
  -e EXTGAP, --extgap EXTGAP
                        extension gap
```

## CBB 752 Homework 1 Trial Command

Assuming our input, score, and output file locations are all in the current directory, after installation, we'd use

```
python -m naivesw.hw1 -i input.txt -s blosum62.txt -t output.txt
```
