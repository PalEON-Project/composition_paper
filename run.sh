#!/bin/bash

# code to run the analyses used in the paper, including the final data product

git clone https://github.com/PalEON-Project/composition

cd composition

# create final data product
cp config_0.4-0 config
./master.sh  # probably best to run individual components of this file and not expect full script to run without issues (plus the fitting takes two weeks)

# now do CV analyses
