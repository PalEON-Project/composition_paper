#!/bin/bash

# code to run the analyses used in the paper, including the final data product

git clone https://github.com/PalEON-Project/composition

cd composition

# create final data product
cp config_0.3-0 config
./master.sh

# now do CV analyses