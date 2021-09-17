#!/bin/bash
## Script that downloads
## some source data and iPython notebooks 
## into a Jupyter image for Python and Machine Learning course

cd /home
# git reflog requires a name and email if user is not in passwd
# even if you're only cloning
wget https://raw.githubusercontent.com/narges-s/rna-seq/main/counttab.csv
wget https://raw.githubusercontent.com/narges-s/rna-seq/main/annottab.cvs
export GIT_COMMITTER_NAME=narges-s
export GIT_COMMITTER_EMAIL=soren.narges@gmail.com
git clone https://github.com/narges-s/rna-seq.git
