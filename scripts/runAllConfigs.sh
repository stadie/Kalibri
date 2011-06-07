#!/bin/bash

# Read config file names
#echo "%%% Enter directory name: "
#read DIR

# Loop over config files
CFG_DIR=(${1})
for i in `ls -1 ${CFG_DIR}/*.cfg`; do
  echo ""
  echo "%%% Executing file '${i}'"
  ./junk ${i}
  
  RES=(`basename ${i} .cfg`)
  echo "%%% Moving results to 'controlPlots/${RES}'"
  cd controlPlots
  mkdir ${RES}
  mv *.root ${RES}
  cp ../${i} ${RES}
  rm *.ps
  cd ..
  mv ${i} ${CFG_DIR}/done
  rm *.tex
done