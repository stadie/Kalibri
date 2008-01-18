#!/bin/bash
for file in $( ls $1/* )
do
  echo 'ln -s '$file
  ln -s $file
done
