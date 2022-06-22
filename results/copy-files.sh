#!/bin/sh
mkdir /Users/marlena/Documents/FileZilla/xs-recent/enhanced/$1

scp -r mnorwood@bayes.biostat.washington.edu:/home/students/mnorwood/hutch/xs-recent/enhanced/$1 /Users/marlena/Documents/FileZilla/xs-recent/enhanced/
