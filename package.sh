#!/bin/bash
shopt -s extglob
rm pecar-anze.zip
zip -j pecar-anze.zip src/homework$1/!(_*) rpt/homework$1/pecar-anze.pdf
