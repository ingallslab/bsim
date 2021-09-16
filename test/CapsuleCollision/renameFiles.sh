#!/bin/bash

#This script replaces "FIND" with "NEW" in all the file names and inside each file
#requires rename package

FIND="CrossFeeding"
NEW="CapsuleCollision"

rename s/${FIND}/${NEW}/ *.java
sed -i 's/${FIND}/${NEW}/g' *.java
