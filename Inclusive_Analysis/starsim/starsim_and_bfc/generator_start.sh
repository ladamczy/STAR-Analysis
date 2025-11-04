#!/bin/bash

for i in {1..1000};
do
    star-submit-template -template generator_template.xml -entities number=$i
done
