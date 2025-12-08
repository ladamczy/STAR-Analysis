#!/bin/bash

rm -rf send.package

# setting default value so that only one batch is made in default case
length=$1
length="${length:-0}"

for i in $( eval echo {0..$length} );
do
    star-submit-template -template bfc_template.xml -entities number=$i
done

mv bfcchain*.* ./star_scheduler_logs/
rm schedTemplateExp.xml
rm send.zip
