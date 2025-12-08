#!/bin/bash

# setting default values
length=$1
length="${length:-0}"
run_number=$2
run_number="${run_number:-0}"

# finishing if there is none (0) events
if [ "$length" -eq "0" ]; then
    exit 0
fi

# modulo to cut it into pieces of n events
n=100
full_events=$(($length / $n))
partial_events=$(($length % $n))

# doing loop of batches of full-sized groups of events
if [ "$full_events" -ne "0" ]; then
    for i in $( eval echo {0..$(($full_events-1))} );
    do
        star-submit-template -template bfc_template.xml -entities number=$i,run_number=$run_number
    done
fi

# doing last, not-full batch
if [ "$partial_events" -ne "0" ]; then
    star-submit-template -template bfc_template.xml -entities number=$full_events,run_number=$run_number
fi

mv bfcchain*.* ./star_scheduler_logs/
rm schedTemplateExp.xml
