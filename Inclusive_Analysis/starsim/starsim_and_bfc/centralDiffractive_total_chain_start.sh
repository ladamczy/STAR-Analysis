#!/bin/bash

#creating output directory
current_date=$(date -Iseconds)
output_directory=/star/u/adamwatroba/STAR-Analysis/Inclusive_Analysis/starsim/starsim_and_bfc/results/$current_date
mkdir -p $output_directory
echo Output directory:
echo $output_directory

# setting default values
length=$1
length="${length:-0}"
run_number=$2
run_number="${run_number:-0}"
seed=$3
seed="${seed:-0}"

#show help if there is 0 arguments provided
if [ "$#" -eq "0" ]; then
    echo Usage:
    echo ./centralDiffractive_total_chain_template.sh [number of events] [run number \(optional, default 18091010\)] ["seed" \(optional, default 0\)]
fi

# finishing if there is none (0) events
if [ "$length" -eq "0" ]; then
    exit 0
fi

#creating file with info about simulation
newfile=$output_directory/simulation_data.txt
echo Number of events: > $newfile
echo $length >> $newfile
echo Run number \(if 0 then actual one is 18091010\): >> $newfile
echo $run_number >> $newfile
echo Seed: >> $newfile
echo $seed >> $newfile

# modulo to cut it into pieces of n events
n=100
full_events=$(($length / $n))
partial_events=$(($length % $n))

# doing loop of batches of full-sized groups of events
if [ "$full_events" -ne "0" ]; then
    for i in $( eval echo {0..$(($full_events-1))} );
    do
        star-submit-template -template centralDiffractive_total_chain_template.xml -entities number=$i,events_number=$n,run_number=$run_number,seed=$seed,output_directory=$output_directory
    done
fi

# doing last, not-full batch
if [ "$partial_events" -ne "0" ]; then
    star-submit-template -template centralDiffractive_total_chain_template.xml -entities number=$full_events,events_number=$partial_events,run_number=$run_number,seed=$seed,output_directory=$output_directory
fi

mv strange_generator*.* ./star_scheduler_logs/
rm schedTemplateExp.xml
