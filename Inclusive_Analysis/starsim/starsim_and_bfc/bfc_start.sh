#!/bin/bash

#creating output directory
output_directory=/star/data05/scratch/adamwatroba/bfc_output
mkdir -p $output_directory

# setting default values
input=$1
input="${input:-0}"
run_number=$2
run_number="${run_number:-0}"

#show help if there is 0 arguments provided
if [ "$#" -eq "0" ]; then
    echo Usage:
    echo ./bfc_start.sh [path to: .fzd file / .list file / \"expression with \* ending in .fzd\" \(in quotation marks\)] [run number \(optional, default 18091010\)]
fi

#creating filelist from asterisks
if [[ $input == *"\*"* ]]; then
  echo $input | tr " " "\n" | tee tempfilelist.list >/dev/null
#creating filelist from a single file
elif [[ $input == *".fzd" ]]; then
  echo $input > tempfilelist.list
#creating filelist from... a filelist
elif [[ $input == *".list" ]]; then
  cp $input tempfilelist.list
#error cause unknown file
else
  echo Wrong input file!
  exit 1
fi

#creating file with info about bfc
newfile=$output_directory/bfc_data.txt
echo List of files used: > $newfile
cat tempfilelist.list >> $newfile
echo Run number \(if 0 then actual one is 18091010\): >> $newfile
echo $run_number >> $newfile

#creating full path filelist because otherwise we get problems
filelist=$(pwd)
filelist=$filelist/tempfilelist.list
echo Path to tempfilelist.list:
echo $filelist

star-submit-template -template bfc_template.xml -entities filelist=$filelist,run_number=$run_number

mv bfcchain*.* ./star_scheduler_logs/
rm schedTemplateExp.xml
rm -f tempfilelist.list
