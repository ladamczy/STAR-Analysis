#!/bin/bash

#creating output directory
current_date=$(date -Iseconds)
output_directory=/star/u/adamwatroba/STAR-Analysis/Inclusive_Analysis/starsim/starsim_and_bfc/results/$current_date
mkdir -p $output_directory
echo Output directory:
echo $output_directory

# setting default values (with checking whether they are actual arguments and not -* ones)
#if we find a -* one, we just skip checking the rest of them
length=0
run_number=0
seed=0
#1
if [[ "$1" != -* ]]; then
    length=$1
    # length="${length:-0}"
    OPTIND=2
fi
#2
if [ "$OPTIND" -lt "2" ]; then
    #do nothing because that means the -* arguments started earlier
    :
elif [[ "$2" != -* ]]; then
    run_number=$2
    # run_number="${run_number:-0}"
    OPTIND=3
fi
#3
if [ "$OPTIND" -lt "3" ]; then
    #do nothing because that means the -* arguments started earlier
    :
elif [[ "$3" != -* ]]; then
    seed=$3
    # seed="${seed:-0}"
    OPTIND=4
fi

# setting values based on additional parameters
extension="root"
filter=""
NO_RUNNING=0
while getopts "e:f:t" flag
do
    case $flag in
        #case for checking all extensions i'd like to copy
        e)  if [[ "$OPTARG" == "Mu" ]]; then
                extension="MuDst.root"
            fi
            ;;
        #case for checking active filters
        f)  case $OPTARG in
                "K0S")
                    filter="K0S"
                    ;;
                "Lambda0")
                    filter="Lambda0"
                    ;;
                "Kstar")
                    filter="Kstar"
                    ;;
                "phi")
                    filter="Phi"
                    ;;
                *)
                    ;;
            esac
            ;;
        #case for testing the arguments without running
        t)  NO_RUNNING=1
            echo -e "\nTHIS IS A TEST; IT WILL *NOT* BE SCHEDULED\n"
            ;;
        #default case
        *)  ;;
    esac
done

#writing down the choices
echo "You chose the following arguments:"
echo -e "length:\t\t$length"
echo -e "run_number:\t$run_number"
echo -e "seed:\t\t$seed"
echo -e "extension:\t$extension"
if [[ -z "$filter" ]]; then
    echo -e "filter:\t\tstrange particles (K0S, Lambda0, K*, phi)"
else
    echo -e "filter:\t\t$filter"
fi

#show help if there is 0 arguments provided
if [ "$#" -eq "0" ]; then
    echo Usage:
    echo ./centralDiffractive_total_chain_start.sh [number of events] [run number \(optional, default 18091010\)] ["seed" \(optional, default 0\)] [optional arguments]
    echo
    echo Optional arguments:
    echo
    echo -e "-e: Mu \t\t\t copies .MuDst.root files only"
    echo -e "    anything else:\t copies .root files"
    echo
    echo -e "-f: K0S \t\t applies filter for K0S particles only"
    echo -e "    Lambda0 \t\t applies filter for Lambda0 particles only"
    echo -e "    Kstar \t\t applies filter for K*(892) particles only"
    echo -e "    phi \t\t applies filter for phi(1020) particles only"
    echo -e "    anything else:\t applies filter for K0S, Lambda0, K*(892) and phi(1020) particles only"
fi

# finishing if the first argument is not a number (through the power of regex)
# or if there is none (0) events
# or if -t option is specified
re='^[0-9]+$'
if ! [[ $length =~ $re ]] ; then
    exit 1
fi
if [ "$length" -eq "0" ]; then
    exit 1
fi
if [ "$NO_RUNNING" -eq "1" ]; then
    exit 1
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
        star-submit-template -template centralDiffractive_total_chain_template.xml -entities number=$i,events_number=$n,run_number=$run_number,seed=$seed,output_directory=$output_directory,extension=$extension,filter=$filter
    done
fi

# doing last, not-full batch
if [ "$partial_events" -ne "0" ]; then
    star-submit-template -template centralDiffractive_total_chain_template.xml -entities number=$full_events,events_number=$partial_events,run_number=$run_number,seed=$seed,output_directory=$output_directory,extension=$extension,filter=$filter
fi

mv strange_generator*.* ./star_scheduler_logs/
rm schedTemplateExp.xml
