#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Usage: ./script.sh argument1 argument2"
    exit 1
fi
arg1="$1"
arg2="$2"
arg3="solved_alarm.bif"
./strt "$arg1" "$arg2" "$arg3"