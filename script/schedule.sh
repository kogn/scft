#!/bin/bash

pro_name=$1
task_num=$2
cycle_num=$3

i=1
while (( $i <= $cycle_num ))
do
  pro_num=$(ps -A | grep $pro_name |wc -l)

  if (( $pro_num < <ask_num )); then
    echo $i
    ./main 2 2 
    sleep 1s
    i=$(($i+1))
  else
    echo 'sleeping 600s'
    sleep 600s
  fi
done
