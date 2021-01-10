#!/bin/bash

OUT=snps
G=false
F=false
V=false

while getopts "hr:b:o:ngva" option; do
  case ${option} in
  h) echo ""
     echo "Usage: ";;
  r) REF=${OPTARG};;
  b) BAM=${OPTARG};;
  o) OUT=${OPTARG};;
  n) NORM=true;;
  g) G=true;;
  f) F=true;;
  v) V=true;;
  a) G=true
     F=true
     V=true;;
  esac
done 
shift "$(($OPTIND -1))"

if [[ -f ${REF} ]]
then
  echo "Found ${REF}"
else
  echo "Reference sequence not found"
  exit 0
fi

if [[ -f ${BAM} ]]
then
  echo "Found ${BAM}"
else
  echo "BAM file not found"
  exit 0
fi

if [ "$g" -eq true ] && [ "$v" -eq true ] && [ "$f" -eq true ]
then
  ./runVCs.sh -r ${REF} -b ${BAM} -o ${OUT} 
elif [ "$g" -eq true ] && [ "$v" -eq true ]
then
  ./runVCs.sh -r ${REF} -b ${BAM} -o ${OUT} -f
elif [ "$g" -eq true ] && [ "$f" -eq true ]
then
  ./runVCs.sh -r ${REF} -b ${BAM} -o ${OUT} -v
elif [ "$v" -eq true ] && [ "$f" -eq true ]
then
  ./runVCs.sh -r ${REF} -b ${BAM} -o ${OUT} -g
elif [ "$v" -eq true ]
then
  ./runVCs.sh -r ${REF} -b ${BAM} -o ${OUT} -g -f
elif [ "$f" -eq true ]
then
  ./runVCs.sh -r ${REF} -b ${BAM} -o ${OUT} -g -v
elif [ "$g" -eq true ]
then
  ./runVCs.sh -r ${REF} -b ${BAM} -o ${OUT} -f -v
fi

if [ "$norm" -eq true ]
then
  DICT=$(echo "${REF%.*}").dict
  ./normAndCombine.sh -r ${REF} -s ${DICT}
fi


