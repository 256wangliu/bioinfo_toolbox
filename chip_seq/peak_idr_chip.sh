#!/bin/bash

# Copyright 2015 Junko Tsuji

# Curate detected peaks passed IDR threshould and make master peak set
# for a group.

#### Software path
export PATH=$PATH:/home/tsujij/soft/bedtools2-2.20.1/bin
IDR_TOOL_DIRECTORY="/home/tsujij/soft/idrCode"

#### Usage
function usage {
cat <<EOF
Usage: $0 [options] [-N | -B] -G <group>

Curate detected peaks passed IDR threshould and make master peak set for a group.

Required argument:
  -N           Narrow peak IDRs.
  -B           Broad peak IDRs.
  -G <group>   Directory containing narrowPeak or broadPeak files that belong to
               the same category, such as cancer or control.
Options:
  -h           Show this help message and exit.
  -I           IDR cutoff (default: 0.05 for narrow peaks and 0.1 for broad peaks).
  -o <prefix>  Prefix for output files (default: <group>).
EOF
}

#### Parameters
R=()
P=()
unset OFPREFIX
unset GROUP
unset PEAKOPT
unset IDR

#### Read arguments and options
[[ $# -eq 0 ]] && usage && exit 1;
while getopts "hG:o:I:NB" OPT
do
  case ${OPT} in
    "G") GROUP=${OPTARG}    ;;
    "o") OFPREFIX=${OPTARG} ;;
    "I") IDR=${OPTARG}      ;;
    "B") PEAKOPT="T"        ;;
    "N") PEAKOPT="F"        ;;
    "h") usage && exit 1    ;;
     * ) usage && exit 1    ;;
  esac
done
[[ $# -lt 3 ]] && echo "$0: need more arguments" && exit 1;

#### Check arguments
if [ ! -d "${GROUP}" ]
  then echo "$0: can't read ${GROUP}" && exit 1; fi
if [ -z "${IDR}" ]
then
  if [ "${PEAKOPT}" == "T" ]
    then IDR="0.1"; else IDR="0.05"; fi
else
  if [ $( echo ${IDR} | sed 's/^[+0-9][0-9]*[.][0-9]*//' | wc -c ) -ne 1 ]
    then echo "$0: can't interpret -I ${IDR}" && exit 1; fi
fi
if [ -z "${OFPREFIX}" ]
  then OFPREFIX=`basename ${GROUP}`; fi

#### Compute IDRs for all pairs of samples in a group
CURRENT_DIRECTORY=`pwd`
GROUP=`readlink -f "${GROUP}"`
cd ${IDR_TOOL_DIRECTORY}
for F in `ls ${GROUP}`; do R=("${R[@]}" "${GROUP}/${F}"); done
for (( I=0; I < ${#R[@]}-1; I++ ))
do
  for (( J=I+1; J < ${#R[@]}; J++ ))
  do
    R1=${R[$I]}
    R2=${R[$J]}
    R1_PREFIX=`basename ${R1} | awk -F "." '{print $1}'`
    R2_PREFIX=`basename ${R2} | awk -F "." '{print $1}'`
    R1_VS_R2="${GROUP}/${R1_PREFIX}_VS_${R2_PREFIX}"
    P=("${P[@]}", "${R1_VS_R2}")

    Rscript batch-consistency-analysis.r ${R1} ${R2} -1 ${R1_VS_R2} 0 ${PEAKOPT} p.value

    #### Extract peaks passed IDR threshould: chr, start, end, p.value
    sed 1d "${R1_VS_R2}-overlapped-peaks.txt" | sed -r 's/"//g' | \
    awk -v co=${IDR} 'BEGIN{OFS="\t"}{ if($11 <= co){print $2,$3,$4,$5} }' >> "${GROUP}/${R1_PREFIX}.pool"
    sed 1d "${R1_VS_R2}-overlapped-peaks.txt" | sed -r 's/"//g' | \
    awk -v co=${IDR} 'BEGIN{OFS="\t"}{ if($11 <= co){print $6,$7,$8,$9} }' >> "${GROUP}/${R2_PREFIX}.pool"
  done
done

#### Draw consistency plot for each comparison
DATAS=`echo ${P[@]} | sed -r 's/,//g'`
NPAIRS=`ls -1 ${GROUP}/*-npeaks-aboveIDR.txt | wc -l`
Rscript batch-consistency-plot.r ${NPAIRS} ${GROUP}/${OFPREFIX} ${DATAS}

#### Make final peak set for each sample
for (( I=0; I < ${#R[@]}; ++I ))
do
  PREFIX=`echo ${R[$I]} | awk -F "." '{print $1}'`
  POOLED=${PREFIX}.pool
  FINAL=${PREFIX}.final
  sort -k1,1 -k2,2n -k3,3n ${POOLED} | uniq | \
  bedtools intersect -a ${R[$I]} -b stdin -wa | \
  sort -k1,1 -k2,2n -k3,3n - | uniq > ${FINAL}
  rm ${POOLED}
done

#### Clean up all files
cd ${GROUP}
rm *-em.sav *-uri.sav *-npeaks-aboveIDR.txt *-Rout.txt *-overlapped-peaks.txt
