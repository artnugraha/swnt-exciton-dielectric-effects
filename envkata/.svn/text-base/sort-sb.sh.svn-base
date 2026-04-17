#!/bin/sh

MaxSubBand=6

rm -f family.out sort-sb.out pattern.out
ifort family.f90 -o family.out
ifort sort-sb.f90 -o sort-sb.out
ifort pattern.f90 -o pattern.out

# make family number 2n+m
for i in qex-kat-s*.dat
do
  cat $i | ./family.out | sort > ${i}2
done

rm -f fort.* e??s?.dat

# sort by "self-energy" - "binding energy"
type=0
while [ $type -lt 3 ]
do
  cat qex-kat-s${type}-*.dat2 | ./sort-sb.out
  for i in fort.*
    do
    cat $i | sort -k 10 > temp
    mv temp $i
    cline=1
    while [ $cline -le $MaxSubBand ]
      do
      cat $i 2> /dev/null | head -n $cline | tail -n 1 >> e${cline}${cline}s${type}.dat
      cline=`expr $cline + 1`
    done
  done
  rm -f fort.*

  echo "finish: s$type"
  type=`expr $type + 1`
done

# write data to files
type=0
while [ $type -lt 3 ]
  do
  cline=1
  while [ $cline -le $MaxSubBand ]
  do
    cat e${cline}${cline}s${type}.dat | sort | ./pattern.out > e${cline}${cline}s${type}.dat2
    cat e${cline}${cline}s${type}.dat2 | awk '{print $4,$6}' > e${cline}${cline}s${type}.dat2de
    cline=`expr $cline + 1`
  done
  type=`expr $type + 1`
done

rm -f family.out sort-sb.out pattern.out
