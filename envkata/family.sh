#!/bin/sh

ifort family.f90 -o family.out
ifort pattern.f90 -o pattern.out

for i in qex-kat-*.dat
do
  cat ${i} | ./family.out | sort | ./pattern.out > ${i}2
  cat ${i}2 | awk '{print $4,$6,$9,$7,$10}' > ${i}3
  cat ${i}2 | awk '{print $4,$6}' > ${i}3de
  cat ${i}2 | awk '{print $4,$9}' > ${i}3ds
  cat ${i}2 | awk '{print $4,$7}' > ${i}3db
  cat ${i}2 | awk '{print $4,$10}' > ${i}3dc
done

rm -f family.out pattern.out
