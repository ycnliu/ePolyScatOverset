#!/usr/bin/python

import sys

Timings={}
for line in sys.stdin:
   if line[:10] == "Time Now =" :
      TimeNowData = line.split()
      DeltaTime = float(TimeNowData[7])
      Label = " ".join(TimeNowData[8:])
      if Label in Timings :
         Timings[Label]+=DeltaTime
      else :
         Timings[Label]=DeltaTime

TimingsSorted=sorted([(key,val) for key,val in Timings.items()], key=lambda item: item[1], reverse=True)
for item in TimingsSorted :
   print "{0:>10.2f} {1:<40s}".format(item[1], item[0])

