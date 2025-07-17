#!/bin/bash
#
#
for FILE ; do
echo ----------------------------------------------
echo For file $FILE
$pe/bin/FindRot.exe <$FILE
done
#