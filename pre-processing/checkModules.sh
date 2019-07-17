#!/bin/sh

MODULE="$1"

echo "$MODULE"

cat module list | grep "$MODULE"


#if module list | grep "$MODULE" &> /dev/null ; then
#	echo "$MODULE is loaded!"
#	exit 0
 #   else
#	 echo "$MODULE is not loaded!"
#	 exit 1
#	fi




















