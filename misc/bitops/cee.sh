#! /bin/bash

i=1
while [ $i -lt 1000 ];
do
   let K=$i
	# echo $K
   ./psileft $K 90
	let i=$i+1
done

