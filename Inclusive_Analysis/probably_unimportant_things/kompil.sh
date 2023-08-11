#!/usr/bin/bash

mkdir build
cp test.root ../build
cd build
cmake ../
make
./Ana ../test.list ./Ana_Output.root

#czekanie na klawisz
echo "Press any key to continue"
while [ true ] ; do
read -t 3 -n 1
if [ $? = 0 ] ; then
exit ;
#else
#echo "waiting for the keypress"
fi
done