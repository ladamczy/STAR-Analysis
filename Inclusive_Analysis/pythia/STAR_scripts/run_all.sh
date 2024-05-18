#!/bin/bash

starver pro

root4star -q -b elastic.C\($1\)
root4star -q -b singleDiffractive.C\($1\)
root4star -q -b doubleDiffractive.C\($1\)
root4star -q -b centralDiffractive.C\($1\)
root4star -q -b nonDiffractive.C\($1\)