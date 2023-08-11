 #!/bin/bash

for FILE in $(ls | grep star-upc-f)
do
    ext=${FILE##*star-upc-f} # delete everything up to star-upc-f
    upcase=F$ext # uppercase only F
    cp $FILE $upcase
done
