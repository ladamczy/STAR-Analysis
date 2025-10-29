#!/bin/bash

# Navigate to the star-upc-new directory relative to current location
cd ../star-upc-new

# Check if the directory was found and build script exists
if [ $? -ne 0 ]; then
  echo "Error: Could not navigate to ../star-upc-new"
  exit 1
fi

if [ ! -f "./buildLib.sh" ]; then
  echo "Error: buildLib.sh not found in ../star-upc-new"
  exit 1
fi

# Run the build script
echo "Running buildLib.sh..."
./buildLib.sh

# Optional: return to original directory (uncomment if needed)
# cd - > /dev/null
### End of script
