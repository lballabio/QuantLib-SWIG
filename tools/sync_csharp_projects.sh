#!/bin/bash

# execute this script from the root of an uncompressed QuantLib-SWIG tarball

# get reference lists of distributed files (done with find; this is
# why this script should be run from an uncompressed tarball created
# with 'make dist', not from a working copy.)

find CSharp/csharp -name '*.cs' \
| awk -F'/' '{ print $3 }' \
| sort > csharp.ref.files

# Extract file names from VC# project.

grep -o -E 'Compile *Include *= *".*"' CSharp/csharp/NQuantLib.csproj \
| awk -F'"' '{ print $2 }' \
| sort > csharp.vc1x.files

# Write out differences...

echo 'Visual Studio:'
diff -b csharp.vc1x.files csharp.ref.files

# ...and cleanup
rm -f csharp.*.files

