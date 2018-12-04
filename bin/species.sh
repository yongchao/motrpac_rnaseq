#!/bin/bash
set -eu -o pipefail
gdir=$1
genome=$(basename $gdir)
genome=${genome%%_*}
##remove the version
echo ${genome%%[0-9]*}
