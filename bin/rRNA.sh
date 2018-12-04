#!/bin/bash -x
set -eu -o pipefail 
genome=$(basename $gdir)
genome=${genome%%_*}
##remove the version

genome=${genome%%[0-9]*}
rRNA=/sc/orga/projects/sealfs01a/motrpac/refdata/rRNA/${genome}_rRNA


