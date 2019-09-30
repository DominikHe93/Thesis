#!/bin/bash

echo " bamm parse -b ./mapping_*/*.bam -c output_coverage.tab -m opmean --max_distance 50 -t 4"

 bamm parse -b ./mapping_*/*.bam -c coverage.tab -m opmean --max_distance 50 -t 4

echo "done"


# -b bamfile input -c name to write in tab -m opmean avaerage coverage --max_distance  maximum allowable edit distance from query to reference -t used threats
