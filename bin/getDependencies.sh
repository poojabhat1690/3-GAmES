#!/bin/bash
module load singularity/3.2.1
#### pulling dependencies and slamdunk
echo "pulling dependencies for 3'GAmES"
singularity pull shub://poojabhat1690/dependencies

echo "pulling SLAMdunk"
singularity pull docker://tobneu/slamdunk:v0.3.4


