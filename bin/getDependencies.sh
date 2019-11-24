#!/bin/bash
#module load singularity/3.2.1
#### pulling dependencies and slamdunk
command -v singularity >/dev/null 2>&1 || { echo >&2 "3' GAmES requires singularity version > 3.0, please load this and try again."; exit 1; }


echo "pulling dependencies for 3'GAmES"
singularity pull shub://poojabhat1690/dependencies

echo "pulling SLAMdunk"
singularity pull docker://tobneu/slamdunk:v0.3.4


