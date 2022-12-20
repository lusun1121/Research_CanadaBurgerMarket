#!/bin/sh

# No Z
matlab -nosplash -nodesktop -r spec1_first
matlab -nosplash -nodesktop -r spec1_simulate
matlab -nosplash -nodesktop -r spec1_second
matlab -nosplash -nodesktop -r spec1_bootstrap

# Z (No Spillovers)
matlab -nosplash -nodesktop -r spec2_first
matlab -nosplash -nodesktop -r spec2_simulate
matlab -nosplash -nodesktop -r spec2_second
matlab -nosplash -nodesktop -r spec2_bootstrap

# Z (Spillovers)
matlab -nosplash -nodesktop -r spec3_first
matlab -nosplash -nodesktop -r spec3_simulate
matlab -nosplash -nodesktop -r spec3_second
matlab -nosplash -nodesktop -r spec3_bootstrap
