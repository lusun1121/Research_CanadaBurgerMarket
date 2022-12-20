#!/bin/bash

mkdir -p figures
matlab -nosplash -nodesktop -r counterfactual_dynamics >& counterfactual_dynamics.log
stata -b do model_fit_analysis
stata -b do impulse_response_analysis
pushd figures
for i in *.eps; do epstopdf $i; done
popd
