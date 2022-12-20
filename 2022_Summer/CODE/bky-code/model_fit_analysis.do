* Spillovers project.
* Evaluating model fit and size distribution of firms.
* By Jason Blevins, Ahmed Khwaja, and Nathan Yang.
* September 12, 2013.
************************************************************************************************

clear

set scheme spillovers

* Load the simulated data into STATA.
insheet using "model_fit_data.csv"

* Rename the variables.
rename v1  cityid
rename v2  time
rename v3  actualNAW
rename v4  actualNBK
rename v5  actualNHARV
rename v6  actualNMCD
rename v7  actualNWEND
rename v8  bblNAW
rename v9  bblNBK
rename v10 bblNHARV
rename v11 bblNMCD
rename v12 bblNWEND
rename v13 spilloversNAW
rename v14 spilloversNBK
rename v15 spilloversNHARV
rename v16 spilloversNMCD
rename v17 spilloversNWEND
rename v18 nospilloversNAW
rename v19 nospilloversNBK
rename v20 nospilloversNHARV
rename v21 nospilloversNMCD
rename v22 nospilloversNWEND
rename v23 bblNAWbias
rename v24 bblNBKbias
rename v25 bblNHARVbias
rename v26 bblNMCDbias
rename v27 bblNWENDbias
rename v28 spilloversNAWbias
rename v29 spilloversNBKbias
rename v30 spilloversNHARVbias
rename v31 spilloversNMCDbias
rename v32 spilloversNWENDbias
rename v33 nospilloversNAWbias
rename v34 nospilloversNBKbias
rename v35 nospilloversNHARVbias
rename v36 nospilloversNMCDbias
rename v37 nospilloversNWENDbias
rename v38 bblNAWmse
rename v39 bblNBKmse
rename v40 bblNHARVmse
rename v41 bblNMCDmse
rename v42 bblNWENDmse
rename v43 spilloversNAWmse
rename v44 spilloversNBKmse
rename v45 spilloversNHARVmse
rename v46 spilloversNMCDmse
rename v47 spilloversNWENDmse
rename v48 nospilloversNAWmse
rename v49 nospilloversNBKmse
rename v50 nospilloversNHARVmse
rename v51 nospilloversNMCDmse
rename v52 nospilloversNWENDmse

***********************************************************************
* Labels and Panel Structure
***********************************************************************

label variable time "Year"
label variable cityid "City ID"

gen city = "Abbotsford"
replace city = "Barrie" if cityid == 2
replace city = "Brampton" if cityid == 3
replace city = "Brantford" if cityid == 4
replace city = "Calgary" if cityid == 5
replace city = "Edmonton" if cityid == 6
replace city = "Guelph" if cityid == 7
replace city = "Halifax" if cityid == 8
replace city = "Hamilton" if cityid == 9
replace city = "Kelowna" if cityid == 10
replace city = "Kingston" if cityid == 11
replace city = "Kitchener" if cityid == 12
replace city = "London" if cityid == 13
replace city = "Moncton" if cityid == 14
replace city = "Montreal" if cityid == 15
replace city = "Niagra Falls" if cityid == 16
replace city = "Oshawa" if cityid == 17
replace city = "Ottawa" if cityid == 18
replace city = "Peterborough" if cityid == 19
replace city = "Quebec City" if cityid == 20
replace city = "Regina" if cityid == 21
replace city = "Saint John" if cityid == 22
replace city = "Saskatoon" if cityid == 23
replace city = "St. Johns" if cityid == 24
replace city = "Sudbury" if cityid == 25
replace city = "Thunder Bay" if cityid == 26
replace city = "Toronto" if cityid == 27
replace city = "Vancouver" if cityid == 28
replace city = "Victoria" if cityid == 29
replace city = "Windsor" if cityid == 30
replace city = "Winnipeg" if cityid == 31

sort cityid time
xtset cityid time

********************************************************************************
* Market share and HHI
********************************************************************************

* Actual data.
gen actualmarketshareAW = actualNAW/(actualNAW + actualNBK + actualNHARV + actualNMCD + actualNWEND)
gen actualmarketshareBK = actualNBK/(actualNAW + actualNBK + actualNHARV + actualNMCD + actualNWEND)
gen actualmarketshareHARV = actualNHARV/(actualNAW + actualNBK + actualNHARV + actualNMCD + actualNWEND)
gen actualmarketshareMCD = actualNMCD/(actualNAW + actualNBK + actualNHARV + actualNMCD + actualNWEND)
gen actualmarketshareWEND = actualNWEND/(actualNAW + actualNBK + actualNHARV + actualNMCD + actualNWEND)

* Plain BBL.
gen bblmarketshareAW = bblNAW/(bblNAW + bblNBK + bblNHARV + bblNMCD + bblNWEND)
gen bblmarketshareBK = bblNBK/(bblNAW + bblNBK + bblNHARV + bblNMCD + bblNWEND)
gen bblmarketshareHARV = bblNHARV/(bblNAW + bblNBK + bblNHARV + bblNMCD + bblNWEND)
gen bblmarketshareMCD = bblNMCD/(bblNAW + bblNBK + bblNHARV + bblNMCD + bblNWEND)
gen bblmarketshareWEND = bblNWEND/(bblNAW + bblNBK + bblNHARV + bblNMCD + bblNWEND)

* With Z.
gen spilloversmarketshareAW = spilloversNAW/(spilloversNAW + spilloversNBK + spilloversNHARV + spilloversNMCD + spilloversNWEND)
gen spilloversmarketshareBK = spilloversNBK/(spilloversNAW + spilloversNBK + spilloversNHARV + spilloversNMCD + spilloversNWEND)
gen spilloversmarketshareHARV = spilloversNHARV/(spilloversNAW + spilloversNBK + spilloversNHARV + spilloversNMCD + spilloversNWEND)
gen spilloversmarketshareMCD = spilloversNMCD/(spilloversNAW + spilloversNBK + spilloversNHARV + spilloversNMCD + spilloversNWEND)
gen spilloversmarketshareWEND = spilloversNWEND/(spilloversNAW + spilloversNBK + spilloversNHARV + spilloversNMCD + spilloversNWEND)

* With Z but no spillovers.
gen nospilloversmarketshareAW = nospilloversNAW/(nospilloversNAW + nospilloversNBK + nospilloversNHARV + nospilloversNMCD + nospilloversNWEND)
gen nospilloversmarketshareBK = nospilloversNBK/(nospilloversNAW + nospilloversNBK + nospilloversNHARV + nospilloversNMCD + nospilloversNWEND)
gen nospilloversmarketshareHARV = nospilloversNHARV/(nospilloversNAW + nospilloversNBK + nospilloversNHARV + nospilloversNMCD + nospilloversNWEND)
gen nospilloversmarketshareMCD = nospilloversNMCD/(nospilloversNAW + nospilloversNBK + nospilloversNHARV + nospilloversNMCD + nospilloversNWEND)
gen nospilloversmarketshareWEND = nospilloversNWEND/(nospilloversNAW + nospilloversNBK + nospilloversNHARV + nospilloversNMCD + nospilloversNWEND)

* Calculate the HHI index for each specification.
gen actualhhi = actualmarketshareAW^2 + actualmarketshareBK^2 + actualmarketshareHARV^2 + actualmarketshareMCD^2 + actualmarketshareWEND^2
gen bblhhi = bblmarketshareAW^2 + bblmarketshareBK^2 + bblmarketshareHARV^2 + bblmarketshareMCD^2 + bblmarketshareWEND^2
gen spillovershhi = spilloversmarketshareAW^2 + spilloversmarketshareBK^2 + spilloversmarketshareHARV^2 + spilloversmarketshareMCD^2 + spilloversmarketshareWEND^2
gen nospillovershhi = nospilloversmarketshareAW^2 + nospilloversmarketshareBK^2 + nospilloversmarketshareHARV^2 + nospilloversmarketshareMCD^2 + nospilloversmarketshareWEND^2

********************************************************************************
* Differences
********************************************************************************

gen spilloversNAW_diff = spilloversNAW - actualNAW
gen spilloversNBK_diff = spilloversNBK - actualNBK
gen spilloversNHARV_diff = spilloversNHARV - actualNHARV
gen spilloversNMCD_diff = spilloversNMCD - actualNMCD
gen spilloversNWEND_diff = spilloversNWEND - actualNWEND

********************************************************************************
* Model comparisons using store counts.
********************************************************************************

* Actual data.
bysort time: egen tactualNAW = mean(actualNAW)
bysort time: egen tactualNBK = mean(actualNBK)
bysort time: egen tactualNHARV = mean(actualNHARV)
bysort time: egen tactualNMCD = mean(actualNMCD)
bysort time: egen tactualNWEND = mean(actualNWEND)
bysort time: egen tactualhhi = mean(actualhhi)

* Plain BBL.
bysort time: egen tbblNAW = mean(bblNAW)
bysort time: egen tbblNBK = mean(bblNBK)
bysort time: egen tbblNHARV = mean(bblNHARV)
bysort time: egen tbblNMCD = mean(bblNMCD)
bysort time: egen tbblNWEND = mean(bblNWEND)
bysort time: egen tbblhhi = mean(bblhhi)

* With Z.
bysort time: egen tspilloversNAW = mean(spilloversNAW)
bysort time: egen tspilloversNBK = mean(spilloversNBK)
bysort time: egen tspilloversNHARV = mean(spilloversNHARV)
bysort time: egen tspilloversNMCD = mean(spilloversNMCD)
bysort time: egen tspilloversNWEND = mean(spilloversNWEND)
bysort time: egen tspillovershhi = mean(spillovershhi)

* With Z but no spillovers.
bysort time: egen tnospilloversNAW = mean(nospilloversNAW)
bysort time: egen tnospilloversNBK = mean(nospilloversNBK)
bysort time: egen tnospilloversNHARV = mean(nospilloversNHARV)
bysort time: egen tnospilloversNMCD = mean(nospilloversNMCD)
bysort time: egen tnospilloversNWEND = mean(nospilloversNWEND)
bysort time: egen tnospillovershhi = mean(nospillovershhi)

label variable tactualNAW "Actual"
label variable tactualNBK "Actual"
label variable tactualNHARV "Actual"
label variable tactualNMCD "Actual"
label variable tactualNWEND "Actual"
label variable tactualhhi "Actual"

label variable tbblNAW "No Z"
label variable tbblNBK "No Z"
label variable tbblNHARV "No Z"
label variable tbblNMCD "No Z"
label variable tbblNWEND "No Z"
label variable tbblhhi "No Z"

label variable tspilloversNAW "Z (with spillovers)"
label variable tspilloversNBK "Z (with spillovers)"
label variable tspilloversNHARV "Z (with spillovers)"
label variable tspilloversNMCD "Z (with spillovers)"
label variable tspilloversNWEND "Z (with spillovers)"
label variable tspillovershhi "Z (with spillovers)"

label variable tnospilloversNAW "Z (no spillovers)"
label variable tnospilloversNBK "Z (no spillovers)"
label variable tnospilloversNHARV "Z (no spillovers)"
label variable tnospilloversNMCD "Z (no spillovers)"
label variable tnospilloversNWEND "Z (no spillovers)"
label variable tnospillovershhi "Z (no spillovers)"

* Plot the simulated data against actual data.
twoway (line tactualNAW time, sort) (line tbblNAW time, sort) (line tspilloversNAW time, sort) (line tnospilloversNAW time, sort), ytitle("Outlets") yscale(range(0 30)) ylabel(0(5)30) saving(figures/modelfit_aw, replace)
twoway (line tactualNBK time, sort) (line tbblNBK time, sort) (line tspilloversNBK time, sort) (line tnospilloversNBK time, sort), ytitle("Outlets") yscale(range(0 30)) ylabel(0(5)30) saving(figures/modelfit_bk, replace)
twoway (line tactualNHARV time, sort) (line tbblNHARV time, sort) (line tspilloversNHARV time, sort) (line tnospilloversNHARV time, sort), ytitle("Outlets") yscale(range(0 30)) ylabel(0(5)30) saving(figures/modelfit_harv, replace)
twoway (line tactualNMCD time, sort) (line tbblNMCD time, sort) (line tspilloversNMCD time, sort) (line tnospilloversNMCD time, sort), ytitle("Outlets") yscale(range(0 30)) ylabel(0(5)30) saving(figures/modelfit_mcd, replace)
twoway (line tactualNWEND time, sort) (line tbblNWEND time, sort) (line tspilloversNWEND time, sort) (line tnospilloversNWEND time, sort), ytitle("Outlets") yscale(range(0 30)) ylabel(0(5)30) saving(figures/modelfit_wend, replace)
twoway (line tactualhhi time, sort) (line tbblhhi time, sort) (line tspillovershhi time, sort) (line tnospillovershhi time, sort), ytitle("HHI") saving(figures/modelfit_hhi, replace)

graph use figures/modelfit_aw.gph
graph export figures/modelfit_aw.eps, as(eps) replace

graph use figures/modelfit_bk.gph
graph export figures/modelfit_bk.eps, as(eps) replace

graph use figures/modelfit_harv.gph
graph export figures/modelfit_harv.eps, as(eps) replace

graph use figures/modelfit_mcd.gph
graph export figures/modelfit_mcd.eps, as(eps) replace

graph use figures/modelfit_wend.gph
graph export figures/modelfit_wend.eps, as(eps) replace

graph use figures/modelfit_hhi.gph
graph export figures/modelfit_hhi.eps, as(eps) replace

********************************************************************************
* Statistical Analysis of Model Fit
********************************************************************************

gen bblbias = bblNAWbias + bblNBKbias + bblNHARVbias + bblNMCDbias + bblNWENDbias
gen spilloversbias = spilloversNAWbias + spilloversNBKbias + spilloversNHARVbias + spilloversNMCDbias + spilloversNWENDbias
gen nospilloversbias = nospilloversNAWbias + nospilloversNBKbias + nospilloversNHARVbias + nospilloversNMCDbias + nospilloversNWENDbias

gen bblmse = bblNAWmse + bblNBKmse + bblNHARVmse + bblNMCDmse + bblNWENDmse
gen spilloversmse = spilloversNAWmse + spilloversNBKmse + spilloversNHARVmse + spilloversNMCDmse + spilloversNWENDmse
gen nospilloversmse = nospilloversNAWmse + nospilloversNBKmse + nospilloversNHARVmse + nospilloversNMCDmse + nospilloversNWENDmse

collapse (mean) bblbias spilloversbias nospilloversbias bblmse spilloversmse nospilloversmse, by(time)

label variable bblbias "No Z"
label variable bblmse "No Z"

label variable spilloversbias "Z (with spillovers)"
label variable spilloversmse "Z (with spillovers)"

label variable nospilloversbias "Z (no spillovers)"
label variable nospilloversmse "Z (no spillovers)"

twoway (tsline bblbias, lstyle(p2)) (tsline spilloversbias, lstyle(p3)) (tsline nospilloversbias, lstyle(p4)), ytitle("Mean Bias") saving(figures/modelfit_bias, replace)
graph use figures/modelfit_bias.gph
graph export figures/modelfit_bias.eps, as(eps) replace

twoway (tsline bblmse, lstyle(p2)) (tsline spilloversmse, lstyle(p3)) (tsline nospilloversmse, lstyle(p4)), ytitle("Mean Squared Error") saving(figures/modelfit_mse, replace)
graph use figures/modelfit_mse.gph
graph export figures/modelfit_mse.eps, as(eps) replace
