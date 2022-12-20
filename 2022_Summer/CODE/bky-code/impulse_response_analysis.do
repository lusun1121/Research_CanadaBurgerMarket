* Spillovers project.
* Impulse response analysis for Z's and store counts.
* By Jason Blevins, Ahmed Khwaja, and Nathan Yang.
* January 26, 2014.
************************************************************************************************

clear
insheet using "impulse_response_I_data.csv"
do impulse_response_common
tsline diffNMCD diffNBK diffNAW diffNHARV diffNWEND if tin(2000, 2005), ytitle("Difference in Net New Outlets") saving(figures/iri, replace)
graph export figures/iri.eps, as(eps) replace

clear
insheet using "impulse_response_D_data.csv"
do impulse_response_common
tsline diffNMCD diffNBK diffNAW diffNHARV diffNWEND if tin(2000, 2005), ytitle("Difference in Net New Outlets") saving(figures/ird, replace)
graph export figures/ird.eps, as(eps) replace

clear
insheet using "impulse_response_Z_data.csv"
do impulse_response_common
tsline diffNMCD diffNBK diffNAW diffNHARV diffNWEND if tin(2000, 2005), ytitle("Difference in Net New Outlets") saving(figures/irz, replace)
graph export figures/irz.eps, as(eps) replace
list diff* if time == 2005

clear
insheet using "impulse_response_ZN_data.csv"
do impulse_response_common
tsline diffNMCD diffNBK diffNAW diffNHARV diffNWEND if tin(2000, 2005), ytitle("Difference in Net New Outlets") saving(figures/irzn, replace)
graph export figures/irzn.eps, as(eps) replace
list diff* if time == 2005
