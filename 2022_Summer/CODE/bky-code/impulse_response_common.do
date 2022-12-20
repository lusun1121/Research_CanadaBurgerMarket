set scheme spillovers

* Rename the variables.
rename v1  cityid
rename v2  time
rename v3  actualNAW
rename v4  actualNBK
rename v5  actualNHARV
rename v6  actualNMCD
rename v7  actualNWEND
rename v8  spilloversNAW
rename v9  spilloversNBK
rename v10 spilloversNHARV
rename v11 spilloversNMCD
rename v12 spilloversNWEND
rename v13 irAW_NAW
rename v14 irAW_NBK
rename v15 irAW_NHARV
rename v16 irAW_NMCD
rename v17 irAW_NWEND
rename v18 irBK_NAW
rename v19 irBK_NBK
rename v20 irBK_NHARV
rename v21 irBK_NMCD
rename v22 irBK_NWEND
rename v23 irHARV_NAW
rename v24 irHARV_NBK
rename v25 irHARV_NHARV
rename v26 irHARV_NMCD
rename v27 irHARV_NWEND
rename v28 irMCD_NAW
rename v29 irMCD_NBK
rename v30 irMCD_NHARV
rename v31 irMCD_NMCD
rename v32 irMCD_NWEND
rename v33 irWEND_NAW
rename v34 irWEND_NBK
rename v35 irWEND_NHARV
rename v36 irWEND_NMCD
rename v37 irWEND_NWEND
rename v38 spilloversPIAW
rename v39 spilloversPIBK
rename v40 spilloversPIHARV
rename v41 spilloversPIMCD
rename v42 spilloversPIWEND
rename v43 irAW_PIAW
rename v44 irAW_PIBK
rename v45 irAW_PIHARV
rename v46 irAW_PIMCD
rename v47 irAW_PIWEND
rename v48 irBK_PIAW
rename v49 irBK_PIBK
rename v50 irBK_PIHARV
rename v51 irBK_PIMCD
rename v52 irBK_PIWEND
rename v53 irHARV_PIAW
rename v54 irHARV_PIBK
rename v55 irHARV_PIHARV
rename v56 irHARV_PIMCD
rename v57 irHARV_PIWEND
rename v58 irMCD_PIAW
rename v59 irMCD_PIBK
rename v60 irMCD_PIHARV
rename v61 irMCD_PIMCD
rename v62 irMCD_PIWEND
rename v63 irWEND_PIAW
rename v64 irWEND_PIBK
rename v65 irWEND_PIHARV
rename v66 irWEND_PIMCD
rename v67 irWEND_PIWEND

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
* Impulse response analysis of store counts
********************************************************************************

* Actual data.
bysort time: egen tactualNAW = total(actualNAW)
bysort time: egen tactualNBK = total(actualNBK)
bysort time: egen tactualNHARV = total(actualNHARV)
bysort time: egen tactualNMCD = total(actualNMCD)
bysort time: egen tactualNWEND = total(actualNWEND)

* Simulated store counts
bysort time: egen tspilloversNAW = total(spilloversNAW)
bysort time: egen tspilloversNBK = total(spilloversNBK)
bysort time: egen tspilloversNHARV = total(spilloversNHARV)
bysort time: egen tspilloversNMCD = total(spilloversNMCD)
bysort time: egen tspilloversNWEND = total(spilloversNWEND)

* Simulated store counts with shock
bysort time: egen tirAW_NAW = total(irAW_NAW)
bysort time: egen tirBK_NBK = total(irBK_NBK)
bysort time: egen tirHARV_NHARV = total(irHARV_NHARV)
bysort time: egen tirMCD_NMCD = total(irMCD_NMCD)
bysort time: egen tirWEND_NWEND = total(irWEND_NWEND)

* Simulated profits
bysort time: egen tspilloversPIAW = total(spilloversPIAW)
bysort time: egen tspilloversPIBK = total(spilloversPIBK)
bysort time: egen tspilloversPIHARV = total(spilloversPIHARV)
bysort time: egen tspilloversPIMCD = total(spilloversPIMCD)
bysort time: egen tspilloversPIWEND = total(spilloversPIWEND)

* Simulated profits with shock
bysort time: egen tirAW_PIAW = total(irAW_PIAW)
bysort time: egen tirBK_PIBK = total(irBK_PIBK)
bysort time: egen tirHARV_PIHARV = total(irHARV_PIHARV)
bysort time: egen tirMCD_PIMCD = total(irMCD_PIMCD)
bysort time: egen tirWEND_PIWEND = total(irWEND_PIWEND)

* Collapse to a time series
drop if time < 1999
drop if cityid > 1
keep time tspillovers* tir*
sort time
tsset time

* Net new stores built at baseline
gen dtspillovers_NAW = d.tspilloversNAW
gen dtspillovers_NBK = d.tspilloversNBK
gen dtspillovers_NHARV = d.tspilloversNHARV
gen dtspillovers_NMCD = d.tspilloversNMCD
gen dtspillovers_NWEND = d.tspilloversNWEND

* Net new stores built under shock
gen dtir_NAW = d.tirAW_NAW
gen dtir_NBK = d.tirBK_NBK
gen dtir_NHARV = d.tirHARV_NHARV
gen dtir_NMCD = d.tirMCD_NMCD
gen dtir_NWEND = d.tirWEND_NWEND

* Differences
gen diffNAW = dtir_NAW - dtspillovers_NAW
gen diffNBK = dtir_NBK - dtspillovers_NBK
gen diffNHARV = dtir_NHARV - dtspillovers_NHARV
gen diffNMCD = dtir_NMCD - dtspillovers_NMCD
gen diffNWEND = dtir_NWEND - dtspillovers_NWEND

label variable diffNAW "A&W"
label variable diffNBK "Burger King"
label variable diffNHARV "Harvey's"
label variable diffNMCD "McDonald's"
label variable diffNWEND "Wendy's"
