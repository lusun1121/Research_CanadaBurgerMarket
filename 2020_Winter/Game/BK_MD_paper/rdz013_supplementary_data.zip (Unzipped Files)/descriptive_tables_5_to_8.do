/* ------------------------------------------------------ */
/* STATA do file that obtains the descriptive statistics  */
/* in Tables 5, 6, 7, and 8 in the paper				  */
/* "Identi?cation and Estimation of Dynamic Games         */
/*  when Players?Beliefs Are Not in Equilibrium"          */
/*  by Victor Aguirregabiria and Arvind Magesan           */
/* at the Review of Economic Studies, 2019                */
/* ------------------------------------------------------ */

capture log close
set matsize 500

log using D:\RACHON\Dropbox\MYPAPERS\ARVIND_RATIONALIZABILITY\prostata\descriptive_tables_5_to_8.log, replace

* -----------------
* 1. Reading data 
* -----------------
use "D:\RACHON\Dropbox\MYPAPERS\ARVIND_RATIONALIZABILITY\data\toivanen_waterson_nolondon_20022019.dta", clear

* -------------------------------------------------------
* 2. Table 5
*    Descriptive Statistics on Local Markets (Year 1991)
* -------------------------------------------------------
gen prop_0514 = pop_0514/population
gen prop_1529 = pop_1529/population
gen prop_4559 = pop_4559/population
gen prop_6064 = pop_6064/population
gen prop_6574 = pop_6574/population
gen ubrate = ue/(population * 1000)

* TABLE 5: Area (square km)
* -------------------------
sum district_area if year==1991, detail

* TABLE 5: Population (thousands)
* -------------------------------
sum population if year==1991, detail

* TABLE 5: Share of children: Age 5-14
* ------------------------------------
sum prop_0514 if year==1991, detail

* TABLE 5: Share of Young: 15-29
* ------------------------------
sum prop_1529 if year==1991, detail

* TABLE 5: Share of Pensioners: 65-74
* -----------------------------------
sum prop_6574 if year==1991, detail

* TABLE 5: GDP per capita (thousand £)
* ------------------------------------
sum gdp_pc if year==1991, detail

* TABLE 5: Claimants of UB / Population ratio
* -------------------------------------------
sum ubrate if year==1991, detail

* TABLE 5: Avg. Weekly Rent per dwelling (£)
* ------------------------------------------
sum avg_rent if year==1991, detail

* TABLE 5: Council tax (£)
* ------------------------
sum ctax if year==1991, detail

* TABLE 5: Number of BK stores [Year 1990]
* ----------------------------------------
sum bk_stock if year==1991, detail

* TABLE 5: Number of MD stores [Year 1990]
* ----------------------------------------
sum mcd_stock if year==1991, detail


* -------------------------------------------------------
* 3. Table 6
*    Evolution of the Number of Stores
* -------------------------------------------------------
egen bk_markets_with_stores = sum((bk_stock>0)), by(year)
egen mcd_markets_with_stores = sum((mcd_stock>0)), by(year)

egen bk_new_markets_with_stores = sum((bk_entry>0)*(bk_stock==0)), by(year)
egen mcd_new_markets_with_stores = sum((mcd_entry>0)*(mcd_stock==0)), by(year)

egen bk_total_stores_stock = sum(bk_stock), by(year)
egen mcd_total_stores_stock = sum(mcd_stock), by(year)

egen bk_total_new_stores = sum(bk_entry), by(year)
egen mcd_total_new_stores = sum(mcd_entry), by(year)

gen bk_stores_per_market_cond = bk_total_stores_stock/bk_markets_with_stores
gen mcd_stores_per_market_cond = mcd_total_stores_stock/mcd_markets_with_stores

* TABLE 6: BK #Markets with Stores
* --------------------------------
tab year, sum(bk_markets_with_stores)

* TABLE 6: BK Change in #Markets with Stores
* ------------------------------------------
tab year, sum(bk_new_markets_with_stores)

* TABLE 6: BK # of Stores
* -----------------------
tab year, sum(bk_total_stores_stock)

* TABLE 6: BK Change in # of Stores
* ---------------------------------
tab year, sum(bk_total_new_stores)

* TABLE 6: BK Mean #Stores per Market(Conditional on #Stores>0)
* -------------------------------------------------------------
tab year, sum(bk_stores_per_market_cond)

* TABLE 6: MCD #Markets with Stores
* --------------------------------
tab year, sum(mcd_markets_with_stores)

* TABLE 6: MCD Change in #Markets with Stores
* ------------------------------------------
tab year, sum(mcd_new_markets_with_stores)

* TABLE 6: MCD # of Stores
* -----------------------
tab year, sum(mcd_total_stores_stock)

* TABLE 6: MCD Change in # of Stores
* ---------------------------------
tab year, sum(mcd_total_new_stores)

* TABLE 6: MCD Mean #Stores per Market(Conditional on #Stores>0)
* -------------------------------------------------------------
tab year, sum(mcd_stores_per_market_cond)


* -------------------------------------------------------
* 4. Table 7
*    Transition Probability Matrix for Market Structure
* -------------------------------------------------------
gen bk_current_stores = bk_stock + bk_entry
gen mcd_current_stores = mcd_stock + mcd_entry

* Generating state at period t
gen state_0 = .
replace state_0 = 1 if (bk_stock==0) & (mcd_stock==0) 
replace state_0 = 2 if (bk_stock==0) & (mcd_stock==1) 
replace state_0 = 3 if (bk_stock==0) & (mcd_stock>=2) 
replace state_0 = 4 if (bk_stock==1) & (mcd_stock==0)
replace state_0 = 5 if (bk_stock==1) & (mcd_stock==1) 
replace state_0 = 6 if (bk_stock==1) & (mcd_stock>=2) 
replace state_0 = 7 if (bk_stock>=2) & (mcd_stock==0) 
replace state_0 = 8 if (bk_stock>=2) & (mcd_stock==1) 
replace state_0 = 9 if (bk_stock>=2) & (mcd_stock>=2) 

* Generating state at period t+1
gen state_1 = .
replace state_1 = 1 if (bk_current_stores==0) & (mcd_current_stores==0) 
replace state_1 = 2 if (bk_current_stores==0) & (mcd_current_stores==1) 
replace state_1 = 3 if (bk_current_stores==0) & (mcd_current_stores>=2) 
replace state_1 = 4 if (bk_current_stores==1) & (mcd_current_stores==0)
replace state_1 = 5 if (bk_current_stores==1) & (mcd_current_stores==1) 
replace state_1 = 6 if (bk_current_stores==1) & (mcd_current_stores>=2) 
replace state_1 = 7 if (bk_current_stores>=2) & (mcd_current_stores==0) 
replace state_1 = 8 if (bk_current_stores>=2) & (mcd_current_stores==1) 
replace state_1 = 9 if (bk_current_stores>=2) & (mcd_current_stores>=2) 

* TABLE 7: Matrix of Transition Probabilities
* -------------------------------------------
tab state_0 state_1, row

* ---------------------------------------------------------
* 5. Table 5
*    Reduced Form Probits for the Decision to Open a Store
* ---------------------------------------------------------
gen lpop = ln(population)
gen lpdens = ln(population/district_area)
gen md_pos_stock = (mcd_stock>0) 
gen bk_pos_stock = (bk_stock>0)
gen lgdp = ln(gdp_pc)
gen lrent = ln(avg_rent)
gen ldistbk  = ln(dist_bkhq_minu)
gen ldistmd  = ln(dist_mdhq_minu)

gen md_stock_1 = (mcd_stock==1)
gen md_stock_2 = (mcd_stock==2)
gen md_stock_3 = (mcd_stock>=3)
gen bk_stock_1 = (bk_stock==1)
gen bk_stock_2 = (bk_stock==2)
gen bk_stock_3 = (bk_stock>=3)

tab year, gen(timed)
tab county_code, gen(cound)
tab district_code, gen(dumdis) 

* -------------------------------------------------------- 
* TABLE 8: Probit Estimates  without county fixed effects
*          BURGER KING
* -------------------------------------------------------- 
dprobit bk_entdum bk_stock_1 bk_stock_2 bk_stock_3 md_stock_1 md_stock_2 md_stock_3 lpop lpdens prop_0514 prop_1529 prop_6574 lgdp ubrate ctax  lrent ldistmd ldistbk timed*

* -------------------------------------------------------- 
* TABLE 8: Probit Estimates  WITH county fixed effects
*          BURGER KING
* -------------------------------------------------------- 
dprobit bk_entdum bk_stock_1 bk_stock_2 bk_stock_3 md_stock_1 md_stock_2 md_stock_3 lpop lpdens prop_0514 prop_1529 prop_6574 lgdp ubrate ctax lrent ldistmd ldistbk timed* cound*

* -------------------------------------------------------- 
* TABLE 8: Probit Estimates  WITH district fixed effects
*          BURGER KING
* -------------------------------------------------------- 
dprobit bk_entdum bk_stock_1 bk_stock_2 bk_stock_3 md_stock_1 md_stock_2 md_stock_3 lpop lpdens prop_0514 prop_1529 prop_6574 lgdp ubrate ctax lrent ldistmd ldistbk timed* dumdis*

* -------------------------------------------------------- 
* TABLE 8: Probit Estimates  without county fixed effects
*          MCDONALDS 
* -------------------------------------------------------- 
dprobit mcd_entdum md_stock_1 md_stock_2 md_stock_3 bk_stock_1 bk_stock_2 bk_stock_3 lpop lpdens prop_0514 prop_1529 prop_6574 lgdp ubrate ctax lrent ldistmd ldistbk timed*

* -------------------------------------------------------- 
* TABLE 8: Probit Estimates  WITH county fixed effects
*          MCDONALDS 
* -------------------------------------------------------- 
dprobit mcd_entdum md_stock_1 md_stock_2 md_stock_3 bk_stock_1 bk_stock_2 bk_stock_3 lpop lpdens prop_0514 prop_1529 prop_6574 lgdp ubrate ctax lrent ldistmd ldistbk  timed* cound*

* -------------------------------------------------------- 
* TABLE 8: Probit Estimates  WITH district fixed effects
*          MCDONALDS 
* -------------------------------------------------------- 
dprobit mcd_entdum md_stock_1 md_stock_2 md_stock_3 bk_stock_1 bk_stock_2 bk_stock_3 lpop lpdens prop_0514 prop_1529 prop_6574 lgdp ubrate ctax lrent ldistmd ldistbk timed* dumdis*


log close
