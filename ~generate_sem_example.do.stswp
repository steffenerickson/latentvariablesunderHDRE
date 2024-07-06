//----------------------------------------------------------------------------//
// Title  : HDRE SEM simulation example 
// Purpose: Simulates two variable model under HDRE
// Author : Steffen Erickson
// Date   : 7/3/24
//----------------------------------------------------------------------------//

clear all 
frame reset 
include generate_sem.do

//----------------------------------------------------------------------------//
// Generate Data
//----------------------------------------------------------------------------//
global factor1 s1 s2          // s skills 
global factor2 r1 r2 r3 r4    // r rater 
global factor3 l1 l2 l3 l4    // l lesson                
foreach i of global factor1 {
foreach j of global factor2 {
foreach k of global factor3 {
local temp `i'_`j'_`k'
local all_combos : list all_combos | temp 
}
}
}
global all_combos `all_combos'
loadingpattern $all_combos
global factors = 	r(vars)
mata Lambda_x = edit_param_mat(st_matrix("r(loadingmat)"),1,.5,.99)
mata Theta_x = J(rows(Lambda_x),rows(Lambda_x),0)
mata diag  = J(rows(Lambda_x),1,1) 
mata _diag(Theta_x,diag)
mata Theta_x = edit_param_mat(Theta_x ,0,.5,.99) 
mata length(tokens(st_global("factors")))
mata Phi = J(length(tokens(st_global("factors"))),length(tokens(st_global("factors"))),0)
mata diag  = J(rows(length(tokens(st_global("factors")))),1,1) 
mata _diag(Phi,diag)
mata Phi = edit_param_mat(Phi,0,1,2) 
mata Phi[2,1] = Phi[1,2] =.97 // hard coded covariance cov(s1,s2)
mata kappa =  J(rows(Phi),1,1)
mata kappa = edit_param_mat(kappa,0,1,2) 
mata s1 = generate_sem_cov(Lambda_x,Theta_x,Phi)
mata Sigma = s1.Sigma
mata s2 = generate_sem_mean(Lambda_x,kappa)
mata upsilon_y = s2.upsilon_y
mata upsilon_x = s2.upsilon_x
mata mu = upsilon_y \ upsilon_x
mata st_matrix("Sigma",Sigma)
mata st_matrix("mu",mu)
drawnorm ${all_combos}, cov(Sigma) means(mu) n(1000)

//----------------------------------------------------------------------------//
// Models 
//----------------------------------------------------------------------------//
// Correlated error model
qui ds
local list `r(varlist)'
matrix  Theta_x = J(`: word count `list'',`: word count `list'',0)
forvalues i = 1/`: word count `list'' { // hard coded search
forvalues j = 1/`: word count `list'' {
local first1  = substr("`:word `i' of `list''",4,2)
local second1 = substr("`:word `j' of `list''",4,2)
local first2  = substr("`:word `i' of `list''",7,2)
local second2 = substr("`:word `j' of `list''",7,2)
if "`first1'" == "`second1'"  |  "`first2'" == "`second2'"  mat Theta_x[`i',`j'] = .
}
}
sem (F1 -> s1*) (F2 -> s2*) ,covstructure(e._OEn, fixed(Theta_x))
di _b[/cov(F1,F2)] // estimate
mata Phi[2,1] 	   // true 

// Ommitting the correlated errors 
sem (F1 -> s1*) (F2 -> s2*) 
di _b[/cov(F1,F2)]   // estimate
mata Phi[2,1] 		 // true 

// Observed mean score covariance 
egen s1 = rowmean(s1*)
egen s2 = rowmean(s2*)
corr s1 s2 , cov // estimate
mata Phi[2,1] 	 // true 

//----------------------------------------------------------------------------//
//Simulation  
//----------------------------------------------------------------------------//
forvalues i = 1/20 {
clear 
drawnorm ${all_combos}, cov(Sigma) means(mu) n(1000)
sem (F1 -> s1*) (F2 -> s2*) ,covstructure(e._OEn, fixed(Theta_x))
mat results = (nullmat(results),_b[/cov(F1,F2)])
}
mata  mean(st_matrix("results")')
mata Phi[2,1] 	// true 





















































