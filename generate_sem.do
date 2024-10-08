//----------------------------------------------------------------------------//
// Title  : SEM simulation programs 
// Purpose: Programs for simulating data using structural equations with latent variables 
// Author : Steffen Erickson
// Date   : 6/10/24
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Generates Structural Equation Models using LISERAL notation and 
// Equations from Bollen, 1989
/*
	/* Inputs */ 
	Endogenous 
	matrix Lambda_y...................(p x m) measurement loading matrix 
	matrix Theta_y....................(p x p) measurement error matrix
	matrix Beta.......................(m x m) latent structural coefficients
	matrix Psi........................(m x m) latent variance/covariance
	vector alpha .....................(m x 1) latent mean vector 
	
	Exogenous
	matrix Lambda_x...................(q x n) measurement loading matrix 
	matrix Theta_x....................(q x q) measurement error matrix 
	matrix Gamma......................(m x n) latent structural coefficients
	matrix Phi........................(n x n) latent variance/covariance
	vector kappa......................(n x 1) latent mean vector 
	
	/* Outputs */ 
	Matrix Sigma......................(m x n) observed covariance matrix 
	vector upsilon_y..................(p x 1) observed endogenous mean vector 
	vector upsilon_x..................(q x 1) observed exogenous mean vector 


*/
//----------------------------------------------------------------------------//

mata 
// ---------------- Define Structure -----------------------------------------//
struct myproblem  {
	real matrix 	Lambda_y 
	real matrix 	Theta_y  
	real matrix 	Beta     
	real matrix 	Psi      
	real matrix 	Lambda_x 
	real matrix 	Theta_x  
	real matrix 	Gamma    
    real matrix 	Phi  
	real colvector 	kappa
	real colvector 	alpha 
	struct derived 	scalar d
}
struct derived { 
	real matrix 	Sigma 
	real matrix 	YY 
	real matrix 	XX 
	real matrix 	YX
	real matrix 	XY 
	real colvector  upsilon_y
	real colvector  upsilon_x
	real matrix     IminusBeta 
}
// ---------------- Population Covariance Matrix -----------------------------//			
struct derived generate_sem_cov(real matrix Lambda_x,
								real matrix Theta_x, 
								real matrix Phi, 
								| real matrix Lambda_y,
								real matrix Theta_y,
								real matrix Beta,
								real matrix Gamma, 
								real matrix Psi) 
{
	struct myproblem scalar pr
	initialize_objects_cov(pr.d)
	
	pr.Lambda_x = Lambda_x 
	pr.Theta_x  = Theta_x  
	pr.Phi      = Phi 
	pr.Lambda_y = Lambda_y 
	pr.Theta_y  = Theta_y  
	pr.Beta     = Beta      
	pr.Gamma    = Gamma
	pr.Psi      = Psi     
	pr.d.IminusBeta = I(rows(pr.Beta),rows(pr.Beta))- pr.Beta
	
	if (args() <= 3) {
		get_xx(pr)
		pr.d.Sigma = pr.d.XX
		_makesymmetric(pr.d.Sigma)
	}
	else {
		get_yy(pr)
		get_xx(pr)
		get_yx(pr)
		get_xy(pr)
		pr.d.Sigma = (pr.d.YY,pr.d.YX\pr.d.XY,pr.d.XX)
		_makesymmetric(pr.d.Sigma)
	}
	return(pr.d)
}

void initialize_objects_cov(struct derived scalar d)
{
	d.Sigma = J(0,0,.) 
	d.YY = J(0,0,.) 
	d.XX = J(0,0,.)
	d.YX = J(0,0,.) 
	d.XY = J(0,0,.) 
}

void get_xx(struct myproblem scalar pr)
{
  pr.d.XX = pr.Lambda_x*pr.Phi*pr.Lambda_x' + pr.Theta_x
}
void get_yy(struct myproblem scalar pr)
{
	pr.d.YY = pr.Lambda_y*luinv(pr.d.IminusBeta)*(pr.Gamma*pr.Phi*pr.Gamma' + pr.Psi)*luinv(pr.d.IminusBeta)'*pr.Lambda_y' + pr.Theta_y
}
void get_yx(struct myproblem scalar pr)
{
	pr.d.YX = pr.Lambda_y*luinv(pr.d.IminusBeta)*pr.Gamma*pr.Phi*pr.Lambda_x'
}
void get_xy(struct myproblem scalar pr) // the transpose of yx
{
	pr.d.XY = pr.Lambda_x*pr.Phi*pr.Gamma'*luinv(pr.d.IminusBeta)'pr.Lambda_y'
}

// ---------------- Population Mean Vector -----------------------------------//
struct derived generate_sem_mean(real matrix Lambda_x,
								 real vector kappa,
								| real matrix Lambda_y,
								real matrix Beta,
								real matrix Gamma, 
								real vector alpha) 
{
	struct myproblem scalar pr
	initialize_objects_mean(pr.d)
	row_to_col(kappa)
	row_to_col(alpha)
	
	pr.Lambda_x 	= Lambda_x 
	pr.kappa 		= kappa 
	pr.Lambda_y 	= Lambda_y 
	pr.Beta 		= Beta     
	pr.Gamma 		= Gamma   
	pr.alpha 		= alpha 
	pr.d.IminusBeta = I(rows(pr.Beta),rows(pr.Beta))- pr.Beta
	
	if (args() <= 2) {
		get_upsilon_x(pr)
	}
	else {
		get_upsilon_x(pr)
		get_upsilon_y(pr)
	}
	return(pr.d)
}

void initialize_objects_mean(struct derived scalar d)
{
	d.upsilon_y = J(0,1,.)
	d.upsilon_x = J(0,1,.)
}
void get_upsilon_x(struct myproblem scalar pr)
{
	pr.d.upsilon_x = pr.Lambda_x*pr.kappa
}
void get_upsilon_y(struct myproblem scalar pr)
{
	pr.d.upsilon_y = pr.Lambda_y*luinv(pr.d.IminusBeta)*(pr.alpha + pr.Gamma*pr.kappa)
}

real scalar is_colvec(z) return(anyof(("colvector","scalar"),orgtype(z)))
void row_to_col(real vector v) 
{
	if (is_colvec(v) == 0) v = v'
}
end 

//----------------------------------------------------------------------------//
// Some code to help set matrices with patterns 
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Create a factor loading pattern matrix based on a list of strings 
//----------------------------------------------------------------------------//

capture program drop loadingpattern
prog loadingpattern, rclass
	
	syntax namelist 
	tempname tempframe tempmat loadingmat 
	
	quietly {
		foreach rowname of local namelist {
			local rownames `"`rownames'"`rowname'""'
		}		
		mata: `tempmat' = list_to_mat("`namelist'")
		mkf `tempframe'
		frame `tempframe' {
				getmata (v*) = `tempmat'
				foreach v of var * {
					tab `v', gen(`v'_)
					levelsof `v', local(list)
					local newlist
					foreach name of local list {
						local newlist : list newlist | name
					}
					ds `v'_*
					rename (`r(varlist)') (`newlist')
					drop `v'
				}
				mkmat * , matrix(`loadingmat')
				ds 
				local list `r(varlist)'
				foreach colname of local list {
					local colnames `"`colnames'"`colname'""'
				}		
				
				mat rownames `loadingmat' = `rownames'
				mat colnames `loadingmat' = `colnames'
			}
		return matrix loadingmat = `loadingmat'
		return local vars `list'
	}
end

mata
string matrix list_to_mat(string scalar vars) 
{
	string rowvector 	varvec
	string matrix 		varmat
	real scalar 		factors,i
	
	varvec = tokens(vars)
	factors = rowsum(ustrwordcount(ustrsplit(varvec[1],"_")))
	varmat = J(cols(varvec),factors,"")
	for (i=1;i<=length(varvec);i++) {
		varmat[i,1...] = ustrsplit(varvec[i],"_")
	}
	return(varmat)
}
end

//----------------------------------------------------------------------------//
// Fills paramater matrices with randomly generated numbers within
//  specified ranges 
//----------------------------------------------------------------------------//
mata
real matrix edit_param_mat(real matrix loadingmat,
					real scalar set1to1, // option finds the first nonzero element to replace as 1 and replaces subsequent non-zero elements with random numbers
					real scalar lower, 
					real scalar upper) 
{						 			
	real scalar i,j,first_nonzero_found
	real matrix loadingmat_edited
	
	loadingmat_edited = loadingmat
	for (i=1;i<=cols(loadingmat);i++){
		first_nonzero_found = 0
		for (j=1; j<=rows(loadingmat); j++) {
			if (loadingmat[j,i] != 0) {
				if (first_nonzero_found == 0 & set1to1 == 1) {
					first_nonzero_found = 1
					loadingmat_edited[j,i] = 1
				} 
				else {
					loadingmat_edited[j,i] = runiform(1,1,lower,upper)
				}
			}
		}
	}
	return(loadingmat_edited)
}			  
end

//----------------------------------------------------------------------------//
// sort_desc_length() sort and insert algorithm
//----------------------------------------------------------------------------//
mata 
string matrix sort_desc_length(string matrix M)
{
	real   scalar 		i,j
	string scalar		temp 
	string matrix 		res
	
	res = M 
	for (i=2;i<=length(M);i++){
		temp = res[i]
		j = i - 1
		while (j >= 1 & strlen(temp) > strlen(res[j])) {
			res[j + 1] = res[j]
			j--
			if (j < 1) break 
		}
		res[j+1] = temp 
	}
	return(res)
}
end 
//----------------------------------------------------------------------------//
// Combinations and Permutations 
//----------------------------------------------------------------------------//

//-----------------------//
// Combinations 
//-----------------------//
capture program drop combn
program combn, rclass 
	syntax , MATNAME(string) N(integer) K(integer)
	tempname tempframe a
	
	if `k' > `n' {
		di "k must be less than or equal to n"
		exit
	}
	
	forvalues i =  1/`n' {
		mat `a' = (nullmat(`a') \ `i')
	}
	mkf `tempframe'
	frame `tempframe' {
		svmat `a', name(`matname')
		local k_1 = `k' - 1
		forvalues i = 1/`k_1' {
			local j = `i' + 1
			gen `matname'`j' = `matname'`i'
		}
		fillin *
		drop _fillin
		qui ds 
		local vars `r(varlist)'
		tokenize `vars'
		local n: list sizeof local(vars)
		local i = 1
		local j = 1
		while (`j' < `n') {
			local j = `i' + 1
			drop if ``i'' >= ``j''
			local++i
		}
		mkmat * , matrix(`matname') rowprefix(comb)
	}
	return matrix `matname' = `matname'
end 

//-----------------------//
// Permutations 
//-----------------------//

capture program drop permin
program permin, rclass 
	syntax , MATNAME(string) N(integer) K(integer)
	tempname tempframe a
	
	
	if `k' > `n' {
		di "k must be less than or equal to n"
		exit 
	}
	
	forvalues i =  1/`n' {
		mat `a' = (nullmat(`a') \ `i')
	}
	mkf `tempframe'
	frame `tempframe' {
		svmat `a', name(`matname')
		local k_1 = `k' - 1
		forvalues i = 1/`k_1' {
			local j = `i' + 1
			gen `matname'`j' = `matname'`i'
		}
		fillin *
		drop _fillin
		qui ds 
		local vars `r(varlist)'
		tokenize `vars'
		local i = 2
		local j = `k' - 1
		forvalues g = 1/`j' {
			forvalues x = 1/`g' {
				capture drop if ``x'' == ``i''
			}
			local++i	
		}
		mkmat * , matrix(`matname') rowprefix(perm)
	}
	return matrix `matname' = `matname'
end 


//----------------------------------------------------------------------------//
// Defining Error Structures  
//----------------------------------------------------------------------------//

cap prog drop errorstructure 
program errorstructure , rclass 

	syntax varlist , stub(string) pattern1(numlist min=2 max=2) [pattern2(numlist min=2 max=2) 	pattern3(numlist min=2 max=2) pattern4(numlist min=2 max=2)]
					  
	tempname errormat
	
	matrix  `errormat' = J(`: word count `varlist'',`: word count `varlist'',0)
	forvalues i = 1/`: word count `varlist'' { 
		forvalues j = 1/`: word count `varlist'' {
			if (strpos("`:word `i' of `varlist''","`stub'") > 0 &  strpos("`:word `j' of `varlist''","`stub'") > 0) {
				if ("`pattern4'" != "") {
					local row1  = substr("`:word `i' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					local col1  = substr("`:word `j' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					local row2  = substr("`:word `i' of `varlist''",`:word 1 of `pattern2'',`:word 2 of `pattern2'')
					local col2  = substr("`:word `j' of `varlist''",`:word 1 of `pattern2'',`:word 2 of `pattern2'')
					local row3  = substr("`:word `i' of `varlist''",`:word 1 of `pattern3'',`:word 2 of `pattern3'')
					local col3  = substr("`:word `j' of `varlist''",`:word 1 of `pattern3'',`:word 2 of `pattern3'')
					local row4  = substr("`:word `i' of `varlist''",`:word 1 of `pattern4'',`:word 2 of `pattern4'')
					local col4  = substr("`:word `j' of `varlist''",`:word 1 of `pattern4'',`:word 2 of `pattern4'')
					if ("`row1'" == "`col1'" | "`row2'" == "`col2'" | "`row3'" == "`col3'" | "`row4'" == "`col4'")  mat `errormat'[`i',`j'] = .
				}
				else if ("`pattern3'" != "") {
					local row1  = substr("`:word `i' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					local col1  = substr("`:word `j' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					local row2  = substr("`:word `i' of `varlist''",`:word 1 of `pattern2'',`:word 2 of `pattern2'')
					local col2  = substr("`:word `j' of `varlist''",`:word 1 of `pattern2'',`:word 2 of `pattern2'')
					local row3  = substr("`:word `i' of `varlist''",`:word 1 of `pattern3'',`:word 2 of `pattern3'')
					local col3  = substr("`:word `j' of `varlist''",`:word 1 of `pattern3'',`:word 2 of `pattern3'')
					if ("`row1'" == "`col1'" | "`row2'" == "`col2'" | "`row3'" == "`col3'")  mat `errormat'[`i',`j'] = .
				}
				else if ("`pattern2'" != "") {
					local row1  = substr("`:word `i' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					local col1  = substr("`:word `j' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					local row2  = substr("`:word `i' of `varlist''",`:word 1 of `pattern2'',`:word 2 of `pattern2'')
					local col2  = substr("`:word `j' of `varlist''",`:word 1 of `pattern2'',`:word 2 of `pattern2'')
					if ("`row1'" == "`col1'" | "`row2'" == "`col2'")  mat `errormat'[`i',`j'] = .
				}
				else {
					local row1  = substr("`:word `i' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					local col1  = substr("`:word `j' of `varlist''",`:word 1 of `pattern1'',`:word 2 of `pattern1'')
					if ("`row1'" == "`col1'")   mat `errormat'[`i',`j'] = .
				}
			}
			else {
				if ("`:word `i' of `varlist''" == "`:word `j' of `varlist''") mat `errormat'[`i',`j'] = .
			}
		}
	}

	return matrix errormat = `errormat'
	
end 





