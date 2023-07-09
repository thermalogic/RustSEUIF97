
/*-----------------------------------------------------------------
the coefficients and   exponents of the dimensionless Helmholtz free energy
                   for Region 3 (Table 30, page 30)
     Speciphic Helmholtz free energy.
       * tau :dimensionless temperature [K]
       * delta: dimensionless density [kg/m3]

Author: Maohua Cheng
-------------------------------------------------------------------------------*/
use crate::algo::fast_ipower::sac_pow;
use crate::common::constant::*;

//  Initialize coefphicients and exponents for region 3
pub const n1:f64 = 0.10658070028513e1;
	
pub const IJn:[IJnData;39] = [
	IJnData{I:0, J:0,n: -0.15732845290239E+02},
	IJnData{I:0, J:1, n: 0.20944396974307E+02},
	IJnData{I:0, J:2, n: -0.76867707878716E+01},
	IJnData{I:0, J:7, n: 0.26185947787954E+01},
	IJnData{I:0, J:10, n: -0.28080781148620E+01},

	IJnData{I:0, J:12, n: 0.12053369696517E+01},
	IJnData{I:0, J:23, n: -0.84566812812502E-02},
	IJnData{I:1, J:2, n: -0.12654315477714E+01},
	IJnData{I:1, J:6, n: -0.11524407806681E+01},
	IJnData{I:1, J:15, n: 0.88521043984318E+00},

	IJnData{I:1, J:17, n: -0.64207765181607E+00},
	IJnData{I:2, J:0, n: 0.38493460186671E+00},
	IJnData{I:2, J:2, n: -0.85214708824206E+00},
	IJnData{I:2, J:6, n: 0.48972281541877E+01},
	IJnData{I:2, J:7, n: -0.30502617256965E+01},

	IJnData{I:2, J:22, n: 0.39420536879154E-01},
	IJnData{I:2, J:26, n: 0.12558408424308E+00},
	IJnData{I:3, J:0, n: -0.27999329698710E+00},
	IJnData{I:3, J:2, n: 0.13899799569460E+01},
	IJnData{I:3, J:4, n: -0.20189915023570E+01},

	IJnData{I:3, J:16, n: -0.82147637173963E-02},
	IJnData{I:3, J:26, n: -0.47596035734923E+00},
	IJnData{I:4, J:0, n: 0.43984074473500E-01},
	IJnData{I:4, J:2, n: -0.44476435428739E+00},
	IJnData{I:4, J:4, n: 0.90572070719733E+00},

	IJnData{I:4, J:26, n: 0.70522450087967E+00},
	IJnData{I:5, J:1, n: 0.10770512626332E+00},
	IJnData{I:5, J:3, n: -0.32913623258954E+00},
	IJnData{I:5, J:26, n: -0.50871062041158E+00},
	IJnData{I:6, J:0, n: -0.22175400873096E-01},

	IJnData{I:6, J:2, n: 0.94260751665092E-01},
	IJnData{I:6, J:26, n: 0.16436278447961E+00},
	IJnData{I:7, J:2, n: -0.13503372241348E-01},
	IJnData{I:8, J:26, n: -0.14834345352472E-01},
	IJnData{I:9, J:2, n: 0.57922953628084E-03},

	IJnData{I:9, J:26, n: 0.32308904703711E-02},
	IJnData{I:10, J:0, n: 0.80964802996215E-04},
	IJnData{I:10, J:1, n: -0.16557679795037E-03},
	IJnData{I:11, J:26, n: -0.44923899061815E-04}];

pub fn phi_reg3(tau:f64,delta:f64)->f64
// Fundamental equation for region 3
{
	let mut result:f64  = n1 * delta.ln();
	for k in 0..39{
		result += IJn[k].n * sac_pow(delta, IJn[k].I) * sac_pow(tau, IJn[k].J);
    	//result += IJn[k].n * delta.powi(IJn[k].I) * tau.powi(IJn[k].J);
	}
	return result;
}

pub fn phi_delta_reg3(tau:f64,delta:f64)->f64
// first derivative in delta of fundamental equation for region 3
{
	let mut result:f64= n1 / delta;
	for k in 0..39
	{
	   result += IJn[k].n * IJn[k].I as f64* sac_pow(delta, IJn[k].I - 1) * sac_pow(tau, IJn[k].J);
	   //result += IJn[k].n * IJn[k].I as f64* delta.powi(IJn[k].I - 1) * tau.powi(IJn[k].J);
	}
	return  result;
}

pub fn phi_deltadelta_reg3(tau:f64,delta:f64)->f64
// Second derivative in delta of fundamental equation for region 3
{
	let mut result:f64 = -n1 / delta / delta;
	for k in 0..39
    {
		 result += IJn[k].n *(IJn[k].I * (IJn[k].I - 1)) as f64 * sac_pow(delta, IJn[k].I - 2) * sac_pow(tau, IJn[k].J);
		// result += IJn[k].n *(IJn[k].I * (IJn[k].I - 1)) as f64 * delta.powi(IJn[k].I - 2) * tau.powi(IJn[k].J);
	 
	}
	return result;
}

pub fn phi_tau_reg3(tau:f64,delta:f64)->f64
// phirst derivative in tau of fundamental equation for region 3
{
	let mut result:f64 = 0.0;
	for k in 0..39
    {	
		result += IJn[k].n * sac_pow(delta, IJn[k].I) * IJn[k].J as f64 * sac_pow(tau, IJn[k].J - 1);
		//result += IJn[k].n * delta.powi(IJn[k].I) * IJn[k].J as f64 * tau.powi(IJn[k].J - 1);

    }
	return result;
}

pub fn phi_tautau_reg3(tau:f64,delta:f64)->f64
// Second derivative in tau of fundamental equation for region 3
{
	let mut result:f64 = 0.0;
	for k in 0..39
	{
		result += IJn[k].n * sac_pow(delta, IJn[k].I) * (IJn[k].J * (IJn[k].J - 1)) as f64 * sac_pow(tau, IJn[k].J - 2);
		//result += IJn[k].n * delta.powi(IJn[k].I) * (IJn[k].J * (IJn[k].J - 1)) as f64 * tau.powi(IJn[k].J - 2);
	}
	return result;
}

pub fn phi_deltatau_reg3(tau:f64,delta:f64)->f64
// Second derivative in delta and tau of fundamental equation for region 3
{
	let mut result:f64 = 0.0;
	for k in 0..39{
		result += IJn[k].n * IJn[k].I as f64 * sac_pow(delta, IJn[k].I - 1) * IJn[k].J as f64* sac_pow(tau, IJn[k].J - 1);
		//result += IJn[k].n * IJn[k].I as f64 *delta.powi(IJn[k].I - 1) * IJn[k].J as f64* tau.powi(IJn[k].J - 1);
	}
	return result;
}
