/*
canopy_et.c
A single-function treatment of canopy evaporation and transpiration
fluxes.  

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int canopy_et(const metvar_struct* metv, const epconst_struct* epc, 
epvar_struct* epv, wflux_struct* wf, const cstate_struct* cs,int verbose)
{
	int ok=1;
	double gl_bl, gl_c, gl_s_sun, gl_s_shade;
	double gl_e_wv, gl_t_wv_sun, gl_t_wv_shade, gl_sh;
	double gc_e_wv, gc_sh;
	double gm;
	double tday;
	double tmin;
	double dayl;
	double vpd,vpd_open,vpd_close;
	double psi,psi_open,psi_close,psi_leaf_yesterday;
	double max_mesophyll_path, min_mesophyll_path;
	double Kl_max, Kl_min, Kl;  /* leaf hydraulic conductivity mmol/m2/s/MPa  */
	double Kl_psi_slope, Kl_psi_intercept;
	double psi_p_50,s_psi_p_50, psi_refill, soil_psi_today, leaf_psi_yesterday, leaf_psi_today; /*WJM 0821*/
	double livestem_leaf_c_ratio, new_livestem_leaf_c_ratio, r_max;//WJM 1021
	double m_ppfd_sun, m_ppfd_shade;
	double m_tmin, m_psi, m_co2, m_vpd, m_final_sun, m_final_shade, m_Kl;
	double m_psi_x_yesterday, m_psi_x,PRC;
//	bool allow_refill;
	double proj_lai;
	double canopy_w;
	double gcorr;
	
	double e, cwe, t, trans, trans_sun, trans_shade, e_dayl,t_dayl;
	pmet_struct pmet_in;

	/* assign variables that are used more than once */
	tday =      metv->tday;
	tmin =      metv->tmin;
	vpd =       metv->vpd;
	dayl =      metv->dayl;
	psi =       epv->psi;
	psi_leaf_yesterday = epv->leaf_psi;
	proj_lai =  epv->proj_lai;
	canopy_w =  wf->prcp_to_canopyw;
	psi_open =  epc->psi_open;
	psi_close = epc->psi_close;
	vpd_open =  epc->vpd_open;
	vpd_close = epc->vpd_close;
        max_mesophyll_path = epc->max_mesophyll_path;
        min_mesophyll_path  = epc->min_mesophyll_path;
        psi_p_50 =  epc->psi_p_50; /*WJM 0821*/
        s_psi_p_50 =  epc->s_psi_p_50;
        psi_refill = epc->psi_refill;
        
        m_psi_x_yesterday = epv->m_psi_x;
       // allow_refill = epv->allow_refill;
       
       /*WJM 1021 find out how much live stem there is relative to the prescribed allometry
       scale  (r_max is the maximum daily m_psi_x)
       capped sapwood allocation when sapwood c is twice the prescribed allometry?*/
       livestem_leaf_c_ratio = cs->livestemc/cs->leafc;
       new_livestem_leaf_c_ratio = epc->alloc_newstemc_newleafc*epc->alloc_newlivewoodc_newwoodc;
       r_max = livestem_leaf_c_ratio/new_livestem_leaf_c_ratio;
       

	/* temperature and pressure correction factor for conductances */
	/* gcorr = pow((metv->tday+273.15)/293.15, 1.75) * 101300/metv->pa; */
	gcorr = pow((metv->tday+273.15)/293.15, 1.75) * P_STD/metv->pa;

	/* calculate leaf- and canopy-level conductances to water vapor and
	sensible heat fluxes */

	/* leaf boundary-layer conductance, set to zero if lai is zero WJM 0921 */
		gl_bl = epc->gl_bl * gcorr;

	
	
	/* leaf cuticular conductance, set to zero if lai is zero WJM 0921 */
		gl_c = epc->gl_c * gcorr;

	
	
	/* leaf stomatal conductance: first generate multipliers, then apply them
	to maximum stomatal conductance */
	/* calculate stomatal conductance radiation multiplier:
	*** NOTE CHANGE FROM BIOME-BGC CODE ***
	The original Biome-BGC formulation follows the arguments in 
	Rastetter, E.B., A.W. King, B.J. Cosby, G.M. Hornberger, 
	   R.V. O'Neill, and J.E. Hobbie, 1992. Aggregating fine-scale
	   ecological knowledge to model coarser-scale attributes of 
	   ecosystems. Ecological Applications, 2:55-70.

	gmult->max = (gsmax/(k*lai))*log((gsmax+rad)/(gsmax+(rad*exp(-k*lai))))

	I'm using a much simplified form, which doesn't change relative shape
	as gsmax changes. See Korner, 1995.
	*/
	/* photosynthetic photon flux density conductance control */
	m_ppfd_sun = metv->ppfd_per_plaisun/(PPFD50 + metv->ppfd_per_plaisun);
	m_ppfd_shade = metv->ppfd_per_plaishade/(PPFD50 + metv->ppfd_per_plaishade);
	

	/* leaf hydraulic conductance effect, mmol/m2/s/MPa  from
	Brodribb et al. 2007.  Leaf maximum photosynthetic rate and venation
	linked by hydraulics. Plant Physiology, Vol. 144, pp. 1890–1898 */
        Kl_max = 12674 * pow(min_mesophyll_path , -1.26);
 	Kl_min = 12674 * pow(max_mesophyll_path , -1.26);

        Kl_psi_slope = (Kl_max - Kl_min) / (psi_open - psi_close);
	Kl_psi_intercept = Kl_max - (Kl_psi_slope * psi_open);
	Kl = (Kl_psi_slope * psi) + Kl_psi_intercept;

	if (Kl < Kl_min){
		Kl = Kl_min;
                m_Kl = 0.000001;
	}

	else if (Kl > Kl_max){
		Kl = Kl_max;
		m_Kl = 1.0;
	}

	else {
	    	m_Kl =  Kl / Kl_max;
		
	}

	/* removed to test impact of leaf psi definition */
	psi = psi - (1/Kl) * epc->gl_smax * 40.0;

    //save yesterday's leaf psi for calculating m_psi_x
	leaf_psi_yesterday = epv->leaf_psi;
        epv->leaf_psi = psi;

	/* soil-leaf water potential multiplier */
	 if (psi > psi_open){    /* no water stress, this was unbracketed WJM 0921 */
		m_psi = 1.0;
	}
	else
	if (psi <= psi_close)   /* full water stress */
	{
		 m_psi = 0.0;

	}
	else                   /* partial water stress */
		m_psi = (psi_close - psi) / (psi_close - psi_open);

	/*xylem embolism multiplier:
	PRC=1-PLC is value of normal probability distribution with mu and s
	at any particular leaf psi
	re-wetting only allowed when leaf_psi goes above threshold
	'extra' variables for mathematical clarity in multiplier calculation
	WJM 0821*/
	/*normal distr version*/
	double mu = psi_p_50;
	double s = s_psi_p_50;
	double m_psi_x_refill = 0.5 * erfc((mu - psi_refill)/(sqrt(2.0)*s));
	double m_psi_x_soil = 0.5 * erfc((mu - epv->psi)/(sqrt(2.0)*s));
	/*gamma distribtution version*/
// 	double k = vulnerability_shape;
// 	double theta = vulnerability_scale;
// 	double m_psi_x_refill = 0.5 * erfc((mu - psi_refill)/(sqrt(2.0)*s));
// 	double m_psi_x_soil = 0.5 * erfc((mu - epv->psi)/(sqrt(2.0)*s));
	/*save yesterdays m_x*/
	m_psi_x_yesterday = epv->m_psi_x;
	//today's PRC = 1-PLC if refilling is allowed or if yesterday's was higher
	//if refilling goes with leaf psi
	PRC = 0.5 * erfc((mu - psi)/(sqrt(2.0)*s));
	//if plc also goes with soil psi
	//m_psi_x_today = m_psi_x_soil;
	//worry about underflow with erfc? guaranteed underflow if argument is g.t. 26.55
	
	
/*since m_psi_x is increasing, psi > psi_yesterday -> m_psi_x > m_psi_x_yesterday... etc*/

// refill based on psi not m_psi
// 		if(epv->psi > psi_refill) { //allow refill; cannot use m_psi* here because m_psi goes to 1 before psi_refill for several taxa
// 				m_psi_x = m_psi_x_today;//+0.0001; //passive 'physical' refilling minimum?
// 			//	printf("refill ;");
// 			} else { //not above refilling threshold
// 				if (m_psi_x_today < m_psi_x_yesterday){ //drier than yesterday
// 					m_psi_x = m_psi_x_today;
// 				//	printf("not wetter ;");
// 				} else {//wetter than yesterday, but still below refilling threshold
// 					m_psi_x = m_psi_x_yesterday;//+0.0001; //passive 'physical' refilling minimum?
// 				//	printf("hysteresis ;");
// 				}
// 			}
// 			
// 			if(m_psi_x > max_daily_m_psi_x){
// 				m_psi_x = max_daily_m_psi_x;
// 			}
			//maybe cap at prescribed allometry here

		if(epv->psi > psi_refill) { //allow refill; cannot use m_psi* here because m_psi goes to 1 before psi_refill for several taxa
				m_psi_x = PRC;//+0.0001; //passive 'physical' refilling minimum?
				//printf("refill ;");
			} else { //not above refilling threshold, take the lower of PRC and m_psi_x_yesterday
				if (PRC < m_psi_x_yesterday){ //drier than yesterday, or growth recovery overshot water limitation
					m_psi_x = PRC;
					//printf("not wetter ;");
				} else {//wetter than yesterday, but still below refilling threshold
					m_psi_x = m_psi_x_yesterday;//+0.0001; //passive 'physical' refilling minimum?
					//printf("hysteresis ;");
				}
			}

// cap stomatal conductance according to allometric ratio found above (65)
// 			if(m_psi_x > r_max){
// 				m_psi_x = r_max;
// 				//printf("m_psi_x capped");
// 			}




	/* CO2 multiplier */
	m_co2 = 1.0;

	/* freezing night minimum temperature multiplier */
	if (tmin > 0.0)        /* no effect */
		m_tmin = 1.0;
	else
	if (tmin < -8.0)       /* full tmin effect */
		m_tmin = 0.0;
	else                   /* partial reduction (0.0 to -8.0 C) */
		m_tmin = 1.0 + (0.125 * tmin);
	
	/* vapor pressure deficit multiplier, vpd in Pa */
	if (vpd < vpd_open)    /* no vpd effect */
		m_vpd = 1.0;
	else
	if (vpd > vpd_close)   /* full vpd effect */
		m_vpd = 0.0;
	else                   /* partial vpd effect */
		m_vpd = (vpd_close - vpd) / (vpd_close - vpd_open);
	

        /* printf("%f\t%f\t%f\t%f\t%f\n", Kl_psi_slope,Kl_psi_intercept,Kl_max,Kl,m_Kl); */

	/* apply all multipliers to the maximum stomatal conductance */
	/* changing for using leaf conductivity instead  JDW 2/2/18  */

	//m_final_sun = m_ppfd_sun * m_psi * m_co2 * m_tmin * m_vpd; /*Kl built into m_psi*/
	m_final_sun = m_ppfd_sun * m_psi * m_psi_x * m_co2 * m_tmin * m_vpd; /*Kl built into m_psi*/
    //m_final_sun = m_ppfd_sun * m_psi * m_Kl * m_co2 * m_tmin * m_vpd; /* both m_psi and m_Kl*/
	//m_final_sun = m_ppfd_sun * m_Kl * m_co2 * m_tmin * m_vpd; /* just m_Kl*/
	if (m_final_sun < 0.00000001) m_final_sun = 0.00000001;
	
	//m_final_shade = m_ppfd_shade * m_psi * m_co2 * m_tmin * m_vpd;
	m_final_shade = m_ppfd_shade * m_psi * m_psi_x * m_co2 * m_tmin * m_vpd;
        /* m_final_shade = m_ppfd_shade * m_psi * m_Kl * m_co2 * m_tmin * m_vpd; */
        /* m_final_shade = m_ppfd_shade * m_Kl * m_co2 * m_tmin * m_vpd; */
	if (m_final_shade < 0.00000001) m_final_shade = 0.00000001;

	if(proj_lai >  0.0)
	{
		gl_s_sun = epc->gl_smax * m_final_sun * gcorr;
		gl_s_shade = epc->gl_smax * m_final_shade * gcorr;
	} else {
		gl_s_sun = 0.0;
		gl_s_shade = 0.0;
	}


 				//printf(" l_psi = %f; l_psi_y = %f ; r_mx = %f ; m_x = %f ; PRC = %f ; m_x_y = %f, gls = %f \n "\
 							,epv->leaf_psi, leaf_psi_yesterday, r_max, m_psi_x, PRC, m_psi_x_yesterday, gl_s_sun);
	/* calculate leaf-and canopy-level conductances to water vapor and
	sensible heat fluxes, to be used in Penman-Monteith calculations of
	canopy evaporation and canopy transpiration. */
	
	/* Leaf conductance to evaporated water vapor, per unit projected LAI */
	gl_e_wv = gl_bl;

	/* Leaf conductance to transpired water vapor, per unit projected
	LAI.  This formula is derived from stomatal and cuticular conductances
	in parallel with each other, and both in series with leaf boundary
	layer conductance. */

	/* gl_t_wv_sun = (gl_bl * (gl_s_sun + gl_c)) / (gl_bl + gl_s_sun + gl_c);
	gl_t_wv_shade = (gl_bl * (gl_s_shade + gl_c)) / (gl_bl + gl_s_shade + gl_c); */


	 /* including mesophyll conductance m/s static value for paleo plants JDW 3/26/18*/
	gm = 0.0066;

	/* gl_t_wv_sun = (gl_bl / (1 + (gl_bl / (gl_c + ((gm * gl_s_sun)/(gm + gl_s_sun))))));
        gl_t_wv_shade = (gl_bl / (1 + (gl_bl / (gl_c + ((gm * gl_s_shade)/(gm + gl_s_shade)))))); */
 
 /*add cutoff to check if LAI is zero, WJM 0921*/

        gl_t_wv_sun = gl_bl * gm * (gl_s_sun + gl_c) / 
		    (gl_bl * (gl_s_sun + gl_c) + gm * gl_bl  + gm * (gl_s_sun + gl_c));  

        gl_t_wv_shade = gl_bl * gm * (gl_s_shade + gl_c) /
		    (gl_bl * (gl_s_shade + gl_c) + gm * gl_bl  + gm * (gl_s_shade + gl_c));

        


	/* Leaf conductance to sensible heat, per unit all-sided LAI */
	gl_sh = gl_bl;

	/* Canopy conductance to evaporated water vapor */
	gc_e_wv = gl_e_wv * proj_lai;
	
	/* Canopy conductane to sensible heat */
	gc_sh = gl_sh * proj_lai;
	
	cwe = trans = 0.0;
	
	/* Assign values in pmet_in that don't change */
	pmet_in.ta = tday;
	pmet_in.pa = metv->pa;
	pmet_in.vpd = vpd;
	
	/* Canopy evaporation, if any water was intercepted */
	/* Calculate Penman-Monteith evaporation, given the canopy conductances to
	evaporated water and sensible heat.  Calculate the time required to 
	evaporate all the canopy water at the daily average conditions, and 
	subtract that time from the daylength to get the effective daylength for
	transpiration. */
	if (canopy_w)
	{
		/* assign appropriate resistance and radiation for pmet_in */
		pmet_in.rv = 1.0/gc_e_wv;
		pmet_in.rh = 1.0/gc_sh;
	        pmet_in.irad = metv->swabs;
                                           		
		/* call penman-monteith function, returns e in kg/m2/s */
		if (penmon(&pmet_in, 0, &e))
		{
			bgc_printf(BV_ERROR, "Error: penmon() for canopy evap... \n");
			ok=0;
		}
		
		/* calculate the time required to evaporate all the canopy water */
		e_dayl = canopy_w/e;
		
		if (e_dayl > dayl)  
		{
			/* day not long enough to evap. all int. water */
			trans = 0.0;    /* no time left for transpiration */
			cwe = e * dayl;   /* daylength limits canopy evaporation */
		}
		else                
		{
			/* all intercepted water evaporated */
			cwe = canopy_w;
			
			/* adjust daylength for transpiration */
			t_dayl = dayl - e_dayl;
			 
			/* calculate transpiration using adjusted daylength */
			/* first for sunlit canopy fraction */
			pmet_in.rv = 1.0/gl_t_wv_sun;
			pmet_in.rh = 1.0/gl_sh;
			pmet_in.irad = metv->swabs_per_plaisun;
			if (penmon(&pmet_in, 0, &t))
			{
				bgc_printf(BV_ERROR, "Error: penmon() for adjusted transpiration... \n");
				ok=0;
			}
			trans_sun = t * t_dayl * epv->plaisun;
			
			/* next for shaded canopy fraction */
			pmet_in.rv = 1.0/gl_t_wv_shade;
			pmet_in.rh = 1.0/gl_sh;
			pmet_in.irad = metv->swabs_per_plaishade;
			if (penmon(&pmet_in, 0, &t))
			{
				bgc_printf(BV_ERROR, "Error: penmon() for adjusted transpiration... \n");
				ok=0;
			}
			trans_shade = t * t_dayl * epv->plaishade;
			trans = trans_sun + trans_shade;
		}
	}    /* end if canopy_water */
	
	else /* no canopy water, transpiration with unadjusted daylength */
	{
		/* first for sunlit canopy fraction */
		pmet_in.rv = 1.0/gl_t_wv_sun;
		pmet_in.rh = 1.0/gl_sh;
		pmet_in.irad = metv->swabs_per_plaisun;
		if (penmon(&pmet_in, 0, &t))
		{
			bgc_printf(BV_ERROR, "Error: penmon() for adjusted transpiration... \n");
			ok=0;
		}
		trans_sun = t * dayl * epv->plaisun;
		
		/* next for shaded canopy fraction */
		pmet_in.rv = 1.0/gl_t_wv_shade;
		pmet_in.rh = 1.0/gl_sh;
		pmet_in.irad = metv->swabs_per_plaishade;
		if (penmon(&pmet_in, 0, &t))
		{
			bgc_printf(BV_ERROR, "Error: penmon() for adjusted transpiration... \n");
			ok=0;
		}
		trans_shade = t * dayl * epv->plaishade;
		trans = trans_sun + trans_shade;
	}
		
	/* assign water fluxes, all excess not evaporated goes
	to soil water compartment */
	wf->canopyw_evap = cwe;
	wf->canopyw_to_soilw = canopy_w - cwe;
	wf->soilw_trans = trans;
	
	/* assign leaf-level conductance to transpired water vapor, 
	for use in calculating co2 conductance for farq_psn() */
	epv->gl_t_wv_sun = gl_t_wv_sun;
	epv->gl_t_wv_shade = gl_t_wv_shade;
	
	/* assign verbose output variables if requested */
	if (verbose)
	{
		epv->m_ppfd_sun  = m_ppfd_sun;
		epv->m_ppfd_shade  = m_ppfd_shade;
		epv->m_psi   = m_psi;
		epv->m_co2   = m_co2;
		epv->m_tmin  = m_tmin;
		epv->m_vpd   = m_vpd;
		epv->m_Kl    = m_Kl;
		epv->m_psi_x    = m_psi_x; /*WJM 0821*/
//		epv->allow_refill = allow_refill;
		epv->m_final_sun = m_final_sun;
		epv->m_final_shade = m_final_shade;
		epv->gl_bl   = gl_bl;
		epv->gl_c    = gl_c;
		epv->gl_s_sun   = gl_s_sun;
		epv->gl_s_shade = gl_s_shade;
		epv->gl_e_wv = gl_e_wv;
		epv->gl_sh   = gl_sh;
		epv->gc_e_wv = gc_e_wv;
		epv->gc_sh   = gc_sh;
	}
	
	return (!ok);
}

int penmon(const pmet_struct* in, int out_flag,	double* et)
{
    /*
	Combination equation for determining evaporation and transpiration. 
	
    For output in units of (kg/m2/s)  out_flag = 0
    For output in units of (W/m2)     out_flag = 1
   
    INPUTS:
    in->ta     (deg C)   air temperature
    in->pa     (Pa)      air pressure
    in->vpd    (Pa)      vapor pressure deficit
    in->irad   (W/m2)    incident radient flux density
    in->rv     (s/m)     resistance to water vapor flux
    in->rh     (s/m)     resistance to sensible heat flux

    INTERNAL VARIABLES:
    rho    (kg/m3)       density of air
    CP     (J/kg/degC)   specific heat of air
    lhvap  (J/kg)        latent heat of vaporization of water
    s      (Pa/degC)     slope of saturation vapor pressure vs T curve

    OUTPUT:
    et     (kg/m2/s)     water vapor mass flux density  (flag=0)
    et     (W/m2)        latent heat flux density       (flag=1)
    */

    int ok=1;
    double ta;
    double rho,lhvap,s;
    double t1,t2,pvs1,pvs2,e,tk;
    double rr,rh,rhr,rv;
    double dt = 0.2;     /* set the temperature offset for slope calculation */
   
    /* assign ta (Celsius) and tk (Kelvins) */
    ta = in->ta;
    tk = ta + 273.15;
        
    /* calculate density of air (rho) as a function of air temperature */
    rho = 1.292 - (0.00428 * ta);
    
    /* calculate resistance to radiative heat transfer through air, rr */
    rr = rho * CP / (4.0 * SBC * (tk*tk*tk));
    
    /* resistance to convective heat transfer */
    rh = in->rh;
    
    /* resistance to latent heat transfer */
    rv = in->rv;
    
    /* calculate combined resistance to convective and radiative heat transfer,
    parallel resistances : rhr = (rh * rr) / (rh + rr) */
    rhr = (rh * rr) / (rh + rr);

    /* calculate latent heat of vaporization as a function of ta */
    lhvap = 2.5023e6 - 2430.54 * ta;

    /* calculate temperature offsets for slope estimate */
    t1 = ta+dt;
    t2 = ta-dt;

    /* calculate saturation vapor pressures at t1 and t2 */
    pvs1 = 610.7 * exp(17.38 * t1 / (239.0 + t1));
    pvs2 = 610.7 * exp(17.38 * t2 / (239.0 + t2));

    /* calculate slope of pvs vs. T curve, at ta */
    s = (pvs1-pvs2) / (t1-t2);
    
    /* calculate evaporation, in W/m^2  */
    e = ( ( s * in->irad ) + ( rho * CP * in->vpd / rhr ) ) /
    	( ( ( in->pa * CP * rv ) / ( lhvap * EPS * rhr ) ) + s );
    
    /* return either W/m^2 or kg/m^2/s, depending on out_flag */	
    if (out_flag)
    	*et = e;
    
    if (!out_flag)
    	*et = e / lhvap;
    
    return (!ok);
}
