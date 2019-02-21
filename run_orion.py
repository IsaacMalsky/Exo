#!/usr/bin/env python
#python script created by Phil Arras UVA, modified by HW Chen NU and then modified more by Isaac Malsky
import numpy as np
import os
import scipy
from scipy import loadtxt, optimize
import mysubsprograms as my
import sys

#Constants
msun = 1.9892e33
rsun = 6.9598e10
rearth = 6.371008e8
mjup = 1.8986e30
rjup = 6.9911e9
mearth = 5.97e27
sigma=5.67e-5
au = 1.496e13

initial_mod = "initial_planet.mod"
####################################################
#########         PARAMETERS LISTS         #########
####################################################

mpList=[10]
orbitalList=[0.25]
enFracList=[0.025]
yList = [.24]
zList = [.02]
entropyList = [9.0]
n_frac_list = [.10]


####################################################
#########        IRRAD/EVOL CONDITIONS        ######
####################################################                            
rs = 1.0                      #star radius in rsun
Teff_star = 6000                 #Host Star Temp
BA= 0.20                         #planet Bond albedo
a = 1.0                          #frac_absorbing_radius
ms = 1.0                         #host_star_mass
ec = 1e9                         #eddy coefficient
formation_time = 6e6             #Disk formation time

#I need to get rid of these but they're in the calls to mysubs
initialage = 0
maxage = 0
rf = .1                #escape_rate_reduction_factor


flux_dayside = (sigma*Teff_star**4 * (rs * rsun / orbitalList[0] / au )**2)*(1-BA)    # flux hitting planet's dayside
Teq = (flux_dayside*(1-BA)/4.0/sigma)**0.25    #equalibrium temperature

print ('Teq',Teq)
for mp in mpList:
	pre_reduce_mod = "pre_reduce_" + ".mod"
	inlist_pre_reduce = "inlist_reduce_" + str(mp)
	run_time = my.run_pre_reduce(inlist_pre_reduce, initial_mod, pre_reduce_mod, mp)

	for enFrac in enFracList:
		core_mass = mp*(1-enFrac)
		rho = my.calculate_rho(mp, enFrac)[0]
		radius1 = my.calculate_rho(mp, enFrac)[1]

		pre_core_mod = "pre_core_" + str(mp) + "_" + str(enFrac) + ".mod"
		inlist_pre_core = "inlist_pre_core_" + str(mp) + "_" + str(enFrac)
		run_time = my.run_pre_core(inlist_pre_core, pre_reduce_mod, pre_core_mod, enFrac,core_mass,rho)


		for z in zList:
			for y in yList:
				comp_mod = "comp_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
				inlist_comp = "inlist_comp_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
				run_time = my.run_comp(initial_mod,inlist_comp, comp_mod, z, y)

				corel_mod = "corel_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
				inlist_corel = "inlist_corel_"  + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
				run_time = my.run_corel(inlist_corel, comp_mod, corel_mod)

				reduce_mod = "reduce_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
				inlist_reduce = "inlist_reduce_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
				run_time = my.run_reduce(inlist_reduce, corel_mod, reduce_mod, mp)

				corem_mod = "corem_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + ".mod"
				inlist_corem = "inlist_corem_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z)
				run_time = my.run_corem(inlist_corem, reduce_mod, corem_mod, core_mass, rho)



				for entropy in entropyList:
					if (os.path.isfile('LOGS/' + corem_mod) == True):
						entropy_list, luminosity_list = loadtxt('LOGS/' + corem_mod, unpack=True, skiprows =6, usecols=[12,13])
						currentropy = entropy_list[-1]
						luminosity = luminosity_list[-1]
					else:
						currentropy = -100000000000000
						luminosity = -100000000000000


					if currentropy<float(entropy):
						heating_mod = "heating_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
						inlist_heating = "inlist_heating_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
						run_time = my.run_heating(inlist_heating, corem_mod, heating_mod, entropy, luminosity)

						remove_mod = "remove_heating_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
						remove_heating_profile = "profile_remove_heating" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy)
						inlist_remove_heating = "inlist_remove_heating_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
						run_time = my.run_remove_heating(remove_heating_profile, inlist_remove_heating, heating_mod, remove_mod)


						for orb_sep in orbitalList:
							for n_frac in n_frac_list:
								if (os.path.isfile(remove_mod) == True):

									flux_dayside = (sigma*Teff_star**4 * (rs * rsun / orb_sep / au )**2)*(1-BA)	# flux hitting planet's dayside
									teq = (flux_dayside*(1-BA)/4.0/sigma)**0.25    #equalibrium temperature
									column_depth = my.calculate_column_depth(Teq, remove_heating_profile)
									if (column_depth != -1):
										#8
										irrad_mod = "irrad_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy)+ ".mod"
										inlist_irrad = "inlist_irrad_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)
										irrad_profile = "profile_irrad" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)

										run_time = my.run_irrad(irrad_profile, inlist_irrad, remove_mod, irrad_mod, column_depth, flux_dayside)
										
										#9
										evolve_mod = "evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac) + ".mod"
										evolve_profile = "profile_evolve" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)
										inlist_evolve = "inlist_evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep) + "_" + str(entropy) + "_" + str(n_frac)

                                        #column_depth = my.calculate_column_depth(Teq, flux_dayside, irrad_profile, orb_sep)
										run_time = my.run_evolve(evolve_profile, inlist_evolve, irrad_mod, evolve_mod, n_frac, a, ms, rf, orb_sep, ec, column_depth, flux_dayside, formation_time, teq)
									else:
										pass
								else:
									pass


					else:
						cooling_mod = "cooling_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
						inlist_cooling = "inlist_cooling_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
						run_time = my.run_cooling(inlist_cooling, corem_mod, cooling_mod, entropy)

						remove_mod = "remove_" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy) + ".mod"
						remove_cooling_profile = "profile_remove_cooling" + str(mp) + "_" + str(enFrac) +"_" + str(y) + "_" + str(z) + "_" + str(entropy)
						inlist_remove_cooling = "inlist_remove_cooling_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(entropy)
						run_time = my.run_remove_cooling(inlist_remove_cooling, cooling_mod, remove_mod)


						for orb_sep in orbitalList:
							for n_frac in n_frac_list:
								if (os.path.isfile(remove_mod) == True):
									flux_dayside = (sigma*Teff_star**4 * (rs * rsun / orb_sep / au )**2)*(1-BA)	# flux hitting planet's dayside
									teq = (flux_dayside*(1-BA)/4.0/sigma)**0.25    #equalibrium temperature
									column_depth = my.calculate_column_depth(Teq, remove_cooling_profile)
									if (column_depth != -1):
										#8
										irrad_mod = "irrad_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + ".mod"
										inlist_irrad = "inlist_irrad_" + str(mp) + "_" + str(enFrac)+ "_" + str(y) + "_" + str(z) + "_" + str(orb_sep) + "_" + str(entropy)
										irrad_profile = "profile_irrad" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)

										run_time = my.run_irrad(inlist_irrad, remove_mod, irrad_mod, column_depth, flux_dayside)

										#9
										evolve_mod = "evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac) + ".mod"
										evolve_profile = "profile_evolve" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep)+ "_" + str(entropy) + "_" + str(n_frac)
										inlist_evolve = "inlist_evolve_" + str(mp) + "_" + str(enFrac) + "_" + str(y) + "_" + str(z) + "_" + str(orb_sep) + "_" + str(entropy) + "_" + str(n_frac)

                                        #column_depth = my.calculate_column_depth(Teq, flux_dayside, irrad_profile, orb_sep)
										run_time = my.run_evolve(evolve_profile, inlist_evolve, irrad_mod, evolve_mod, n_frac, a, ms, rf, orb_sep, ec, column_depth, flux_dayside, formation_time, teq)
									else:
										pass
								else:
									pass

os.system('./clear_inlists')						
