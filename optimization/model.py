# Input parameters
r = 0.1
a1 = 0.05
a2 = 0.06
a3 = 0.01
b = 0.025

random = False # Set generation of random parameters

# Fixed geometry parameters
t1 = 0.025
t2 = 0.05
t3 = 0.015
S = 0.003
V1 = 5e-3
V2 = 6e-3
delta = 0.005

# Limits of input parameters
r_min = 0.05
r_max = 0.15
a1_min = 0.01
a1_max = 0.05
a2_min = 0.01
a2_max = 0.05
a3_min = 0.01
a3_max = 0.05
b_min = 0.01
b_max = 0.05

# Material parameters
mu_r_Fe = 500
gamma_Fe = 7.2452e+06
lamda_Fe = 7.1677e+01
rho_Fe = 7.8323e+03
cp_Fe = 4.7741e+02

lamda_Ceramics = 0.23
rho_Ceramics = 2140
cp_Ceramics = 1010
alpha_Ceramics = 5

lamda_Air = 10
rho_Air = 1.205
cp_Air = 1005

# Problem parameters
J_ext = 1.0e7
f = 50.0
t_step = 60.0
t_total = 5*60.0
T_ext = 20.0
T_init = T_ext

import math
import random as rnd

# Model of magnetic field
def magnetic_model(a, adaptivity="disabled", adaptivity_steps = 10, adaptivity_tolerance = 3):
	newdocument("Model of magnetic field", "axisymmetric", "magnetic", 1, 2, adaptivity, adaptivity_steps, adaptivity_tolerance, f, "harmonic", 1.0, 1.0, 0.0)

	addboundary("A = 0", "magnetic_vector_potential", 0, 0)
	addmaterial("Fe", 0, 0, mu_r_Fe, gamma_Fe, 0, 0, 0, 0, 0)
	addmaterial("Fe (skin depth)", 0, 0, mu_r_Fe, gamma_Fe, 0, 0, 0, 0, 0)
	addmaterial("Ceramics", 0, 0, 1, 0, 0, 0, 0, 0, 0)
	addmaterial("Cu", J_ext, 0, 1, 0, 0, 0, 0, 0, 0)
	addmaterial("Air", 0, 0, 1, 0, 0, 0, 0, 0, 0)

	# Charge
	addedge(0, 0, r-3.0*a, 0)
	addedge(r-3.0*a, 0, r, 0)
	addedge(r-3.0*a, 0, r-3.0*a, h1)
	addedge(r, 0, r, h1)
	addedge(r, h1, r-3.0*a, h1)
	addedge(r-3.0*a, h1, 0, h1)
	addedge(0, h1, 0, 0, 0, "A = 0")
	addlabel(r/2.0, (h2+t3)/2.0, 0, 3, "Fe")
	addlabel(r-3*a/2.0, (h2+t3)/2.0, 1e-6, 0, "Fe (skin depth)")

	# Ceramic lid
	addedge(0, h2, r, h2)
	addedge(r, h2, r+t1, h2)
	addedge(r+t1, h2, r+t1, h2+t3)
	addedge(r+t1, h2+t3, 0, h2+t3)
	addedge(0, h2+t3, 0, h2, 0, "A = 0")
	addlabel(r/2.0, h2+t3/2.0, 0, 0, "Ceramics")

	# Ceramic crucible
	addedge(r, h2, r, h1)
	addedge(r+t1, h2, r+t1, -t2)
	addedge(r+t1, -t2, 0, -t2)
	addedge(0, -t2, 0, 0, 0, "A = 0")
	addedge(0, h1, 0, h2, 0, "A = 0")
	addlabel(r/2.0, -t2/2.0, 0, 0, "Ceramics")

	# Inductors
	addrect(r+t1+delta, -t2+a1, b, a2, "none", "Cu")
	addrect(r+t1+delta, -t2+a1+a2+a3, b, a4, "none", "Cu")

	# Save image
	#zoombestfit()
	#saveimage("./Variants/r_%f-a1_%f-a2_%f_a3_%f-b_%f.png" % (r, a1, a2, a3, b))

	radius = 10*t2+h2+t3
	addedge(0, 0, 0, (h2+t3)/2.0-radius, 0, "A = 0")
	addedge(0, (h2+t3)/2.0-radius, radius, (h2+t3)/2.0, 90, "A = 0")
	addedge(radius, (h2+t3)/2.0, 0, (h2+t3)/2.0+radius, 90, "A = 0")
	addedge(0, (h2+t3)/2.0+radius, 0, h2+t3, 0, "A = 0")
	addlabel(radius-0.01, (h2+t3)/2.0, 0, 0, "Air")
	addlabel(r/2.0, h1+(h2-h1)/2.0, 0, 0, "Air")

	zoombestfit()
	solve()

	# Calculation of losses
	Q = []
	volume = volumeintegral(0, 1)
	Q.append(volume["Pj"])
	volume = volumeintegral(0)
	Q.append(volume["Pj"]/volume["V"])
	volume = volumeintegral(1)
	Q.append(volume["Pj"]/volume["V"])

	return Q

# Model of temperature field
def heat_model(a, Q, adaptivity="disabled", adaptivity_steps = 10, adaptivity_tolerance = 3):
	newdocument("Model of temperature field", "axisymmetric", "heat", 1, 2, adaptivity, adaptivity_steps, adaptivity_tolerance, 0, "transient", t_step, t_total, T_init)

	addboundary("Symmetry", "heat_heat_flux", 0, 0, 0)
	addboundary("Convection", "heat_heat_flux", 0, alpha_Ceramics, T_ext)
	addmaterial("Fe", Q[1], lamda_Fe, rho_Fe, cp_Fe)
	addmaterial("Fe (skin depth)", Q[2], lamda_Fe, rho_Fe, cp_Fe)
	addmaterial("Ceramics", 0, lamda_Ceramics, rho_Ceramics, cp_Ceramics)
	addmaterial("Air", 0, lamda_Air, rho_Air, cp_Air)

	# Charge
	addedge(0, 0, r-3.0*a, 0)
	addedge(r-3.0*a, 0, r, 0)
	addedge(r-3.0*a, 0, r-3.0*a, h1)
	addedge(r, 0, r, h1)
	addedge(r, h1, r-3.0*a, h1)
	addedge(r-3.0*a, h1, 0, h1)
	addedge(0, h1, 0, 0, 0, "Symmetry")
	addlabel(r/2.0, (h2+t3)/2.0, 0, 0, "Fe")
	addlabel(r-3*a/2.0, (h2+t3)/2.0, 0, 0, "Fe (skin depth)")

	# Ceramic lid
	addedge(0, h2, r, h2)
	addedge(r, h2, r+t1, h2)
	addedge(r+t1, h2, r+t1, h2+t3, 0, "Convection")
	addedge(r+t1, h2+t3, 0, h2+t3, 0, "Convection")
	addedge(0, h2+t3, 0, h2, 0, "Symmetry")
	addlabel(r/2.0, h2+t3/2.0, 0, 0, "Ceramics")

	# Ceramic crucible
	addedge(r, h2, r, h1)
	addedge(r+t1, h2, r+t1, -t2, 0, "Convection")
	addedge(r+t1, -t2, 0, -t2, 0, "Convection")
	addedge(0, -t2, 0, 0, 0, "Symmetry")
	addedge(0, h1, 0, h2, 0, "Symmetry")
	addlabel(r/2.0, -t2/2.0, 0, 0, "Ceramics")

	addlabel(r/2.0, h1+(h2-h1)/2.0, 0, 0, "Air")

	# Inductors
	addrect(r+t1+delta, -t2+a1, b, a2, "none", "none")
	addrect(r+t1+delta, -t2+a1+a2+a3, b, a4, "none", "none")

	radius = 10*t2+h2+t3
	addedge(0, 0, 0, (h2+t3)/2.0-radius, 0)
	addedge(0, (h2+t3)/2.0-radius, radius, (h2+t3)/2.0, 90)
	addedge(radius, (h2+t3)/2.0, 0, (h2+t3)/2.0+radius, 90)
	addedge(0, (h2+t3)/2.0+radius, 0, h2+t3, 0)
	addlabel(radius-0.01, (h2+t3)/2.0, 0, 0)
	addlabel(r/2.0, h1+(h2-h1)/2.0, 0, 0)

	zoombestfit()
	solve()

	# Calculation of criterial functions
	volume = volumeintegral(0, 1)
	Tavg = volume["T_avg"]
	V = volume["V"]

	Nr = 20
	Nz = 30
	rr = []
	zz = []
	vv = []
	eps = r/1.0e9
	for i in range(Nz-1):
		dz = (h1-2*eps)/(Nz-1)
		for j in range(Nr-1):
			dr = (r-2*eps)/(Nr-1)
			rr.append(dr/2.0 + j*dr)
			zz.append(dz/2.0 + i*dz)
			vv.append((math.pi*(dr + j*dr)**2-math.pi*(j*dr)**2)*dz/2.0)
	
	U = 0.0
	Tavg_p = 0.0
	for i in range(len(rr)):
		point = pointresult(rr[i], zz[i])
		U += ((point["T"] - Tavg)**2) * vv[i]/(V1/2.0)
		Tavg_p += (point["T"]) * vv[i]/(V1/2.0)

	#print "Test of integration: %f" % ((Tavg_p-Tavg)/Tavg)

	results = [V*rho_Fe*cp_Fe*(Tavg-T_init), U]
	return results

# Calculation of skin depth
omega = 2.0*math.pi*f
mu_Fe = mu_r_Fe*4.0*math.pi*1e-7
a = math.sqrt(2.0/(omega*gamma_Fe*mu_Fe))
#print 'a (skin depth): %f m' % (a)

N = 1
for i in range(N):
	# Generate random parameters
	if (random):
		r = (rnd.random()*r_max) + r_min
		a1 = (rnd.random()*a1_max) + a1_min
		a2 = (rnd.random()*a2_max) + a2_min
		a3 = (rnd.random()*a3_max) + a3_min
		b = (rnd.random()*b_max) + b_min

	# Calculation of other dimensions
	h1 = V1/(math.pi*r**2)
	h2 = V2/(math.pi*r**2)
	a4 = S/b-a2

	# Solution of problem
	Q = magnetic_model(a)
	results = heat_model(a, Q)
	#print 'Q (max. function): %f J' % (results[0])
	#print 'U (min. function): %f deg.' % (results[1])

	# Save result to file
	#file = open('Results.csv','a')
	#file.write("%f, %f, %f, %f, %f, %f, %f\n" % (r, a1, a2, a3, b, results[0], results[1]))
	#file = open('Results.dat','a')
	#file.write("%f %f %f %f %f %f %f\n" % (r, a1, a2, a3, b, results[0], results[1]))
	#file.close()
