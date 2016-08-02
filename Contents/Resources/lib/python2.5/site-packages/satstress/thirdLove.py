"""
Code to compute complex Love numbers given a satellite with concentric, spherically
symmetric homogeneous layers. Two layers in contact cannot both be liquids.


Based on a method developed by Zabranova, Hanyk, and Matyska, "Matrix Pseudospectral
Method for Elastic Tides Modeling", modified to use radial equations developed by
Dahlen and Tromp, "Theoretical Global Seismology".

For and explanation of the A, B, C, D, E, F, P, and Q matrices, see Zabranova.

thirdLove will handle an arbitrary number of layers, so long as there are not two
liquid layers in contact. The inputs have been changed so that the viscoelastic 
transformation for the Lame parameters needs to happen outside the love number 
function. The inputs include both the complex Lame parameters, so viscosity is no
longer needed as an input. 
"""

import scipy
import numpy

def thirdLove(M,param,omega):
	"""
	Master function for computation of Love numbers.

	Input:
	M: (int) Number of node points per layer. Total number of nodes is 4M.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	love: (1x6 list of floats) Real and imaginary parts of second degree love 
	       numbers.

		   [Re(h2), Im(h2), Re(k2), Im(k2), Re(l2), Im(l2)]
	"""

	#Read the number of layers from the number of columns in param.
	layerNum = param.shape[1]

	#Return an error if two liquid layers are next to each other.
	for layer in range(0,layerNum-1):

		if (param[1,layer]) == 0 and (param[1,layer+1] == 0):

			raise ValueError('Two adjacent layers cannot both be liquids.')

	#Planetary parameters and constants.
	R = param[0,layerNum-1]*1000.0
	G = 6.67384e-11

	###########################################################################
	#  Node Points  ###########################################################
	###########################################################################

	#Define an empty array of node points first, to be filled. Each row contains
	#nodes for a single layer.
	rSetArray = numpy.zeros((layerNum,M))

	#Fill the matrix for each layer using a Chebyshev grid. The first layer 
	#needs a slightly different formulation because no layers precede it.
	i = numpy.arange(M,0.0,-1.0) 
	rSetArray[0,:] = (param[0,0].real*1000.0/2.0)*(numpy.cos(numpy.pi*(i-1)/(M-1)) + 1)

	for index in range(1,layerNum):

		rSetArray[index,:] = param[0,index-1].real*1000.0 + ((param[0,index].real-param[0,index-1].real)*1000.0/2.0)*(numpy.cos(numpy.pi*(i-1)/(M-1)) + 1)

	###########################################################################
	#  Create P Matrix  #######################################################
	###########################################################################

	#Apply center of core boundary conditions. 
	FTemp = FGen(M)
	P = numpy.hstack((FTemp, numpy.zeros((3,3*(layerNum-1)*M)) ))

	#Cycle through the layers.
	for layer in range(0,layerNum):

		#Create the interior of the layer. Pad with zeros for other layers.
		PInterior = interiorGen(rSetArray[layer,:],M,param,omega)

		zerBelow = numpy.zeros((3*(M-2),3*layer*M),'complex')
		zerAbove = numpy.zeros((3*(M-2),3*(layerNum-layer-1)*M),'complex')

		P = numpy.vstack((P, numpy.hstack((zerBelow, PInterior, zerAbove)) ))

		#Handle the boundary conditions between layers. The free surface 
		#behaves differently.
		if layer == (layerNum-1):

			PBound = boundFree(rSetArray[-1,:],param,omega)
			zer = numpy.zeros((3,3*(layerNum-1)*M),'complex')

			P = numpy.vstack((P, numpy.hstack((zer, PBound)) ))

		else:

			#Apply different boundary conditions depending on the nature of the
			#boundary.

			#Solid-liquid:
			if (param[1,layer] == 1) and (param[1,layer+1] == 0):

				PBound = boundSL(rSetArray[layer,:],rSetArray[layer+1,:],M,param,omega)

			#Liquid-solid
			elif (param[1,layer] == 0) and (param[1,layer+1] == 1):

				PBound = boundLS(rSetArray[layer,:],rSetArray[layer+1,:],M,param,omega)

			#Solid-solid
			else:

				PBound = boundSS(rSetArray[layer,:],rSetArray[layer+1,:],M,param,omega)

			zerBelow = numpy.zeros((6,3*layer*M),'complex')
			zerAbove = numpy.zeros((6,3*(layerNum-layer-2)*M),'complex')

			P = numpy.vstack((P, numpy.hstack((zerBelow, PBound, zerAbove)) ))

	#Generate Q matrix.
	Q = QGen(rSetArray[-1,-1],layerNum,M,param)

	###########################################################################
	#  Matrix Conditioning  ###################################################
	###########################################################################

	#The following will take the logarithm of zeros, but the resulting infs are
	#fine. Temporarily adjust the error flag settings so they don't appear.
	numpy.seterr(divide='ignore')

	#The following makes the linear system described in the P matrix better
	#conditioned, so it can be solved more accurately. The matrix elements
	#are scaled to be on the order of 10^3 at largest. All but the last row
	#combine to zero, so can be scaled without affecting results.
	for index in range(0,3*layerNum*M-1):

		power = numpy.max(numpy.floor(numpy.log10(numpy.abs(P[index,:]))))
		P[index,:] = P[index,:]/pow(10.0,power-3)

	#Scale the last (non-homogeneous) equation and last element of Q to 
	#preserve the results.
	power = numpy.max(numpy.floor(numpy.log10(numpy.abs(P[-1,:]))))
	P[-1,:] = P[-1,:]/pow(10.0,power-3)
	Q[-1,0] = Q[-1,0]/pow(10.0,power-3)

	#Scale Q again for better conditioning. This time also scales the results.
	power = numpy.max(numpy.floor(numpy.log10(numpy.abs(Q[-1,0]))))
	Q[-1,0] = Q[-1,0]/pow(10.0,power-3)

	#Restore error flag settings.
	numpy.seterr(divide='warn')

	###########################################################################
	#  Solve System  ##########################################################
	###########################################################################

	y = numpy.linalg.solve(P, Q)

	#Real and imaginary love number components, cast as floats.
	loveh2Re = float((pow(10.0,power-3)*y[(3*layerNum-2)*M-1]/R).real)
	loveh2Im = float((pow(10.0,power-3)*y[(3*layerNum-2)*M-1]/R).imag)
	lovek2Re = float((-1.0*pow(10.0,power-3)*y[3*layerNum*M-1]/(R*gravCompute(R,param)) - 1).real)
	lovek2Im = float((-1.0*pow(10.0,power-3)*y[3*layerNum*M-1]/(R*gravCompute(R,param)) - 1).imag)
	lovel2Re = float((pow(10.0,power-3)*y[(3*layerNum-1)*M-1]/(numpy.sqrt(6.0)*R)).real)
	lovel2Im = float((pow(10.0,power-3)*y[(3*layerNum-1)*M-1]/(numpy.sqrt(6.0)*R)).imag)

	return[loveh2Re,loveh2Im,lovek2Re,lovek2Im,lovel2Re,lovel2Im]


def AGen(r,param,omega):
	"""
	Computes the A matrix for the point r, using satellite properties specified
	by param and forcing frequency omega.

	Input:
	r: (float) Radial coordinate at which the A matrix will be evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	A: (3x3 float) A matrix defined at the point r.
	"""

	#Compute viscoelastic Lame parameters for the point r.
	lameMu = lameMuTwiddle(r,param,omega)
	lameLambda = lameLambdaTwiddle(r,param,omega)

	#Linear combination of Lame parameters.
	lameBeta = lameLambda + 2.0*lameMu

	A = numpy.array([[lameBeta,    0.0, 0.0],
					 [     0.0, lameMu, 0.0],
					 [     0.0,    0.0, 1.0]]) 
	
	return A


def BGen(r,param,omega):
	"""
	Computes the B matrix for the point r, using satellite properties specified
	by param and forcing frequency omega.

	Input:
	r: (float) Radial coordinate at which the B matrix will be evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	B: (3x3 float) B matrix defined at the point r.
	"""

	#Compute body parameters for the point r.
	lameMu = lameMuTwiddle(r,param,omega)
	lameLambda = lameLambdaTwiddle(r,param,omega)
	rho = rhoCompute(r,param)
	G = 6.67384e-11

	#Linear combination of Lame parameters.
	lameBeta = lameLambda + 2.0*lameMu

	B = numpy.array([[                       2.0*lameBeta/r, -numpy.sqrt(6.0)*(lameLambda+lameMu)/r,  -rho],
					 [numpy.sqrt(6.0)*(lameLambda+lameMu)/r,                           2.0*lameMu/r,   0.0],
					 [                   4.0*numpy.pi*G*rho,                                    0.0, 2.0/r]])
	
	return B


def CGen(r,param,omega):
	"""
	Computes the C matrix for the point r, using satellite properties specified
	by param and forcing frequency omega.

	Input:
	r: (float) Radial coordinate at which the B matrix will be evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	C: (3x3 float) C matrix defined at the point r.
	"""

	#Compute body parameters for the point r.
	lameMu = lameMuTwiddle(r,param,omega)
	lameLambda = lameLambdaTwiddle(r,param,omega)
	rho = rhoCompute(r,param)
	G = 6.67384e-11
	g = gravCompute(r,param)

	#Linear combination of Lame parameters.
	lameBeta = lameLambda + 2.0*lameMu

	C = numpy.array([[4.0*rho*g/r - 4.0*numpy.pi*G*pow(rho,2) - (2.0*lameBeta+6.0*lameMu)/pow(r,2), numpy.sqrt(6.0)*(lameLambda+3.0*lameMu)/pow(r,2) - numpy.sqrt(6)*rho*g/r,                    0.0],
					 [                 numpy.sqrt(6)*2.0*lameBeta/pow(r,2) - numpy.sqrt(6)*rho*g/r,                                                     -6*lameBeta/pow(r,2), -numpy.sqrt(6.0)*rho/r],
					 [                                                        8.0*numpy.pi*G*rho/r,                                      -4.0*numpy.sqrt(6)*numpy.pi*G*rho/r,          -6.0/pow(r,2)]])
	
	return C


def DGen(r,param):
	"""
	Computes the D matrix for the point r, using satellite properties specified
	by param.

	Input:
	r: (float) Radial coordinate at which the D matrix will be evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	Output:
	D: (3x3 float) D matrix defined at the point r.
	"""

	#Compute body parameters for the point r.
	rho = rhoCompute(r,param)

	D = numpy.array([[rho, 0.0, 0.0],
					 [0.0, rho, 0.0],
					 [0.0, 0.0, 0.0]])

	return D


def EGen(r,param,omega):
	"""
	Computes the E matrix for the point r, using satellite properties specified
	by param and forcing frequency omega.

	Input:
	r: (float) Radial coordinate at which the B matrix will be evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	E: (3x3 float) E matrix defined at the point r.
	"""

	#Compute body parameters for the point r.
	lameLambda = lameLambdaTwiddle(r,param,omega)
	lameMu = lameMuTwiddle(r,param,omega)
	rho = rhoCompute(r,param)
	G = 6.67384e-11

	E = numpy.array([[        2.0*lameLambda/r, -numpy.sqrt(6.0)*lameLambda/r,   0.0],
					 [numpy.sqrt(6.0)*lameMu/r,                     -lameMu/r,   0.0],
					 [      4.0*numpy.pi*G*rho,                           0.0, 3.0/r]])

	return E


def FGen(M):
	"""
	Creates the F matrix containing satellite center boundary conditions.

	Input:
	M: (int) Number of node points in first layer.

	Output:
	F: (3x3M float) Matrix appropriate for satellite center boundary conditions.
	"""
	
	F = numpy.zeros((3,3*M))
	F[0,0] = 1.0
	F[1,M] = 1.0
	F[2,2*M] = 1.0

	return F


def QGen(r,layerNum,M,param):
	"""
	Computes the Q matrix for the point r, using satellite properties specified
	by param for a model with M nodes per layer.

	Input:
	r: (float) Radial coordinate at which the B matrix will be evaluated.

	layerNum: (int) Number of layers in the satellite.

	M: (int) Number of node points per layer. Total number of nodes is 4M.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	Output:
	Q: (12Mx1 float) Q matrix defined at the point r.
	"""

	#All rows but the last are zero.
	Q = numpy.zeros((3*layerNum*M,1),'complex')
	Q[-1,0] = -5.0*gravCompute(r,param)

	return Q


def interiorGen(rSet,M,param,omega):
	"""
	Creates the matrix elements corresponding to the interior of a non-homogeneous
	layer.

	Input:
	rSet: (1xM float) Radial coordinates of set of nodes in the layer (including)
		   boundaries.

	M: (int) Number of node points per layer.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	PInterior: (3(M-2)x3M float) Set of matrix elements to be appended to P matrix
			    that contain layer interior equations.
	"""

	#Define PInterior matrix as zeros which will be filled.
	PInterior = numpy.zeros((3*(M-2),3*M),'complex')

	#Cycle through interior layer nodes not on boundaries. 
	for index in range(1,M-1):

		#Determine 0th, 1st, and 2nd derivative quadrature coefficients for the
		#inside nodes.
		[alphaCoef, betaCoef, gammaCoef] = quadCoef(rSet[index],rSet,2)

		#Generate A, B, C, and D matrices for current node.
		ATemp = AGen(rSet[index],param,omega)
		ATemp = numpy.vstack(( numpy.hstack((gammaCoef*ATemp[0,0],gammaCoef*ATemp[0,1],gammaCoef*ATemp[0,2])),
							   numpy.hstack((gammaCoef*ATemp[1,0],gammaCoef*ATemp[1,1],gammaCoef*ATemp[1,2])),
							   numpy.hstack((gammaCoef*ATemp[2,0],gammaCoef*ATemp[2,1],gammaCoef*ATemp[2,2])) ))

		BTemp = BGen(rSet[index],param,omega)
		BTemp = numpy.vstack(( numpy.hstack((betaCoef*BTemp[0,0],betaCoef*BTemp[0,1],betaCoef*BTemp[0,2])),
		 					   numpy.hstack((betaCoef*BTemp[1,0],betaCoef*BTemp[1,1],betaCoef*BTemp[1,2])),
		 					   numpy.hstack((betaCoef*BTemp[2,0],betaCoef*BTemp[2,1],betaCoef*BTemp[2,2])) ))

		CTemp = CGen(rSet[index],param,omega)
		CTemp = numpy.vstack(( numpy.hstack((alphaCoef*CTemp[0,0],alphaCoef*CTemp[0,1],alphaCoef*CTemp[0,2])),
							   numpy.hstack((alphaCoef*CTemp[1,0],alphaCoef*CTemp[1,1],alphaCoef*CTemp[1,2])),
							   numpy.hstack((alphaCoef*CTemp[2,0],alphaCoef*CTemp[2,1],alphaCoef*CTemp[2,2])) ))  

		DTemp = DGen(rSet[index],param)
		DTemp = numpy.vstack(( numpy.hstack((alphaCoef*DTemp[0,0],alphaCoef*DTemp[0,1],alphaCoef*DTemp[0,2])),
							   numpy.hstack((alphaCoef*DTemp[1,0],alphaCoef*DTemp[1,1],alphaCoef*DTemp[1,2])),
							   numpy.hstack((alphaCoef*DTemp[2,0],alphaCoef*DTemp[2,1],alphaCoef*DTemp[2,2])) ))

		#Store information for current node in PInterior.
		PInterior[3*(index-1):3*(index-1)+3,:] = ATemp + BTemp + CTemp + pow(omega,2)*DTemp

	return PInterior


def boundSL(rSetS,rSetL,M,param,omega):
	"""
	Creates the matrix elements corresponding to the boundary between a solid 
	layer below a liquid layer. The liquid layer is assumed to have zero shear
	modulus exactly. The layers are prefixed for (s)olid and (l)iquid.

	Input:
	rSetS: (1xM float) Radial coordinates of set of nodes in the solid layer.

	rSetL: (1xM float) Radial coordinates of set of nodes in the liquid layer.

	M: (int) Number of node points per layer.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	PBound: (6x3M float) Set of matrix elements to be appended to P matrix that
		     contain boundary layer equations.
	""" 

	#Compute body parameters at the boundary.
	G = 6.67384e-11
	g = gravCompute(rSetS[-1],param)

	lameLambdaS = lameLambdaTwiddle(rSetS[-1],param,omega)
	lameMuS = lameMuTwiddle(rSetS[-1],param,omega)
	lameBetaS = lameLambdaS + 2.0*lameMuS
	rhoS = rhoCompute(rSetS[-1],param)

	[alphaCoefS, betaCoefS] = quadCoef(rSetS[-1],rSetS,1)

	#Node 1 is used because the parameter functions return parameters for the 
	#lower layer when at the boundary exactly.
	lameLambdaL = lameLambdaTwiddle(rSetL[1],param,omega)
	lameMuL = lameMuTwiddle(rSetL[1],param,omega)
	lameBetaL = lameLambdaL + 2.0*lameMuL
	rhoL = rhoCompute(rSetL[1],param)

	[alphaCoefL, betaCoefL] = quadCoef(rSetL[0],rSetL,1)

	#Add each boundary condition to the PBound matrix one by one.

	#1. Continuity of radial displacement.
	zer = numpy.zeros((1,6*M))
	zer[0,M-1] = 1.0
	zer[0,3*M] = -1.0

	PBound = zer

	#2. Continuity of radial traction - note zero multiplication simply gives 
	#zero submatrices of the correct size.
	sTemp1 = numpy.hstack((betaCoefS*lameBetaS, betaCoefS*0.0, betaCoefS*0.0))
	sTemp2 = numpy.hstack((alphaCoefS*2.0*lameLambdaS/rSetS[-1], -alphaCoefS*numpy.sqrt(6)*lameLambdaS/rSetS[-1], alphaCoefS*0.0))

	lTemp1 = numpy.hstack((-betaCoefL*lameBetaL, betaCoefL*0.0, betaCoefL*0.0))
	lTemp2 = numpy.hstack((-alphaCoefL*2.0*lameLambdaL/rSetL[0], alphaCoefL*numpy.sqrt(6)*lameLambdaL/rSetL[0], alphaCoefL*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((sTemp1+sTemp2,lTemp1+lTemp2)) ))

	#3. Continuity of shear traction - note that zero multiplication simply 
	#gives zero submatrices of the correct size.
	zer = numpy.zeros((1,3*M))

	sTemp1 = numpy.hstack((betaCoefS*0.0, betaCoefS*lameMuS, betaCoefS*0.0))
	sTemp2 = numpy.hstack((alphaCoefS*numpy.sqrt(6.0)*lameMuS/rSetS[-1], -alphaCoefS*lameMuS/rSetS[-1], alphaCoefS*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((sTemp1+sTemp2, zer[0,:])) ))

	#4. Continuity of gravitational potential.
	zer = numpy.zeros((1,6*M))
	zer[0,3*M-1] = 1
	zer[0,5*M] = -1

	PBound = numpy.vstack((PBound, zer))

	#5. Continuity of normal component of gradient of gravitational potential.
	sTemp1 = numpy.hstack((betaCoefS*0.0, betaCoefS*0.0, betaCoefS))
	sTemp2 = numpy.hstack((alphaCoefS*4.0*numpy.pi*G*rhoS, alphaCoefS*0.0, alphaCoefS*0.0))

	lTemp1 = numpy.hstack((betaCoefL*0.0, betaCoefL*0.0, -betaCoefL))
	lTemp2 = numpy.hstack((-alphaCoefL*4.0*numpy.pi*G*rhoL, alphaCoefL*0.0, alphaCoefL*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((sTemp1+sTemp2,lTemp1+lTemp2)) ))

	#6. Equation for tangential displacement applied to liquid boundary.
	zer = numpy.zeros((1,3*M))

	lTemp1 = numpy.hstack((betaCoefL*numpy.sqrt(6.0)*lameLambdaL/rSetL[0], betaCoefL*0.0, betaCoefL*0.0))
	lTemp2 = numpy.hstack((alphaCoefL*(2.0*numpy.sqrt(6.0)*lameLambdaL/pow(rSetL[0],2) - numpy.sqrt(6.0)*rhoL*g/rSetL[0]), alphaCoefL*( -6.0*lameLambdaL/pow(rSetL[0],2) + rhoL*pow(omega,2) ), -alphaCoefL*numpy.sqrt(6)*rhoL/rSetL[0] ))

	PBound = numpy.vstack((PBound, numpy.hstack((zer[0,:],lTemp1+lTemp2)) ))
	#PBound is now 6x6M.

	return PBound


def boundLS(rSetL,rSetS,M,param,omega):
	"""
	Creates the matrix elements corresponding to the boundary between a liquid 
	layer below a solid layer. The liquid layer is assumed to have zero shear
	modulus exactly. The layers are prefixed for (l)iquid and (s)olid.

	Input:
	rSetL: (1xM float) Radial coordinates of set of nodes in the liquid layer.

	rSetS: (1xM float) Radial coordinates of set of nodes in the solid layer.

	M: (int) Number of node points per layer.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	PBound: (6x3M float) Set of matrix elements to be appended to P matrix that
		     contain boundary layer equations.
	""" 

	#Compute body parameters at the boundary.
	G = 6.67384e-11
	g = gravCompute(rSetL[-1],param)

	lameLambdaL = lameLambdaTwiddle(rSetL[-1],param,omega)
	lameMuL = lameMuTwiddle(rSetL[-1],param,omega)
	lameBetaL = lameLambdaL + 2.0*lameMuL
	rhoL = rhoCompute(rSetL[-1],param)

	[alphaCoefL, betaCoefL] = quadCoef(rSetL[-1],rSetL,1)

	#Node 1 is used because the parameter functions return parameters for the 
	#lower layer when at the boundary exactly.
	lameLambdaS = lameLambdaTwiddle(rSetS[1],param,omega)
	lameMuS = lameMuTwiddle(rSetS[1],param,omega)
	lameBetaS = lameLambdaS + 2.0*lameMuS
	rhoS = rhoCompute(rSetS[1],param)

	[alphaCoefS, betaCoefS] = quadCoef(rSetS[0],rSetS,1)

	#Add each boundary condition to the PBound matrix one by one.

	#1. Continuity of radial displacement.
	zer = numpy.zeros((1,6*M))
	zer[0,M-1] = 1.0
	zer[0,3*M] = -1.0

	PBound = zer

	#2. Continuity of radial traction - note zero multiplication simply gives 
	#zero submatrices of the correct size.
	lTemp1 = numpy.hstack((betaCoefL*lameBetaL, betaCoefL*0.0, betaCoefL*0.0))
	lTemp2 = numpy.hstack((alphaCoefL*2.0*lameLambdaL/rSetL[-1], -alphaCoefL*numpy.sqrt(6)*lameLambdaL/rSetL[-1], alphaCoefL*0.0))

	sTemp1 = numpy.hstack((-betaCoefS*lameBetaS, betaCoefS*0.0, betaCoefS*0.0))
	sTemp2 = numpy.hstack((-alphaCoefS*2.0*lameLambdaS/rSetS[0], alphaCoefS*numpy.sqrt(6)*lameLambdaS/rSetS[0], alphaCoefS*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((lTemp1+lTemp2,sTemp1+sTemp2)) ))

	#3. Continuity of shear traction - note that zero multiplication simply 
	#gives zero submatrices of the correct size.
	zer = numpy.zeros((1,3*M))

	sTemp1 = numpy.hstack((betaCoefS*0.0, betaCoefS*lameMuS, betaCoefS*0.0))
	sTemp2 = numpy.hstack((alphaCoefS*numpy.sqrt(6.0)*lameMuS/rSetS[0], -alphaCoefS*lameMuS/rSetS[0], alphaCoefS*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((zer[0,:], sTemp1+sTemp2)) ))

	#4. Continuity of gravitational potential.
	zer = numpy.zeros((1,6*M))
	zer[0,3*M-1] = 1
	zer[0,5*M] = -1

	PBound = numpy.vstack((PBound, zer))

	#5. Continuity of normal component of gradient of gravitational potential.
	lTemp1 = numpy.hstack((betaCoefL*0.0, betaCoefL*0.0, betaCoefL))
	lTemp2 = numpy.hstack((alphaCoefL*4.0*numpy.pi*G*rhoL, alphaCoefL*0.0, alphaCoefL*0.0))

	sTemp1 = numpy.hstack((betaCoefS*0.0, betaCoefS*0.0, -betaCoefS))
	sTemp2 = numpy.hstack((-alphaCoefS*4.0*numpy.pi*G*rhoS, alphaCoefS*0.0, alphaCoefS*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((lTemp1+lTemp2, sTemp1+sTemp2)) ))

	#6. Equation for tangential displacement applied to liquid boundary.
	zer = numpy.zeros((1,3*M))

	lTemp1 = numpy.hstack((betaCoefL*numpy.sqrt(6.0)*lameLambdaL/rSetL[-1], betaCoefL*0.0, betaCoefL*0.0))
	lTemp2 = numpy.hstack((alphaCoefL*(2.0*numpy.sqrt(6.0)*lameLambdaL/pow(rSetL[-1],2) - numpy.sqrt(6.0)*rhoL*g/rSetL[-1]), alphaCoefL*( -6.0*lameLambdaL/pow(rSetL[-1],2) + rhoL*pow(omega,2) ), -alphaCoefL*numpy.sqrt(6)*rhoL/rSetL[-1] ))

	PBound = numpy.vstack((PBound, numpy.hstack((lTemp1+lTemp2, zer[0,:])) ))
	#PBound is now 6x6M.

	return PBound


def boundSS(rSetL,rSetU,M,param,omega):
	"""
	Creates the matrix elements corresponding to the boundary between two solid 
	layers. The layers are prefixed for (l)ower and (u)pper.

	Input:
	rSetL: (1xM float) Radial coordinates of set of nodes in the lower layer.

	rSetU: (1xM float) Radial coordinates of set of nodes in the upper layer.

	M: (int) Number of node points per layer.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	PBound: (6x3M float) Set of matrix elements to be appended to P matrix that
		     contain boundary layer equations.
	""" 

	#Compute body parameters at the boundary.
	G = 6.67384e-11
	g = gravCompute(rSetL[-1],param)

	lameLambdaL = lameLambdaTwiddle(rSetL[-1],param,omega)
	lameMuL = lameMuTwiddle(rSetL[-1],param,omega)
	lameBetaL = lameLambdaL + 2.0*lameMuL
	rhoL = rhoCompute(rSetL[-1],param)

	[alphaCoefL, betaCoefL] = quadCoef(rSetL[-1],rSetL,1)

	#Node 1 is used because the parameter functions return parameters for the 
	#lower layer when at the boundary exactly.
	lameLambdaU = lameLambdaTwiddle(rSetU[1],param,omega)
	lameMuU = lameMuTwiddle(rSetU[1],param,omega)
	lameBetaU = lameLambdaU + 2.0*lameMuU
	rhoU = rhoCompute(rSetU[1],param)

	[alphaCoefU, betaCoefU] = quadCoef(rSetU[0],rSetU,1)

	#Add each boundary condition to the PBound matrix one by one.

	#1. Continuity of radial displacement.
	zer = numpy.zeros((1,6*M))
	zer[0,M-1] = 1.0
	zer[0,3*M] = -1.0

	PBound = zer

	#2. Continuity of radial traction - note zero multiplication simply gives 
	#zero submatrices of the correct size.
	lTemp1 = numpy.hstack((betaCoefL*lameBetaL, betaCoefL*0.0, betaCoefL*0.0))
	lTemp2 = numpy.hstack((alphaCoefL*2.0*lameLambdaL/rSetL[-1], -alphaCoefL*numpy.sqrt(6)*lameLambdaL/rSetL[-1], alphaCoefL*0.0))

	uTemp1 = numpy.hstack((-betaCoefU*lameBetaU, betaCoefU*0.0, betaCoefU*0.0))
	uTemp2 = numpy.hstack((-alphaCoefU*2.0*lameLambdaU/rSetU[0], alphaCoefU*numpy.sqrt(6)*lameLambdaU/rSetU[0], alphaCoefU*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((lTemp1+lTemp2, uTemp1+uTemp2)) ))

	#3. Continuity of shear traction - note that zero multiplication simply 
	#gives zero submatrices of the correct size.
	lTemp1 = numpy.hstack((betaCoefL*0.0, betaCoefL*lameMuL, betaCoefL*0.0))
	lTemp2 = numpy.hstack((alphaCoefL*numpy.sqrt(6.0)*lameMuL/rSetL[-1], -alphaCoefL*lameMuL/rSetL[-1], alphaCoefL*0.0))

	uTemp1 = numpy.hstack((betaCoefU*0.0, -betaCoefU*lameMuU, betaCoefU*0.0))
	uTemp2 = numpy.hstack((-alphaCoefU*numpy.sqrt(6.0)*lameMuU/rSetU[0], alphaCoefU*lameMuU/rSetU[0], alphaCoefU*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((lTemp1+lTemp2, uTemp1+uTemp2)) ))

	#4. Continuity of gravitational potential.
	zer = numpy.zeros((1,6*M))
	zer[0,3*M-1] = 1
	zer[0,5*M] = -1

	PBound = numpy.vstack((PBound, zer))

	#5. Continuity of normal component of gradient of gravitational potential.
	lTemp1 = numpy.hstack((betaCoefL*0.0, betaCoefL*0.0, betaCoefL))
	lTemp2 = numpy.hstack((alphaCoefL*4.0*numpy.pi*G*rhoL, alphaCoefL*0.0, alphaCoefL*0.0))

	uTemp1 = numpy.hstack((betaCoefU*0.0, betaCoefU*0.0, -betaCoefU))
	uTemp2 = numpy.hstack((-alphaCoefU*4.0*numpy.pi*G*rhoU, alphaCoefU*0.0, alphaCoefU*0.0))

	PBound = numpy.vstack((PBound, numpy.hstack((lTemp1+lTemp2, uTemp1+uTemp2)) ))

	#6. Continuity of tangential displacement.
	zer = numpy.zeros((1,6*M))
	zer[0,2*M-1] = 1
	zer[0,4*M] = -1
	
	PBound = numpy.vstack((PBound, zer))
	#PBound is now 6x6M.

	return PBound


def boundFree(rSet,param,omega):
	"""
	Creates the matrix elements corresponding to the free boundary of the body. 

	Input:
	rSet: (1xM float) Radial coordinates of set of nodes in the upmost layer.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	PBound: (3x3M float) Set of matrix elements to be appended to P matrix that
		     contain free boundary layer equations.
	""" 

	#Compute differential quadrature coefficients at free surface.
	[alphaCoefFree, betaCoefFree] = quadCoef(rSet[-1],rSet,1)

	#Generate A and E matrices for free surface node.
	ATemp = AGen(rSet[-1],param,omega)
	ATemp = numpy.vstack(( numpy.hstack((betaCoefFree*ATemp[0,0],betaCoefFree*ATemp[0,1],betaCoefFree*ATemp[0,2])),
						   numpy.hstack((betaCoefFree*ATemp[1,0],betaCoefFree*ATemp[1,1],betaCoefFree*ATemp[1,2])),
						   numpy.hstack((betaCoefFree*ATemp[2,0],betaCoefFree*ATemp[2,1],betaCoefFree*ATemp[2,2])) ))

	ETemp = EGen(rSet[-1],param,omega)
	ETemp = numpy.vstack(( numpy.hstack((alphaCoefFree*ETemp[0,0],alphaCoefFree*ETemp[0,1],alphaCoefFree*ETemp[0,2])),
						   numpy.hstack((alphaCoefFree*ETemp[1,0],alphaCoefFree*ETemp[1,1],alphaCoefFree*ETemp[1,2])),
						   numpy.hstack((alphaCoefFree*ETemp[2,0],alphaCoefFree*ETemp[2,1],alphaCoefFree*ETemp[2,2])) ))

	PBound = ATemp + ETemp

	return PBound


def lameMuTwiddle(r,param,omega):
	"""
	Returns the complex shear modulus of the point r. The core is currently
	treated as an elastic body, as viscoelastic effects are small.

	Input:
	r: (float) Radial coordinate at which the complex shear modulus will be 
	evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	lameMu: (float) Complex shear modulus of the point r.
	"""

	#While loop determines the layer of point r.
	layer = 0

	while r > param[0,layer]*1000.0:

		layer = layer+1

	lameMu = param[3,layer]

	return lameMu


def lameLambdaTwiddle(r,param,omega):
	"""
	Returns the complex Lame lambda parameter of the point r. The core is 
	currently treated as an elastic body, as viscoelastic effects are small.

	Input:
	r: (float) Radial coordinate at which the complex Lame lambda parameter will 
	be evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	lameLambda: (float) Complex Lame lambda parameter of the point r.
	"""

	#While loop determines the layer of point r.
	layer = 0

	while r > param[0,layer]*1000.0:

		layer = layer+1

	lameLambda = param[4,layer]

	return lameLambda


def rhoCompute(r,param):
	"""
	Returns the density of the point r.

	Input:
	r: (float) Radial coordinate at which the denisity will be evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	rho: (float) Density of the point r.
	"""

	#While loop determines the layer of point r.
	layer = 0

	while r > param[0,layer]*1000.0:

		layer = layer+1

	rho = param[2,layer]

	return rho


def gravCompute(r,param):
	"""
	Computes the gravitational acceleration g at the point r, using satellite 
	properties specified by param.

	Input:
	r: (float) Radial coordinate at which the gravitational acceleration will be 
	evaluated.

	param: (5x4 float) Matrix containing parameters for layers. Column 0 is the 
		    innermost layer, and the rightmost column is the outer layer.

			[Layer outer radius (km)
			 State of layer (0 for liquid, 1 for solid)
			 Density (kg/m^3)
			 Viscoelastic shear modulus (Pa)
			 Viscoelastic Lame's first parameter (Pa)]

	omega: (float) Forcing frequency (rad/s).

	Output:
	g: (float) Gravitational acceleration at the point r.
	"""

	#Physical constants.
	G = 6.67384e-11

	#While loop determines the layer of point r.
	layer = 0

	while r > param[0,layer]*1000.0:

		layer = layer+1

	#The following uses a simple formula the avoids integration, valid for
	#uniform layers.
	if layer == 0:

		g = 4.0/3.0*numpy.pi*rhoCompute(r,param)*G*r

	else:

		#Innermost and outermost layers first.
		gNumerator = rhoCompute(param[0,0]*1000.0,param)*pow(param[0,0]*1000.0,3) + rhoCompute(r,param)*( pow(r,3) - pow(param[0,layer-1]*1000.0,3) )

		#Then contributions of all intermediate layers.
		for layerIndex in range(1,layer):

			gNumerator = gNumerator + rhoCompute(param[0,layerIndex]*1000.0,param)*( pow(param[0,layerIndex]*1000.0,3) - pow(param[0,layerIndex-1]*1000.0,3) )

		g = 4.0/3.0*numpy.pi*G*gNumerator/pow(r,2)

	return g


def quadCoef(r,rSet,order):
	"""
	Finds the differential quadrature coefficients needed to compute derivatives
	up to order order, at point r, using points specified by rSet. Uses an 
	algorithm developed by Fornberg, "Practical Guide to Pseudospectral Methods".
	Alpha and beta are intermediate algorithm values unrelated to those used in
	stress computation. 

	Input:
	r: (float) Point at which derivatives which will be computed.

	rSet: (Mx1 float) Set of points in the layer to use in quadrature.

	order: (int) Order of highest derivative needed. 

	Output:
	coefOut: (orderx1 list of Mx1 floats) Differential quadrature coefficients to
				to be multiplied by each element of rSet to approximate derivatives
				at r. coefOut[0] corresponds to 0th derivative, coefOut[1] is first,
				etc. 
	"""

	M = rSet.size

	#Coefficient matrix to store intermediate values.
	coefMat = numpy.zeros((M,M,order+1))
	coefMat[0,0,0] = 1;

	#Begin algorithm.
	alpha = 1;

	for i in range(1,M):

		beta = 1

		for j in range(0,i):

			beta = beta*(rSet[i] - rSet[j])

			for k in range(0,min(i,order)+1):

				if k == 0:

					coefMat[i,j,k] = ( (rSet[i]-r)*coefMat[i-1,j,k] )/(rSet[i] - rSet[j])

				else:

					coefMat[i,j,k] = ( (rSet[i]-r)*coefMat[i-1,j,k] - k*coefMat[i-1,j,k-1] )/(rSet[i] - rSet[j])

		for k in range(0,min(i,order)+1):

			if k == 0:

				coefMat[i,i,k] = alpha*( -(rSet[i-1] - r)*coefMat[i-1,i-1,k] )/beta

			else:

				coefMat[i,i,k] = alpha*( k*coefMat[i-1,i-1,k-1] - (rSet[i-1] - r)*coefMat[i-1,i-1,k] )/beta

		alpha = beta

	#Algorithm complete, but the final values are stored in coefMat[-1,:,:]. Make these 
	#as seperate elements in a list.
	coefOut = []

	for orderIndex in range(0,order+1):

		coefOut.append(coefMat[-1,:,orderIndex])

	return coefOut
