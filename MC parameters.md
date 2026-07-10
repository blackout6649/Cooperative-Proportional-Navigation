## INPUT

tau\_m 	= U(0.005, 0.05) 		# mechanical lag time constant, in \[sec]
tau\_f 	= U(0.01, 0.1) 		# filter lag time constant, in \[sec]
tau\_sk 	= U(0.01, 0.2) 		# seeker lag time constant, in \[sec]

M\_x0\_vec 	= U^3i(-500, 500) 	# missile initial positions; mx1 cell; in \[m]
						% M\_x0\_vec{i} = \[x0, y0, z0] for missile i

M\_V\_vec\_3D 	= U^i(100, 150) 		# missile velocities; mx1 cell; in \[m/sec]
						% M\_V\_vec{i} = V\_tot for missile i

M\_gamma0\_vec 	= U^3i(60, 120) 	# missile initial headings; mx1 cell; in \[deg]
						% M\_gamma0\_vec{i} = \[phi0, theta0, psi0] for missile i

## ADDITIONAL

sigma\_gust 	= U^3(\[2, 2, 1], \[12, 12, 5]) 	# gust std per axis, in \[m/sec]

L\_gust 	= U^3(\[100, 100, 50], \[400, 400, 200]) 	# Gauss-Markov length scale per axis, in \[m]

v\_gust(t) 	= first-order Gauss-Markov process 	# 3 x n\_sim gust profile, in \[m/sec]
						% generated from sigma\_gust, L\_gust, V\_ref, dt

sigma\_RM 	= U^3(\[0, 0, 0], \[5, 5, 5]) 	# RM measurement-noise std per axis, in \[m]

sigma\_VM 	= U^3(\[0, 0, 0], \[2, 2, 2]) 	# VM measurement-noise std per axis, in \[m/sec]

## OUTPUT

N\_CPN(t) 		# time-varying navigation gain; mx1 cell of time-series

epsilon(t) 		# cooperative time-to-go error; mx1 cell of time-series

v\_gust(t) 		# gust time-series used in the run; 3 x n\_sim

R\_hit 			# hit radius used in the run, in \[m]

hit 			# hit flag per run

P\_hit 			# hit probability over all MC runs

t\_impact 		# impact time for hit runs, in \[sec]

hit\_angle 		# final hit attitude \[psi, theta, phi] from gammaM(:, end), in \[rad]

