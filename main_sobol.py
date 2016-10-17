from SALib.sample import saltelli
from SALib.analyze import sobol
import api as api
import numpy as np
from joblib import Parallel, delayed
from os import remove
import uuid
import div_tools as dt
from time import time

def v2p(param):
#	p = param.reshape((22,3))
#	p = np.r_[p[0:9,:],[np.zeros(3)],p[9:15,:],[np.zeros(3)],p[15:,:]]
	p = param.reshape((param.size/3,3))
	return p

def f(p,pos):
	'''
		Function for the sobol analysis
	'''
	fname = 'temp_' + str(uuid.uuid4())
	api.write_param(param = v2p(p),filename_save = fname)
	output_ = api.output(param=fname)
	remove(fname)
	
	return(pos,output_[3])

problem = api.generate_problem(delta=0.05)
# Generate samples
param_values = saltelli.sample(problem, 50, calc_second_order=False)
nparam=param_values.shape[0]

# Run model
n_par = 1
if n_par ==1: t0=time()
if n_par>1: 
	results = Parallel(n_jobs=n_par,verbose=10)(delayed(f)(param_values[i,:],i) for i in xrange(nparam))
else:
	results=np.zeros((nparam,2))
	for i in xrange(nparam):
		results[i] = f(param_values[i,:],i)
		if i%10==0:
			np.savetxt('save_mc.runs',results)
			print('runs: ' + str(i) + ' / ' + str(nparam))
			t = time() - t0
			time_left = (nparam - i) * t/(i+1)
			print('time so far: ' + str(int(t/60)) + 'mn')
			print('time before end is approx: ' + str(int(time_left/60)) + 'mn')

print 'comp are finished'

#save results
results.sort()
dt.save(obj = results, filename = 'data_sobol.dat')
results2 = np.array([r[1] for r in results])
#remove explosion if any
# !!!! to be improve !!
expl = np.where(results2 == -1)[0]
if expl.size>0:results2[expl] = (results2[expl-1] + results2[expl+1])/2
# Perform analysis
Si = sobol.analyze(problem, results2, print_to_console=False,calc_second_order = False)
ST = np.reshape(Si['ST'],(Si['ST'].size/3,3))
# Print the first-order sensitivity indexes
print Si['S1']
# Print the total-order sensitivity indexes
print Si['ST']
