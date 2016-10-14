from SALib.sample import saltelli
from SALib.analyze import sobol
import api as api
import numpy as np
from joblib import Parallel, delayed
from os import remove
import uuid
import div_tools as dt

def v2p(param):
	p = param.reshape((22,3))
	p = np.r_[p[0:9,:],[np.zeros(3)],p[9:15,:],[np.zeros(3)],p[15:,:]]
	return p

def f(p,pos):
	'''
		an objective function
	'''
	fname = 'temp_' + str(uuid.uuid4())
	api.write_param(param = v2p(p),filename_save = fname)
	output_ = api.output(param=fname)
	remove(fname)
	
	return(pos,output_[3])

problem = api.generate_problem(delta=0.05)
# Generate samples
param_values = saltelli.sample(problem, 10, calc_second_order=False)
nparam=param_values.shape[0]
# Run model
#Y = np.zeros([param_values.shape[0]])

n_par = 3
if n_par>1: 
	results = Parallel(n_jobs=n_par,verbose=10)(delayed(f)(param_values[i,:],i) for i in xrange(nparam))
else:
	results=np.zeros((nparam,2))
	for i in xrange(nparam):
		results[i] = f(param_values[i,:],i)
		if i%10==0:
			np.savetxt('save_mc.runs',results)
			print('runs: ' + str(i) + ' / ' + str(nparam))

print 'comp are finished'
results.sort()
dt.save(obj = results, filename = 'data_sobol.dat')
results2 = np.array([r[1] for r in results])
#for i, X in enumerate(param_values):
#	if i%25==0:
#		np.savetxt('save_mc.runs',Y)
#	if i%100==0:
#		print('runs: ' + str(i) + ' / ' + str(param_values.shape[0]))
#	Y[i] = f(X)

# Perform analysis
Si = sobol.analyze(problem, results2, print_to_console=False,calc_second_order = False)

# Print the first-order sensitivity indexes
print Si['S1']
# Print the total-order sensitivity indexes
print Si['ST']
