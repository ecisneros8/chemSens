# python api for esteban code
import subprocess as sp
import numpy as np
import fileinput

def rparam(mech='sandiego'):
	if mech == 'sandiego':
		r_param = np.zeros((24,3))
	else:
		r_param = np.zeros((11,3))
	r_param[:,0] = np.random.uniform(1,20,r_param.shape[0])
	r_param[:,1] = np.random.uniform(-3,3,r_param.shape[0])
	r_param[:,2] = np.random.uniform(-100,2500,r_param.shape[0])
	if mech == 'sandiego':
		pass
	else:
		r_param[3,:] = 0.
	return r_param

def write_param(param = False, mech = 'sandiego', filename_load = False, filename_save = 'param.parameters'):
	if filename_load is not False:
		import shutil
		shutil.copy(filename,'param.parameters')
		return
	if param is False:
		param = rparam(mech = mech)
	np.savetxt(filename_save,param)

def wp(param = False, mech = 'sandiego', filename_load = False, filename_save = 'param.parameters'):
	write_param(param = param,mech = mech, filename_load = filename_load, filename_save = filename_save)

def recompile():
	sp.call(['make', 'clean'])
	sp.call(['make', 'exec'])

def change_mechanism(mech = 'sandiego'):
	if mech.lower() == 'boivin':
		sp.check_output(['./newsrc.sh', 'boivin'])
	elif mech.lower() == 'sandiego':
		sp.check_output(['./newsrc.sh', 'sanDiego'])
	else :
		print('SanDiego mechanism is selected')
		sp.check_output(['./newsrc.sh', 'sanDiego'])
	recompile()

def cm(mech = 'sandiego'):
	change_mechanism(mech = mech)

def run(temp = 1200, h = 1., param = 'param.parameters'):
	try:
		s = sp.check_output(['./exec', str(temp), str(h), param],)
	except sp.CalledProcessError as e:
#		print e.output
		s = 'Blow out \nneg time and temp are returned \nignition t:    -1.000E0 \nfinal temp:    -1.000E0'
	return s

def output(temp = 1200, h = 1., param = 'param.parameters'):
	s=run(temp = temp, h=h, param = param)
	ind_ = s.find('final temp:')
	ftemp = float(s[ind_+15:ind_+25])
	ind_ = s.find('ignition t:')
	it = float(s[ind_+15:ind_+25])
	return temp,ftemp,h,it


def generate_problem(delta= 0.05,mech='sandiego'):
	param = np.loadtxt(mech+'.parameters')
	if mech == 'boivin':param = np.delete(param,[3],axis=0)
	bounds = np.zeros((param.size,2))
	ecart = np.tile(delta*np.abs(np.max(param,axis=0) - np.min(param,axis=0)),param.shape[0])
	bounds[:,0] = param.flatten() - ecart
	bounds[:,1] = param.flatten() + ecart
	problem = {'num_vars': param.size,
		'bounds': bounds}
	return problem
