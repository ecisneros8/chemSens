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
		r_param[9,:] = 0.
		r_param[16,:] = 0.
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


def old_change_parameters(mech = 'sandiego', param = False):
	if param is False:
		print('parameters needed')
		return
	fi = fileinput.FileInput("rhsjac.cpp", inplace=1)
	if mech == 'sandiego':
		indmin = 502
		for line in fi:
			if fi.lineno()>indmin-1:
				n = fi.lineno() - indmin
				if n>8: n+=1 # no kfwd9?
				if n>15: n+=1 # no kfwd16?
				if n < param.shape[0]:
					line = '    kfwd[' + str(n) + '] =  exp(' 
					if param[n,0]>0: line += ' + '
					line +=  "{:.8e}".format(param[n,0]) + ' ' 
					if param[n,1]>0: line += ' + '
					line+= "{:.6e}".format(param[n,1]) + ' * tlog ' 
					if param[n,2]>0: line += ' + '
					line+= "{:.6e}".format(param[n,2]) + ' * rt);\n'
			print line,

	else: 
		indmin = 473
		for line in fi:
			if fi.lineno()>indmin-1:
				n = fi.lineno() - indmin
				if n>2: n+=1 # no kfwd3?
				if n < param.shape[0]:
					line = str(n) + '    kfwd[' + str(n) + '] =  exp(' 
					if param[n,0]>0: line += ' + '
					line +=  "{:.8e}".format(param[n,0]) + ' ' 
					if param[n,1]>0: line += ' + '
					line+= "{:.6e}".format(param[n,1]) + ' * tlog ' 
					if param[n,2]>0: line += ' + '
					line+= "{:.6e}".format(param[n,2]) + ' * rt);\n'
			print line,

def old_cp(mech = 'sandiego', param = False):
	old_change_parameters(mech = mech, param = param)

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
	if mech == 'sandiego':param = np.delete(param,[9,16],axis=0)
	else:param = np.delete(param,[3],axis=0)
	bounds = np.zeros((param.size,2))
	ecart = np.tile(delta*np.abs(np.max(param,axis=0) - np.min(param,axis=0)),param.shape[0])
	bounds[:,0] = param.flatten() - ecart
	bounds[:,1] = param.flatten() + ecart
	problem = {'num_vars': param.size,
		'bounds': bounds}
	return problem
