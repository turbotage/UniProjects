# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:34:50 2021

@author: TurboTage
"""

import torch


class Model:
	def __init__(self, model_expr, dependent, data, init_guess):
		self.__model_expr = model_expr
		self.__model = self.__get_model(model_expr)
		
		self.__nparams = init_guess.size()[0]
		self.__npixels = init_guess.size()[1]
		
		# Initial values are saved
		self.__init_dependent = dependent
		self.__init_data = data
		self.__init_guess = init_guess
		
		
	def __get_model(self, model_expr):
		# Replace math expr in order of pytorch docs
		model_expr = model_expr.replace('abs', 'torch.abs')
		model_expr = model_expr.replace('acos', 'torch.acos')
		model_expr = model_expr.replace('acosh', 'torch.acosh')
		model_expr = model_expr.replace('asin', 'torch.sin')
		model_expr = model_expr.replace('asinh', 'torch.asinh')
		model_expr = model_expr.replace('atan', 'torch.atan')
		model_expr = model_expr.replace('atanh', 'torch.atanh')
		model_expr = model_expr.replace('atan2', 'torch.atan2')
		#model_expr = model_expr.replace('conj', 'torch.conj')
		model_expr = model_expr.replace('cos', 'torch.cos')
		model_expr = model_expr.replace('cosh', 'torch.cosh')
		model_expr = model_expr.replace('exp', 'torch.exp')
		model_expr = model_expr.replace('log', 'torch.log')
		model_expr = model_expr.replace('pow', 'torch.pow')
		model_expr = model_expr.replace('rsqrt', 'torch.rsqrt')
		model_expr = model_expr.replace('sin', 'torch.sin')
		model_expr = model_expr.replace('sinh', 'torch.sinh')
		model_expr = model_expr.replace('sqrt', 'torch.sqrt')
		model_expr = model_expr.replace('tan', 'torch.tan')
		model_expr = model_expr.replace('tanh', 'torch.tanh')
		
		# Paramaters
		for i in range(50):
			model_expr = model_expr.replace('P' + str(i), 'params[:][' + str(i) + ']')
		
		# Variables
		for i in range(50):
			model_expr = model_expr.replace('X' + str(i), 'dependent[' + str(i) + ']')
		
		#print(model_expr)        
		return lambda dependent, params: eval(model_expr)
		
	def __get_func(self, dependent):
		return lambda params: self.__model(dependent, params)
	
	def __get_jacobians(self, f, params):
		JT = torch.autograd.functional.jacobian(f, params).transpose(0,2)
		J = JT.transpose(1,2)
		return J, JT
  
	def __get_diff(self, func, params, data):
		diff = (data - func(params)).transpose(0,1).reshape(data.size()[1], 1, data.size()[0]).transpose(1,2)
		return diff
  
	def __get_step(self, JTJ, JT_diff, damping):
		# AX=B
		
		eye = torch.eye(JTJ.size()[1]).repeat(JTJ.size()[0],1).reshape(JTJ.size()[0],JTJ.size()[1],JTJ.size()[1]).cuda()
		
		perturb = JTJ * eye * damping.reshape(JTJ.size()[0],1,1)
		
		JTJ = JTJ + perturb
		
		# Handle singularity, inf
		good_behaved = (torch.logical_and(torch.det(JTJ) > 0.001, torch.det(JTJ) < 1e25))
		
		param_step = torch.zeros(JTJ.size()[1], JTJ.size()[0]).cuda()
		
		JTJ_n = JTJ[good_behaved]
		JT_diff_n = JT_diff[good_behaved]
		#print(JTJ_n.shape)
		#print(JT_diff_n.shape)
		
		# Solve via Cholskey decomposition, should be better than LU for this problem
		u = torch.cholesky(JTJ_n)
		#print(torch.cholesky_solve(JT_diff_n,u).transpose(0,2)[0].shape)
		#print(param_step.shape)
		#print(good_behaved.shape)
		
		param_step[:,good_behaved] = torch.cholesky_solve(JT_diff_n,u).transpose(0,2)[0]
		
		return param_step, torch.logical_not(good_behaved)
	
	def solve(self, threshold = 0.0001, max_iter=10000):
		uphil_step_mul = 1.2
		downhill_step_mul = 1.5
		
		guess = self.__init_guess
		data = self.__init_data
		dependent = self.__init_dependent
		
		func = self.__get_func(dependent)
		f = lambda p: func(p).sum(dim=1)
		
		
		non_convergence_list = torch.arange(0,self.__npixels).cuda()
		
		damping = torch.ones(1,self.__npixels).cuda() * 1e6
		old_damping = damping
		#mul = torch.ones(1,self.__npixels).cuda() * 1
		
		has_stepped = (torch.ones(1,self.__npixels) < 0).cuda()
		
		# Initial Step
		J, JT = self.__get_jacobians(f, guess)
		diff = self.__get_diff(func, guess, data)
		
		JT_diff = torch.bmm(JT, diff)
		JTJ = torch.bmm(JT,J)
		
		
		old_S = torch.zeros(1,self.__npixels).cuda()
		S = (diff * diff).sum(dim=1).transpose(0,1)
		
		#step0, ill_behaved0 = self.__get_step(JTJ, JT_diff, damping * 0)
		step1, ill_behaved1 = self.__get_step(JTJ, JT_diff, damping / downhill_step_mul)
		step2, ill_behaved2 = self.__get_step(JTJ, JT_diff, damping)
		
		diff = self.__get_diff(func, guess + step1, data)
		S1 = (diff * diff).sum(dim=1).transpose(0,1)
		
		diff = self.__get_diff(func, guess + step2, data)
		S2 = (diff * diff).sum(dim=1).transpose(0,1)
		
		damping = damping * ((S1 < S) / downhill_step_mul + torch.logical_and(S1 >= S, S2 <= S) + torch.logical_and(S1 >= S, S2 > S) * uphil_step_mul)
		
		old_step = torch.zeros(self.__nparams, self.__npixels).cuda()
		step = (S1 < S) * step1 + torch.logical_and(S1 >= S, S2 <= S) * step2
		
		converges = ((torch.abs((old_step - step) / guess) < (threshold * 10 / (1+damping))).sum(dim=0) + (torch.abs((old_S - S) / S) < (threshold / (1+damping))) + has_stepped)[0,:]
		
		params = guess
		#print("initial step done")
		
		iteration = 0
		while (converges.sum() != guess.size()[1] * (guess.size()[0] + 2)) and iteration < max_iter:
			iteration += 1
			
			has_stepped = damping <= old_damping
			
			guess = guess + step
			
			J, JT = self.__get_jacobians(f, guess)
			diff = self.__get_diff(func, guess, data)
			JT_diff = torch.bmm(JT, diff)
			JTJ = torch.bmm(JT,J)
			
			old_S = S
			S = (diff * diff).sum(dim=1).transpose(0,1)
			
			#step0, ill_behaved0 = self.__get_step(JTJ, JT_diff, damping * 0)
			step1, ill_behaved1 = self.__get_step(JTJ, JT_diff, damping / downhill_step_mul)
			step2, ill_behaved2 = self.__get_step(JTJ, JT_diff, damping)
			
			
			diff = self.__get_diff(func, guess + step1, data)
			S1 = (diff * diff).sum(dim=1).transpose(0,1)
			
			diff = self.__get_diff(func, guess + step2, data)
			S2 = (diff * diff).sum(dim=1).transpose(0,1)
			
			old_damping = damping
			damping = damping * ((S1 < S) / downhill_step_mul + torch.logical_and(S1 >= S, S2 <= S) + torch.logical_and(S1 >= S, S2 > S) * uphil_step_mul)
			
			old_step = step
			step = (S1 < S) * step1 + torch.logical_and(S1 >= S, S2 <= S) * step2
			
			#if ( (damping < 0.1).sum() > 0):
			#    print("Small damping")
			
			#print(damping.shape)
			#(old_damping.shape)
			converges = ((torch.abs((old_step - step) / guess) < (threshold * 10 / (1+damping))).sum(dim=0) + (torch.abs((old_S - S) / S) < (threshold / (1+damping))) + has_stepped)[0,:]
			
			has_stepped = damping <= old_damping
			
			convergence_percentage = float(guess.size()[1] - (converges < (self.__nparams + 2)).sum().cpu().numpy()) / float(guess.size()[1])
			
			if convergence_percentage * self.__npixels > 5000 or convergence_percentage > 0.50:
				not_converging = converges < (self.__nparams + 2)
				converging = torch.logical_not(not_converging)
				
				# copy converging pixels
				params[:,non_convergence_list[converging]] = guess[:,converging]
				
				
				# rebuild solver to solve with non-converging pixels
				
				dependent = dependent[:,:,not_converging]
				data = data[:,not_converging]
				guess = guess[:,not_converging]
				
				step1 = step1[:,not_converging]
				step2 = step2[:,not_converging]
				step = step[:,not_converging]
				
				S1 = S1[:,not_converging]
				S2 = S2[:,not_converging]
				S = S[:,not_converging]
				
				damping = damping[:,not_converging]
				old_damping = old_damping[:,not_converging]
				#mul = mul[:,not_converging]
				
				converges = converges[not_converging]
				
				func = self.__get_func(dependent)
				f = lambda p: func(p).sum(dim=1)
				
				# rebuild non convergence list
				non_convergence_list = non_convergence_list[not_converging]
				
				
				#print((1-float(self.__npixels - converges.sum().cpu().numpy()) / float(self.__npixels)) * float(self.__npixels))
				#print(iteration)
				#print((ill_behaved1 + ill_behaved2).sum() )
				
				#if guess.size()[1] < 20:
					#print(guess)
		
		
		params[:,non_convergence_list[converges > (self.__nparams + 1)]] = guess[:,converges > (self.__nparams + 1)]
		
		convergence_percentage = float(self.__npixels - (converges < (self.__nparams + 2)).sum().cpu().numpy()) / float(self.__npixels)
		
		return params, convergence_percentage, iteration

# END OF MODEL


# START OF SOLVER
   


# END OF SOLVER



# PARAMETER GUESSING

def get_uniform_guess(data,uniform_param_guess):
	guess = torch.zeros(len(uniform_param_guess),data.size()[1]).cuda()
	for i in range(len(uniform_param_guess)):
		guess[i,:] = uniform_param_guess[i]
		
	return guess



# END PARAMETER GUESSING


