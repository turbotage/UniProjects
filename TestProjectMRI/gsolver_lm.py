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
	
	def solve(self, threshold = 0.0001, max_iter=1000):
		
		
		
		guess = self.__init_guess
		data = self.__init_data
		dependent = self.__init_dependent
		
		func = self.__get_func(dependent)
		f = lambda p: func(p).sum(dim=1)
		
		non_convergence_list = torch.arange(0,self.__npixels).cuda()
		
		
		
		
		
			"""
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
	
			"""
	
	def dogleg(res, J, JT, delta):
		p = torch.zeros(self.__npixels, self.__nparams)
		step = torch.zeros(self.__npixels, self.__nparams)
		
		Jn2 = (J * J).sum(dim=1)
		Jn = torch.sqrt(Jn2)
		
		D = torch.diag_embed(1 / Jn)
		
		Js = torch.bmm(J,D)
		JsT = Js.transpose(1,2)
		
		Hs = torch.bmm(JsT,Js)
		gs = torch.bmm(JsT,res)
		
		u = torch.cholesky(Hs)
		q = torch.cholesky_solve(gs,u)
		
		pGN = torch.bmm(D,q)
		
		case1 = torch.norm(pGN,dim=2) <= delta
		p[case1] = 
		
		
		
		
		
		
	def levenberg_marquardt_powell(self, delta0, mu, eta, threshold, max_iter):
		
		
		
		
		while True:
			
		
	
	
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


