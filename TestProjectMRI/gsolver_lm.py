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
	
	def __get_jacobians(self, fun, params):
		JT = torch.autograd.functional.jacobian(fun, params).transpose(0,2)
		J = JT.transpose(1,2)
		return J, JT
  
	def __get_residuals(self, func, params, data):
		resT = (data - func(params)).transpose(0,1).reshape(data.size()[1], 1, data.size()[0])
		res = resT.transpose(1,2)
		return res, resT  
	
	def __converges(self, Jp, res, tol):
		return (torch.norm(Jp,dim=1) <= tol * torch.norm(res, dim=1))[:,0]
	        
	def __dogleg(self, res, J, delta):
		
		Jn2 = (J * J).sum(dim=1)
		Jn = torch.sqrt(Jn2)
		
		D = torch.diag_embed(1 / Jn)
		
		Js = torch.bmm(J,D)
		JsT = Js.transpose(1,2)
		
		Hs = torch.bmm(JsT,Js)
		gs = torch.bmm(JsT,res)
		
		u = torch.cholesky(Hs)
		q = torch.cholesky_solve(-gs,u)
		
		pGN = torch.bmm(D,q)
		
		p = torch.zeros(pGN.size()[0],pGN.size()[1],pGN.size()[2]).cuda()
		step = torch.zeros(pGN.size()[0]).cuda()
		
		mask1 = torch.norm(pGN,dim=1)[:,0] <= delta
		
		
		if mask1.sum() != 0:
			p[mask1] = pGN[mask1]
			step[mask1] = 0
		
		invD = torch.diag_embed(Jn)
		invD2 = torch.diag_embed(Jn2)
		
		invD2gs = torch.bmm(invD2,gs)
		invD2gsT = invD2gs.transpose(1,2)
		
		g = torch.bmm(invD,gs)
		gT = g.transpose(1,2)
		
		lambdaStar = torch.bmm(gT,g) / (torch.bmm(invD2gsT,invD2gs))
		
		CP = -lambdaStar*g
		
		mask2 = torch.logical_and(torch.norm(CP,dim=1)[:,0] > delta, torch.logical_not(mask1))
		if mask2.sum() != 0:
			p[mask2] = -(g[mask2] / torch.norm(g[mask2],dim=1)[:,None]) * delta[mask2][:,None,None]
			step[mask2] = 2
		
		not_mask = torch.logical_not(torch.logical_or(mask1,mask2))
		if not_mask.sum() != 0:
			A = (CP[not_mask] - pGN[not_mask]).pow(2).sum(dim=1)
			B = (2 * CP[not_mask] * (pGN[not_mask]-CP[not_mask])).sum(dim=1)
			C = CP[not_mask].pow(2).sum(dim=1) - delta[not_mask].pow(2).reshape(delta[not_mask].size()[0],1)
			
			k = (-B + torch.sqrt(B**2-4*A*C)/(2*A)).reshape(C.size()[0],1,1)
			
			p[not_mask] = CP[not_mask] + k*(pGN[not_mask] - CP[not_mask])
			step[not_mask] = 1
		
		return p, pGN, step
		
		
		
	def __levenberg_marquardt_powell(self, guess, data, dependent, delta, mu, eta, tol, max_iter):
		
		func = self.__get_func(dependent)
		fun = lambda p: func(p).sum(dim=1)
		
		res, resT = self.__get_residuals(func, guess, data)
		J, JT = self.__get_jacobians(fun, guess)
		f = 0.5 * torch.bmm(resT,res)
		
		
		#if True:
		#	srJ = 
		
		conv_perc = 0
		iteration = 0
		while True:
			p, pGN, step = self.__dogleg(res, J, delta)
			
			Jp = torch.bmm(J,p)
			JpT = Jp.transpose(1,2)
			
			converges = torch.logical_and((step == 0), self.__converges(torch.bmm(J,pGN), res, tol))
			
			t = guess + p.transpose(0,2)[0]
			rt, rtT = self.__get_residuals(func, t, data)
			ft = 0.5 * torch.bmm(rtT, rt)
			
			predicted = -torch.bmm(rtT,Jp)-0.5*torch.bmm(JpT,Jp)
			actual = f-ft
			
			rho = (actual / predicted).reshape(actual.size()[0])
			
			
			mask = rho <= mu
			not_mask = torch.logical_not(mask)
			if mask.sum() != 0:
				delta[mask] = 0.5 * delta[mask]
			
			pGN_Norm = torch.norm(pGN,dim=1)[:,0]
			mask = torch.logical_and(mask, delta > pGN_Norm)
			if mask.sum() != 0:
				delta[mask] = delta[mask] / (2 ** torch.ceil(torch.log2(delta[mask] / pGN_Norm[mask])))
			
			if not_mask.sum() != 0:
				print("changed guess: ", not_mask.sum())
				guess[:,not_mask] = t[:,not_mask]
			
			mask = rho >= eta
			if mask.sum() != 0:
				delta[mask] = delta[mask] * 2
			
			
			res, resT = self.__get_residuals(func, guess, data)
			J, JT = self.__get_jacobians(fun, guess)
			f = 0.5 * torch.bmm(resT,res)
			
			iteration += 1
			if iteration >= max_iter or converges.sum() == self.__npixels:
				conv_perc = float(converges.sum()) / float(self.__npixels)
				break
			
		
		return guess, conv_perc, iteration
		
	
	def solve(self, tol = 0.0001, max_iter=1000):
		guess = self.__init_guess
		data = self.__init_data
		dependent = self.__init_dependent
		
		delta0 = torch.norm(guess, dim=0)
		mu = 0.25
		eta = 0.75
		
		solution = self.__levenberg_marquardt_powell(guess, data, dependent, delta0, mu, eta, tol, max_iter)
		
		return solution
		
	
	
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


