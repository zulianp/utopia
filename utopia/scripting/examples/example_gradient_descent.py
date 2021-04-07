import math
import numpy as np
import torch
from torch.optim.optimizer import Optimizer
from typing import Any, Callable, Dict, Iterable, Optional, Tuple, Union
from torch import Tensor


Params = Union[Iterable[Tensor], Iterable[Dict[str, Any]]]
LossClosure = Callable[[], float]
OptLossClosure = Optional[LossClosure]
OptFloat = Optional[float]


class MyOptimizerBase(Optimizer):     
  def __init__(self,params, defaults):
        super(MyOptimizerBase, self).__init__(params, defaults)


  def optimize_history(self, fun, x, tensor_size, plt_fun):
    num_iter=self.param_groups[0]['num_iter']
    g_tol=self.param_groups[0]['g_norm_tol']

    print("it     f     ||g||")

    iters = np.zeros((tensor_size, num_iter))
    iters[:,0] = x.detach().numpy()    
    for it in range(0, num_iter):
        self.zero_grad()
        f = fun(x)
        f.backward(create_graph=True, retain_graph=True)
        torch.nn.utils.clip_grad_norm_(x, 1.0)
        self.step()
        iters[:, it] = x.detach().numpy()    

        norm_p=0;
        for p in self.param_groups[0]['params']:
            norm_p=norm_p+np.linalg.norm(p.grad.data.detach().numpy())

        if(norm_p < g_tol):
            print("Solver converged at: ", it, "  iteration.")
            break

        print(it, "     ", f.detach().numpy(), "     ",  norm_p)
    

    plt_fun(iters)


  def optimize(self, fun, x):
    num_iter=self.param_groups[0]['num_iter']
    g_tol=self.param_groups[0]['g_norm_tol']

    print("it        f           ||g||")


    for it in range(0, num_iter):
        self.zero_grad()
        f = fun(x)
        f.backward(create_graph=True, retain_graph=True)
        torch.nn.utils.clip_grad_norm_(x, 1.0)
        self.step()
        
        norm_p=0;
        for p in self.param_groups[0]['params']:
            norm_p=norm_p+np.linalg.norm(p.grad.data.detach().numpy())

        # this should not be here.... at the end of the day 
        if(norm_p < g_tol):
            print("Solver converged at: ", it, "  iteration.")
            break

        print(it, "     ", f.detach().numpy(), "     ",  norm_p)
    



class GD(MyOptimizerBase):     
  def __init__(
        self,
        params: Params,
        lr: float = 1e-3,
        num_iter: int=100,
        g_norm_tol: float =1e-4
    ) -> None:


        defaults = dict(
            lr=lr,
            num_iter=num_iter,
            g_norm_tol=g_norm_tol, 
        )
        super(GD, self).__init__(params, defaults)



  def step(self, closure: OptLossClosure = None) -> OptFloat:
        loss = None
        if closure is not None:
            loss = closure()
            # print("loss ", loss)

        for group in self.param_groups:

            # print('group: ', group)

            lr=group['lr']

            num_iter=group['num_iter']


            for p in group['params']:
                if p.grad is None:
                    continue

                grad = p.grad.data

                # to get gradient 
                # grad_d_n = grad.detach().numpy()
                # print("grad_d_n: ", grad_d_n)                

                # add gradient 
                p.data.add_(grad, alpha=-lr)

        return loss