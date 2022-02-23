#!/usr/bin/env python3
import math
import numpy as np
import torch
# import utopia as u 
from torch.optim.optimizer import Optimizer
from typing import Any, Callable, Dict, Iterable, Optional, Tuple, Union
from torch import Tensor, nn


Params = Union[Iterable[Tensor], Iterable[Dict[str, Any]]]
LossClosure = Callable[[], float]
OptLossClosure = Optional[LossClosure]
OptFloat = Optional[float]



class MyLinearRegression(nn.Module):
    def __init__(
        self,
        dtype: torch.dtype = torch.float
    ) -> None:
        super().__init__()
        self.a = nn.Parameter(torch.randn(1, requires_grad=True, dtype=dtype))
        self.b = nn.Parameter(torch.randn(1, requires_grad=True, dtype=dtype))

    def forward(self, x):
        return self.a*x + self.b



class GD:
    def __init__(
        self,
        model: nn.Module,
        lr: float = 1e-3,
        epochs: int = 100,
        tol: float = 1e-2
    ) -> None:
        self.model = model
        self.lr = lr
        self.epoch = 0
        self.max_epochs = epochs
        self.tol = tol
        
    def step(self, 
        x: Tensor,
        y: Tensor,
        loss_fn: LossClosure = nn.MSELoss(reduction='mean')
    ) -> None:
        if self.epoch >= self.max_epochs:
            return
         
        self.model.train()

        y_pred = model(x)

        loss_val = loss_fn(y, y_pred)
        if loss_val < self.tol:
            return

        loss_val.backward(retain_graph=True)
        print(loss_val.item())

        with torch.no_grad():
            for param in self.model.parameters():
                param -= self.lr * param.grad

        for param in self.model.parameters():
            param.requires_grad = True
        
        self.epoch += 1
    
    def train(self,
        x: Tensor,
        y: Tensor,
        loss_fn: LossClosure = nn.MSELoss(reduction='mean')
    ) -> None:
        for epoch in range(self.max_epochs):
            self.step(x, y, loss_fn)
        


dtype = torch.float
device = "cuda" if torch.cuda.is_available() else 'cpu'
device = torch.device(device)


x = torch.randn(100, requires_grad=True, dtype=torch.float, device=device)
y = 1 + 2*x + .1*torch.randn(100, requires_grad=True, dtype=torch.float, device=device)

model = MyLinearRegression().to(device)

newgd = GD(model)

newgd.train(x, y)
