import numpy as np # linear algebra
import torch
import torch.nn as nn
import time
import math, time
from conversion_functions import *
import utopia as ut 



########################## Classes ########################## 
class LSTM(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, output_dim):
        super(LSTM, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_dim).requires_grad_()
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_dim).requires_grad_()
        out, (hn, cn) = self.lstm(x, (h0.detach(), c0.detach()))
        out = self.fc(out[:, -1, :]) 
        return out
        

    
class GRU(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, output_dim):
        super(GRU, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        
        self.gru = nn.GRU(input_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_dim).requires_grad_()
        out, (hn) = self.gru(x, (h0.detach()))
        out = self.fc(out[:, -1, :]) 
        return out
    

########################## Methods ##########################
def read_stock_data(filepath):
    data= np.genfromtxt(filepath, delimiter=',', skip_header=1, dtype=np.float32)
    open = data[:,1]
    high = data[:,2]
    low = data[:,3]
    close = data[:,4]
    volume = data[:,5]
        
    return open, high, low, close, volume


def split_data(stock, lookback):
    data_raw = stock # convert to numpy array
    data = []
    
    # create all possible sequences of length seq_len
    for index in range(len(data_raw) - lookback): 
        data.append(data_raw[index: index + lookback])
    
    data = np.array(data);
    test_set_size = int(np.round(0.2*data.shape[0]));
    train_set_size = data.shape[0] - (test_set_size);
    
    x_train = data[:train_set_size,:-1,:]
    y_train = data[:train_set_size,-1,:]
    
    x_test = data[train_set_size:,:-1]
    y_test = data[train_set_size:,-1,:]
  
    return [x_train, y_train, x_test, y_test]


def normalize_data(data):
    data_max, data_min = data.max(), data.min()
    norm_data = (data - data_min)/(data_max - data_min)
    
    return data_max, data_min, norm_data

def denormalize_data(norm_data, min, max): 
    return norm_data*(max-min) + min
 

def evaluate_predictions(y_true, y_pred):
    assert y_true.shape == y_pred.shape
    return ((y_true - y_pred) ** 2).mean()

########################## Body ########################## 
ut.init()
filepath = './data/ibm.csv'
open, high, low, close, volume = read_stock_data(filepath)
close_max, close_min, close = normalize_data(close.reshape(-1,1))

# a financial year has ~ 252 days, so 21 represents the days in one month
lookback = 21 
x_train, y_train, x_test, y_test = split_data(close, lookback)


x_train = torch.from_numpy(x_train).type(torch.Tensor)
x_test = torch.from_numpy(x_test).type(torch.Tensor)
y_train_lstm = torch.from_numpy(y_train).type(torch.Tensor)
y_test_lstm = torch.from_numpy(y_test).type(torch.Tensor)
y_train_gru = torch.from_numpy(y_train).type(torch.Tensor)
y_test_gru = torch.from_numpy(y_test).type(torch.Tensor)

input_dim = 1
hidden_dim = 32
num_layers = 2
output_dim = 1
num_epochs = 100

model = LSTM(input_dim=input_dim, hidden_dim=hidden_dim, output_dim=output_dim, num_layers=num_layers)
criterion = torch.nn.MSELoss(reduction='mean') # TODO: utopia? 
# optimiser = torch.optim.Adam(model.parameters(), lr=0.01)  TODO: utopia?  

hist = np.zeros(num_epochs)
start_time = time.time()
lstm = []


for t in range(num_epochs):
    y_train_pred = model(x_train)

    loss = criterion(y_train_pred, y_train_lstm)
    print("Epoch ", t, "MSE: ", loss.item())
    hist[t] = loss.item()

    #optimiser.zero_grad()
    loss.backward()
    #optimiser.step()
    with torch.no_grad():
        #  make gradient step x_{n + 1} = x_n - lr*gradient
        i = 0
        for param in model.parameters():
          
            grad = param.grad.data
            grad_flat = grad.flatten()
       
            param_flat = param.data.flatten()
         
            if i > 0:
                momentum = param_flat
                momentum_flat = momentum.flatten()
                print(momentum_flat.shape)

            p_ut = pytorch_to_utopia(param_flat)
            grad_ut = pytorch_to_utopia(grad_flat)
            
            # # gradient step
            if i == 0:
                p_ut.axpy(-0.01, grad_ut)

            else: 
                p_ut.axpy(-0.01, grad_ut)
                momentum_ut = pytorch_to_utopia(momentum_flat)
                momentum_ut.scale(0.9)
                p_ut.axpy(-1, momentum_ut)
              

            utopia_to_pytorch(p_ut, param_flat)
     
            i += 1
            
        # allow gradient on model parameters
    for param in model.parameters():
        param.requires_grad = True
    
training_time = time.time()-start_time
print("Training time: {}".format(training_time))

predict = denormalize_data(y_train_pred.detach().numpy(), close_min, close_max)
original = denormalize_data(y_train_lstm.detach().numpy(), close_min, close_max)

y_test_pred = model(x_test)

y_train_pred = denormalize_data(y_train_pred.detach().numpy(), close_min, close_max)
y_train = denormalize_data(y_train_lstm.detach().numpy(), close_min, close_max)
y_test_pred = denormalize_data(y_test_pred.detach().numpy(), close_min, close_max)
y_test = denormalize_data(y_test_lstm.detach().numpy(), close_min, close_max)


train_score = evaluate_predictions(y_train[:,0], y_train_pred[:,0])
print('LSTM train Score: %.2f' % (train_score))
test_score = evaluate_predictions(y_test[:,0], y_test_pred[:,0])
print('LSTM test Score: %.2f' % (test_score))
ut.finalize()





