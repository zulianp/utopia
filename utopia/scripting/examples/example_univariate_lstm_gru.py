import numpy as np # linear algebra
import torch
import torch.nn as nn
import time
import math, time
import utopia as ut 
import sys



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


def sequence_data(stock, lookback):
    data_raw = stock 
    data = []
    
    for index in range(len(data_raw) - lookback): 
        data.append(data_raw[index: index + lookback])
    
    data = np.array(data);
    test_set_size = int(np.round(0.2*data.shape[0]));
    train_set_size = data.shape[0] - (test_set_size);
    
    x_train = data[:train_set_size,:-1,:]
    y_train = data[:train_set_size,-1,:]
    
    x_test = data[train_set_size:,:-1]
    y_test = data[train_set_size:,-1,:]

    x_train = torch.from_numpy(x_train).type(torch.Tensor)
    x_test = torch.from_numpy(x_test).type(torch.Tensor)
    y_train = torch.from_numpy(y_train).type(torch.Tensor)
    y_test = torch.from_numpy(y_test).type(torch.Tensor)

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

def create_LSTM(input_dim, hidden_dim,  output_dim, num_layers):
    return LSTM(input_dim=input_dim, hidden_dim=hidden_dim, output_dim=output_dim, num_layers=num_layers)


def create_GRU(input_dim, hidden_dim,  output_dim, num_layers):
    return GRU(input_dim=input_dim, hidden_dim=hidden_dim, output_dim=output_dim, num_layers=num_layers)


def command_line_args(args):
    if len(args) != 8:
        print("8 inputs are reuquired")
        return
    return [str(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])]  


def finalize_evaluation(model, y_train, y_test, y_train_pred, close_min, close_max, name):
    y_test_pred = model_lstm(x_test)

    y_train_pred = denormalize_data(y_train_pred.detach().numpy(), close_min, close_max)
    y_train = denormalize_data(y_train.detach().numpy(), close_min, close_max)
    y_test_pred = denormalize_data(y_test_pred.detach().numpy(), close_min, close_max)
    y_test = denormalize_data(y_test.detach().numpy(), close_min, close_max)


    train_score = evaluate_predictions(y_train[:,0], y_train_pred[:,0]) 
    test_score = evaluate_predictions(y_test[:,0], y_test_pred[:,0])
    return  train_score, test_score
 

def train(model, criterion, y_train_lstm, x_train, lr, tol, name):
    start_time = time.time()
    hist = np.zeros(1000)
    lstm = []
    j = 0
    y_train_pred = model(x_train)
    loss = criterion(y_train_pred, y_train_lstm)
    while True:
        
        y_train_pred = model(x_train)
        loss = criterion(y_train_pred, y_train_lstm)

        
        if j > 2 and hist[j-1] == hist[j] or loss < tol or j == 100:
            training_time = time.time()-start_time
            training_time = int(training_time)

            return hist, y_train_pred, j, training_time


        print("Epoch ", j, "MSE: ", loss.item())
        hist[j] = loss.item()

        loss.backward()

        with torch.no_grad():
            i = 0
            for param in model.parameters():
          
                grad = param.grad.data
                grad_flat = grad.flatten()
       
                param_flat = param.data.flatten()
                param_flat_np = param_flat.detach().numpy().flatten()
         
                if i > 0:
                    momentum = param_flat
                    momentum_flat = momentum.flatten()

                param_np = param_flat.detach().numpy().flatten()
                param_ut = ut.Vector()
                param_ut.numpy_to_utopia(param_np)
                grad_np = grad_flat.detach().numpy().flatten()
                grad_ut = ut.Vector()
                grad_ut.numpy_to_utopia(grad_np)
            
                if i == 0:
                    param_ut.axpy(-lr, grad_ut)

                else:   
                    param_ut.axpy(-lr, grad_ut)
                    momentum_np = momentum_flat.detach().numpy().flatten()
                    momentum_ut = ut.Vector()
                    momentum_ut.numpy_to_utopia(momentum_np)
                    momentum_ut.scale(0.9)
                    param_ut.axpy(-1, momentum_ut)
              
                param_flat_np = np.ascontiguousarray(param_flat_np, dtype='float')
                param_ut.utopia_to_numpy(param_flat_np)
                param_flat = torch.from_numpy(param_flat_np)

                i += 1
    
            j += 1 
      
            
        for param in model.parameters():
            param.requires_grad = True
        
     
    return hist, y_train_pred    


# For example, it can be called with 
# python3 example_univariate_lstm_gru.py ./ibm.csv 1 10 2 1 0.1 0.15
########################## Body ###################################### 
ut.init()

args = sys.argv
filepath, input_dim, hidden_dim, num_layers, output_dim, lr, tol = command_line_args(args)
stock_name = filepath[2:-4]

open, high, low, close, volume = read_stock_data(filepath)
close_max, close_min, close = normalize_data(close.reshape(-1,1))

lookback = 21 
x_train, y_train, x_test, y_test = sequence_data(close, lookback)

model_lstm = create_LSTM(input_dim, hidden_dim,  output_dim, num_layers)
model_gru = create_LSTM(input_dim, hidden_dim,  output_dim, num_layers)
criterion = torch.nn.MSELoss(reduction='mean')

hist_lstm, y_train_pred_lstm, num_iterations_lstm, training_time_lstm = train(model_lstm, criterion, y_train, x_train, lr, tol, "LSTM")
hist_gru, y_train_pred_gru, num_iterations_gru, training_time_gru = train(model_gru, criterion, y_train, x_train, lr, tol, "GRU")

train_score_lstm, test_score_lstm = finalize_evaluation(model_lstm, y_train, y_test, y_train_pred_lstm, close_min, close_max, "LSTM")
train_score_gru, test_score_gru = finalize_evaluation(model_lstm, y_train, y_test, y_train_pred_gru, close_min, close_max, "GRU")

print("---------------------------------------------------------------------------------------------------")
print("CHOSEN  PARAMETRS")
print("Stock prices: ",stock_name, " | input dimension: ", input_dim, " | hidden dimension: ", hidden_dim, " | number of layers: ", num_layers)
print("learning rate: ", lr, " | tolerance: ", tol)
print("---------------------------------------------------------------------------------------------------")
print("Training time for LSTM: ", training_time_lstm, "s in", num_iterations_lstm, "iterations")
print("Train score for LSTM: ", train_score_lstm, ". Test score: ", test_score_lstm)
print("---------------------------------------------------------------------------------------------------")
print("Training time for GRU: ", training_time_gru, "s in", num_iterations_gru, "iterations")
print("Train score for GRU: ", train_score_gru, ". Test score: ", test_score_gru)
print("---------------------------------------------------------------------------------------------------")



# with open('out.txt', 'a') as f:
#     print("---------------------------------------------------------------------------------------------------", file=f)
#     print("CHOSEN  PARAMETRS", file=f)
#     print("Stock prices: ",stock_name, " | input dimension: ", input_dim, " | hidden dimension: ", hidden_dim, " | number of layers: ", num_layers, file=f)
#     # print("learning rate: ", lr, " | tolerance: ", tol, file=f)
#     # print("Training time for LSTM: ", training_time_lstm, "s in", num_iterations_lstm, "iterations",  file=f)
#     # print("Train score for LSTM: ", train_score_lstm, ". Test score: ", test_score_lstm, file=f)
#     print("---------------------------------------------------------------------------------------------------", file=f)

ut.finalize()

