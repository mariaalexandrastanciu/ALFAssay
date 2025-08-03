# Created by alexandra at 28/07/2023
from torch import nn
import numpy as np
import math
import torch.nn.functional as F
import torch as t
import utils as u


class NeuralNetworkModel(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_layers, dev="cpu",
                 name="FragmentationPatterns", batch_size=1):
        super(NeuralNetworkModel, self).__init__()
        self.dev = dev
        self.name = name
        self.batch_size = batch_size
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.hidden_layer = hidden_layers


        super().__init__()  #100, 50, 20, 50, 100
        self.autoencoder = nn.Sequential(nn.Linear(input_dim, 100), nn.LayerNorm(100), nn.ReLU(),
                                             nn.Linear(100, 50), nn.LayerNorm(50), nn.ReLU(),
                                             nn.Linear(50, 20), nn.LayerNorm(20),  nn.ReLU(),
                                             nn.Linear(20, 50), nn.LayerNorm(50), nn.ReLU(),
                                             nn.Linear(50, 100), nn.LayerNorm(100), nn.ReLU(),
                                             nn.Linear(100, input_dim) ).to(dev)

        hidden_dim1 = input_dim + 100 #+ round(input_dim/10)   # 600
        hidden_dim2 = hidden_dim1 #round(hidden_dim1/3) # 300
        hidden_dim3 = hidden_dim1 #round(hidden_dim2*3) # 100
        self.forward_block = nn.Sequential(nn.Linear(input_dim, hidden_dim1), nn.LayerNorm(hidden_dim1),  nn.ReLU(),
                                                nn.Linear(hidden_dim1, hidden_dim2), nn.LayerNorm(hidden_dim2), nn.ReLU(),
                                                nn.Linear(hidden_dim2, hidden_dim3),  nn.LayerNorm(hidden_dim3), nn.ReLU(),
                                                nn.Linear(hidden_dim3, 1)
                                                ).to(dev)

    def forward(self, x):

        bias_correction = self.autoencoder(x.view(-1, self.input_dim))
        corrected_x = x.view(-1, self.input_dim)-bias_correction
        # tumor_fraction = self.forward_block(corrected_x.view(-1, self.input_dim)).squeeze()

        ## no correction
        tumor_fraction = self.forward_block(x.view(-1, self.input_dim)).squeeze()
        # result = per_sample.squeeze()


        return tumor_fraction, bias_correction, corrected_x




def save_grad(self, name):  # stuff you need to plot. leave it there
        def hook(grad): self.grads[name] = grad

        return hook

