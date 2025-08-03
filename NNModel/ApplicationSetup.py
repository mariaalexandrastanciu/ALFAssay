# Created by alexandra at 28/07/2023

"""
all application setup should be in here
"""
import torch


def application_setup():
    torch.set_default_dtype(torch.float64)
    pass


# dev_to_use=torch.device('cpu') #change here to use GPU or CPU
def if_cuda():
    if torch.cuda.is_available():
        dev = "cuda"

    else:
         dev = "cpu"
    return dev
