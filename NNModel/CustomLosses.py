# Created by alexandra at 21/08/2023
import torch as t


class MSE:
    def __init__(self):
        super(MSE, self).__init__()
        self.mse = t.nn.MSELoss()

    def __call__(self, predicted_y, true_y):
        # mask = (true_y == 0)
        # mse_tensor = self.mse(predicted_y[~mask], true_y[~mask])
        mse_tensor = self.mse(predicted_y, true_y)
        mse_value = [float(mse_tensor.data.cpu().numpy())]
        return mse_tensor, mse_value

    # def mse(self, predicted_y, true_y):
    #     """ function that calculates the mean square error """
    #     predicted_y = t.transpose(predicted_y, 1, 2)
    #     true_y = t.transpose(true_y, 1, 2)
    #
    #     # mask = (true_backbone == 0)
    #     mse_tensor = self.mse(predicted_y, true_y)
    #     mse_value = [float(mse_tensor.data.cpu().numpy())]
    #     return mse_tensor, mse_value