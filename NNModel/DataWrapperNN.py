# Created by alexandra at 08/09/2023
from torch.utils.data import Dataset, DataLoader
import torch as t
import numpy as np
from . import CONSTANTS as c


def dataWrapper(batch_size, num_workers, data, shuffle=False):
    dataset = MultiTaskDataset(data)
    data_generator = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, sampler=None,
                                num_workers=num_workers, collate_fn=multiTaskCollate)
    return data_generator


def batch_padding(dataset, batch_size):
    """
    function used to add random samples from dataset to have complete batches
    params: dataset(ndarray) - the dataset to be separates in batches
            batch_size(int) - model parameter representing the batch size
    return: dataset with (batch_size - len(dataset) % batch_size) number of samples added
    """
    mod = len(dataset) % batch_size
    if mod > 0:
        nr_of_samples = batch_size - mod
        for i in range(nr_of_samples):
            rnd = np.random.randint(len(dataset) - 1)
            # dataset = np.append(dataset, dataset[rnd])
            dataset.append(dataset[rnd])
    return dataset

def multiTaskCollate(batch):
    # TODO:
    nr_samples = len(batch)
    # print(batch)
    fragmentSize = []
    labels = []
    ctDNADetection = []
    ctDNAbyVAF = []
    ichorTF = []
    sampleNames = []
    for sample in batch:
        # sampleNames += sample[0]
        labels += [sample[1][1]]
        ctDNADetection += [sample[1][2]]
        ctDNAbyVAF += [sample[1][3]]
        ichorTF += [sample[1][4]]
        sampleNames += [sample[1][0]]
        fragmentSize += [t.from_numpy(np.array(sample[0], dtype=np.float64))]
    fragmentSize = t.stack(fragmentSize)
    labels = t.stack(labels)
    ctDNADetection = t.stack(ctDNADetection)
    ctDNAbyVAF = t.stack(ctDNAbyVAF)
    ichorTF = t.stack(ichorTF)


    return fragmentSize, sampleNames,  labels, ctDNADetection, ctDNAbyVAF, ichorTF

class MultiTaskDataset(Dataset):

    def __init__(self, x):
        self.fragmentSizes = []
        self.Y = []
        # tasks = x[1]
        # fragmentsFeatures = x[0]
        ## this function return two arrays, one with the unique sample name and another with the indexes, I am
        ## interested in the indexes, but also I am removing the last index, as for some reason besides the sample
        ## names I have also "PatientId";I do tasks[:, 0] because the unique function works only forone column
        # y_indexes = np.unique(tasks[:, PatientId], return_index=True)[1][:-1]
        for sample in x:
            sampleName = [sample[1].squeeze()[c.PatientId]]
            label = t.tensor(sample[1].squeeze()[c.Label], dtype=t.float64)
            ctDNADetection = t.tensor(float(sample[1].squeeze()[c.ctDNADetected]), dtype=t.float64)
            ctDNAbyVAF = t.tensor(sample[1].squeeze()[c.VAFg0p001], dtype=t.float64)
            ichorTF = t.tensor(sample[1].squeeze()[c.ichorTF], dtype=t.float64)
            self.fragmentSizes += [sample[0]]
            self.Y += [[sampleName, label, ctDNADetection, ctDNAbyVAF, ichorTF]]

    def __len__(self):
        return len(self.fragmentSizes)

    def __getitem__(self, idx):
        return self.fragmentSizes[idx], self.Y[idx]
