# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:30:16 2022

@author: Shujaat
"""



import numpy as np
from tensorflow.keras.models import model_from_json
import tensorflow.keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Flatten, LSTM
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Conv1D, MaxPooling1D,AveragePooling1D
from tensorflow.keras.optimizers import SGD
#from tensorflow.keras.layers.wrappers import Bidirectional, TimeDistributed
from tensorflow.keras import regularizers
from tensorflow.keras import optimizers
from tensorflow.keras.layers import Input, BatchNormalization
from tensorflow.keras.models import Model
from sklearn import metrics
import tensorflow as tf
from tensorflow.keras import regularizers
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from tensorflow.keras import initializers
from tensorflow.keras.layers import Activation, Dense, Add
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
import pandas as pd
from tensorflow.keras.utils import to_categorical
from sklearn.model_selection import learning_curve
from sklearn import metrics
from sklearn.metrics import auc
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from sklearn.model_selection import StratifiedKFold
from keras.layers import GaussianNoise

def get_model():

    input_shape = (81,4)
    inputs = Input(shape = input_shape)
    gl=(GaussianNoise(0.2))(inputs)
    convLayer = Conv1D(filters = 128, kernel_size = 5,activation = 'relu',kernel_regularizer = regularizers.l2(1e-4), bias_regularizer = regularizers.l2(1e-6),input_shape = input_shape)(gl)#(inputs)
    poolingLayer = MaxPooling1D(pool_size = 4, strides=2)(convLayer)
    dropoutLayer = Dropout(0.5)(poolingLayer)
    convLayer2 = Conv1D(filters = 256, kernel_size = 7,activation = 'relu',kernel_regularizer = regularizers.l2(1e-4), bias_regularizer = regularizers.l2(1e-6))(dropoutLayer)
    poolingLayer2 = MaxPooling1D(pool_size = 4, strides=2)(convLayer2)
    dropoutLayer2 = Dropout(0.5)(poolingLayer2)
    flattenLayer = Flatten()(dropoutLayer2)
    denseLayer2 = Dense(64, activation = 'relu',kernel_regularizer = regularizers.l2(1e-4),bias_regularizer = regularizers.l2(1e-6))(flattenLayer)#(dropoutLayer4)
    dropoutLayer3 = Dropout(0.5)(denseLayer2)
    outLayer = Dense(1, activation='sigmoid')(dropoutLayer3)
    model2 = Model(inputs = inputs, outputs = outLayer)
    #optimizer= SGD(momentum = 0.96, lr = 0.0077)
    model2.compile(loss='binary_crossentropy',optimizer= "adam", metrics=['binary_accuracy']);

    return model2
#%%

modelProMN = get_model()
modelProMN.load_weights('D:/Research papers/iProm-Sigma 54/model_sigme54_weights.h5')

#modelSecond=get_model()
#modelSecond.load_weights('D:/Research papers/Phage Promoter Review/Final Model and Results/SecondLayerWeights.h5')
#%%


def encode_seq(s):
    Encode = {'A':[1,0,0,0],'C':[0,1,0,0],'G':[0,0,1,0],'T':[0,0,0,1],'N':[0,0,0,0]}
    return np.array([Encode[x] for x in s])
#%%

X1 = {}
accumulator=0

Strng2=["AAGATAGGCGTTGACTTGATGGGTCTTTAGGTGTAGGCTTTAGGTGTTGGCTTTAG"]
def iProm_phage(Strng):
    p=0
    n=0
    prediction=""
    X1 = {}
    my_hottie = encode_seq((Strng))
    out_final=my_hottie
    out_final = np.array(out_final)
    X1[accumulator]=out_final
      #out_final=list(out_final)
    X1[accumulator] = out_final    
    X1 = list(X1.items()) 
    an_array = np.array(X1)
    an_array=an_array[:,1]    
    transpose = an_array.T
    transpose_list = transpose.tolist()
    X1=np.transpose(transpose_list)
    X1=np.transpose(X1)
    pr=modelProMN.predict(X1)
    pr=pr.round()
    if(pr==1):
        print('Query Sequence is promoter')
        p=1
        # prS70=modelSecond.predict(X1)
        # prS70=prS70.round()
        # if(prS70==1):
        #     print('& Phage Promoter')
        #     p=1
        # else:
        #     print('& Host Promoter ')
        #     n=1
    else:
       print('Query Sequence is non promoter')
       n=1
    return p,n

predictions=[]
def predict_seq(seqs):
    for i in range(len(seqs)):
        pred=iProm_phage(seqs[i])
        predictions.append(pred)
    return predictions

#%%
       
sequences = [] 
for record in SeqIO.parse("D:/Research papers/iProm-Sigma 54/Dataaset Pro54/Prom54 Promoter Dataset.txt", "fasta"):
    sequences.append((record.seq.upper()))

indexes=[]
PosSeq=[]
lengths=[]
seqs=[]
lent=81
for i in range (len(sequences)):
    c=sequences[i]
    b=c._data
    #b=b.decode(encoding="utf-8")
    if(len(b)==81):
        #newNeg.append(b)
        PosSeq.append(b)
    else:
        ac=list(b)
        for j in range (len(ac),lent):
            ac.append('N')
        b2 = ''.join(ac)
        seqs.append(b2)
PositiveSeq=PosSeq+seqs
PositiveSeq=PositiveSeq[50:200]
#%%
score=[]    
for i in range(len(PositiveSeq)):
    e=iProm_phage(PositiveSeq[i])
    score.append(e)

