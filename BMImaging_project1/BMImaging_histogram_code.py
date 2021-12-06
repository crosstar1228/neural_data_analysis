import os
import numpy as np
import nibabel as nib
from matplotlib import pyplot as plt

q = 0
w = 0
for i in range(5):
    q=q+1
    for j in range(2):
        w=w+1
        T1image=nib.load('0%d_0%d_T1ce.nii'%(q,w ))
        T1=T1image.get_data()
        T2image=nib.load('0%d_0%d_T2.nii'%(q,w))
        T2=T2image.get_data()
        mask_image=nib.load('0%d_0%d_mask.nii.gz'%(q, w))
        mask=mask_image.get_data()
        if w==2 :
            w = w-2     #01_01 부터 05_02까지 T1CE,T2, mask image 불러온 후 voxel data 값 불러오기

        len_axis=max(T1.max()+1,T2.max()+1)
        hist=np.zeros((len_axis,len_axis)) #histogram을 위한 영행렬 생성
        len_x=len(T1[:])
        len_y=len(T1[1][:])
        len_z=len(T1[1][1][:]) #for 문 구동을 위해 T1 intensity의 data값 길이
        for i in range(len_x):
            for j in range(len_y):
                for k in range(len_z):
                    if mask[i][j][k]!=0:
                        a=T1[i][j][k]
                        b=T2[i][j][k]
                        hist[a,b] = hist[a,b] + 1 #histogram 해당하는 행렬 구현

        x = []
        y = []
        for i in range(len_axis):
            for j in range(len_axis):
                for a in range(int(hist[i][j])):
                    x.append(i)
                for b in range(int(hist[i][j])):
                    y.append(j) #hist2d 함수 이용를 위해 x, y 성분에 각각 list 형태로 분리

        plt.hist2d(x, y, bins=(300,300),range=[[0,len_axis],[0,len_axis]])
        plt.show() #히스토그램 plotting