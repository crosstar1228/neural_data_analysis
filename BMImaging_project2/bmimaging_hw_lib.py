# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import mpmath


'''
Written by Seung-Wook Kim on May 4, 2018.

Reference:
   https://www.mristudio.org/wiki/faq
'''

class fiber():
    def __init__(self):
        self.nLength=0
        self.nSelStatus=[]
        self.rgbColor=[]
        self.nSelLenStart=0
        self.nSelLenEnd=0
        self.pxyzChain=[]

class Fibers():
    def __init__(self):
        self.FiberNr=0
        self.FiberLenMax=0
        self.FiberLenMean=0
        self.ImgWidth=0
        self.ImgHeight=0
        self.ImgSlice=0
        self.PixelSizeWidth=0
        self.PixelSizeHeight=0
        self.SliceThickness=0
        self.SliceOrientation=0
        self.SliceSequencing=0
        self.Fiber=[]

'''
use read_fiber to load the fibers from .dat file.
output.Fiber contains the fiber informations.
'''
def read_fiber(file):
    new_fiber = Fibers()
    fid=open(file)

    for i in range(2):
        chk1=np.fromfile(fid,np.dtype(np.int32),1)
    
    new_fiber.FiberNr = np.fromfile(fid,np.dtype(np.int32),1)
    new_fiber.FiberLenMax = np.fromfile(fid,np.dtype(np.int32),1)
    new_fiber.FiberLenMean = np.fromfile(fid,np.dtype(np.float32),1)
    new_fiber.ImgWidth = np.fromfile(fid,np.dtype(np.int32),1)
    new_fiber.ImgHeight = np.fromfile(fid,np.dtype(np.int32),1)
    new_fiber.ImgSlice = np.fromfile(fid,np.dtype(np.int32),1)
    new_fiber.PixelSizeWidth = np.fromfile(fid,np.dtype(np.float32),1)
    new_fiber.PixelSizeHeight = np.fromfile(fid,np.dtype(np.float32),1)
    new_fiber.SliceThickness = np.fromfile(fid,np.dtype(np.float32),1)
    new_fiber.SliceOrientation = np.fromfile(fid,np.dtype(np.int32),1)
    new_fiber.SliceSequencing = np.fromfile(fid,np.dtype(np.int32),1)
    
    for i in range(19):
        chk2=np.fromfile(fid,np.dtype(np.int32),1)
        
    for i in range(new_fiber.FiberNr[0]):
        n=fiber()
        n.nLength = np.fromfile(fid,np.dtype(np.int32),1)
        for j in range(3):
            chk3=np.fromfile(fid,np.dtype(np.int32),1)
        
        for j in range(n.nLength[0]):
            pxyz = np.fromfile(fid,np.dtype(np.float32),3)
            n.pxyzChain.append(pxyz)
            
        new_fiber.Fiber.append(n)

    return new_fiber

'''
use load_img to load .img files you saved from DTI Studio
'''
def load_img(file):
    fid=open(file,'rb')
    data=np.fromfile(fid,np.dtype('float32'))
    fa=data.reshape((70,128,128))
    
    return fa

'''
use load_nifti to load given aal.nii.gz file.
if you don't use this function and load the image yourself, 
then the image direction and dimension might be wrong with others.
'''
def load_nifti(file):
    aal=nib.load(file)
    aimg=aal.get_data()
    aimg=np.rot90(aimg,3)
    
    cimg=np.zeros([70,128,128],dtype=np.int32)
    for i in range(70):
        cimg[i,:,:]=aimg[:,:,i]
        
    return cimg

for q in range(1,3):
    for w in range(1,6):#01_01~02_05??????
        fadata=load_img('../scripts/DWI_data/0%d_0%d_Image.img'%(q,w))
        aaldata=load_nifti('../scripts/DWI_data/0%d_0%d_aal.nii.gz'%(q,w))
        fiberdata=read_fiber('../scripts/DWI_data/0%d_0%d_Fiber.dat'%(q,w)) # fa,fiber,aal data ????????????

        aal_value = open('../scripts/DWI_data/aal_values.txt')#aal value txt ????????????
        aal_list=[]
        for line in aal_value :
            a=int(line.strip())
            aal_list.append(a) #????????????list??? ??????


        start_list=[0]*len(fiberdata.Fiber)
        end_list=[0]*len(fiberdata.Fiber) #fiber??????,????????? intensity??? list
        start_point = [0] * len(fiberdata.Fiber)
        end_point = [0] * len(fiberdata.Fiber)#fiber?????? ??????????????? ??????,????????? ????????? list
        length_list=[0]*len(fiberdata.Fiber) #?????? list
        fa_sum_list=[0]*len(fiberdata.Fiber) #fa?????? ??? list
        fa_avg_list=[0]*len(fiberdata.Fiber) #fa?????? ?????? list ??? fa ???/??????

        '''fiber????????? intensity, ?????? ????????? for??? '''
        for I in range(len(fiberdata.Fiber) ):
            for i in range(len(fiberdata.Fiber[I].pxyzChain)):
                AAL = aaldata[int(fiberdata.Fiber[I].pxyzChain[i][2])][int(fiberdata.Fiber[I].pxyzChain[i][1])][int(fiberdata.Fiber[I].pxyzChain[i][0])]
                idx=0
                for j in range(90):
                    if AAL == aal_list[j]:
                        idx=AAL
                if idx!=0:
                    start_list[I]=idx #intensity ??????
                    start_point[I]=i  #?????? ??????
                    break

        '''fiber?????? intensity, ?????? ????????? for??? '''
        for k in range(len(fiberdata.Fiber)):
            for i in range(len(fiberdata.Fiber[k].pxyzChain)-1,-1,-1):
                AAL2 = aaldata[int(fiberdata.Fiber[k].pxyzChain[i][2])][int(fiberdata.Fiber[k].pxyzChain[i][1])][int(fiberdata.Fiber[k].pxyzChain[i][0])]
                idx2=0
                for j in range(90):
                    if AAL2 == aal_list[j]:
                        idx2=AAL2
                if idx2!=0:
                    end_list[k]=idx2 #intensity ??????
                    end_point[k]=i  #?????? ??????
                    break

        '''?????????,?????? ?????? ???????????? fiber?????? ????????? for??? '''
        for m in range(len(fiberdata.Fiber)):
            length=0
            for n in range(start_point[m],end_point[m]):
                z_length=abs(fiberdata.Fiber[m].pxyzChain[n][2]-fiberdata.Fiber[m].pxyzChain[n+1][2])
                y_length=abs(fiberdata.Fiber[m].pxyzChain[n][1]-fiberdata.Fiber[m].pxyzChain[n+1][1])
                x_length=abs(fiberdata.Fiber[m].pxyzChain[n][0]-fiberdata.Fiber[m].pxyzChain[n+1][0])#z,x,y?????? ?????? ??????
                length+=mpmath.sqrt((1.718*x_length)**2+(1.718*y_length)**2+(2*z_length)**2)
            if float(length)>=20: #20mm??????????????? ??????
                length_list[m]=float(length)

        '''fa ??? ??? fa?????? ????????? for??? '''
        for m in range(len(fiberdata.Fiber)):
            fa_sum = 0
            for n in range(len(fiberdata.Fiber[m].pxyzChain)):
                fa_sum+=fadata[int(fiberdata.Fiber[m].pxyzChain[n][2])][int(fiberdata.Fiber[m].pxyzChain[n][1])][int(fiberdata.Fiber[m].pxyzChain[n][0])]
            fa_sum_list[m]=fa_sum #fa sum
            if length_list[m]!=0:
                fa_avg_list[m]=fa_sum_list[m]/length_list[m]#length??? ?????? fa??????

        '''fa value ?????? ????????? '''
        matrix = np.zeros([90,90],dtype=np.float32) #90*90 ????????? ??????
        for n in range(len(fiberdata.Fiber)):
            for i in range(90):
                if start_list[n]==aal_list[i]:
                    I=i
                    break
                else:
                    continue #intensity ??? ???????????? ?????? ?????? ?????? ???????????? ??? ??????
            for j in range(90):
                if end_list[n]==aal_list[j]:
                    J=j
                    break
                else:
                    continue  #intensity ??? ???????????? ?????? ?????? ?????? ???????????? ??? ??????
            favalue=fa_avg_list[n]
            if favalue!=0:
                a=matrix[I][J]
                a+=favalue
                matrix[I][J]=a #????????? fa??? ????????????
                matrix[J][I]=a #?????? symmetric ?????? ???????????????
        for i in range(90):
            a=matrix[i][i]
            a=0
            matrix[i][i]=a #?????? ?????? ?????? 0?????? ???????????????
        plt.imshow(matrix)
        plt.colorbar()
        plt.clim(0,80) #???????????? ??????????????? ????????? ????????? 80??????
        plt.show() #connectivity diagram ??????

        '''degree distribution ???????????? ????????? ??????'''
        fa_degree=[]
        for i in range(90):
            line_sum = 0
            for j in range(90):
                line_sum += matrix[i][j]
            fa_degree.append(line_sum) #?????? ?????? ????????? append

        histogram = plt.figure()
        plt.hist(fa_degree, bins=5, range=[0, 400])
        plt.yticks(range(0, 55, 5)) #????????? ???????????? ????????? ????????? ????????????
        plt.title("degree distribution")
        plt.xlabel("degree")
        plt.ylabel('number of node')
        plt.show()  #??????????????? ??????









