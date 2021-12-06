import matplotlib.pyplot as plt
import nibabel as nib
import bct
from scipy.stats import pearsonr
import numpy as np

for A in range(1, 3):
    for B in range(1, 6):
        aaldata = nib.load('../Scripts/Data/0%d_0%d_aal.nii.gz' % (A, B))
        fmridata = nib.load('../Scripts/Data/0%d_0%d_fmri.nii.gz' % (A, B))

        aal_value = open('../scripts/Data/aal_values.txt')  # aal value txt 불러오기
        aal_list = []
        for line in aal_value:
            a = int(line.strip())
            aal_list.append(a)

        ROIlist = []
        for i in range(90):
            ROIlist.append([])

        time_avg_list = []
        random_list = []
        random_list2 = []

        aal = aaldata.get_data()
        fmri = fmridata.get_data()

        for i in range(128):
            for j in range(128):
                for k in range(35):
                    for J in range(90):
                        AAL = aal[i][j][k]

                        if AAL == aal_list[J]:
                            ROIlist[J].append((i, j, k))  # 각 aal value 해당하는 index 값의 모임.

        avg = 0
        for j in range(90):
            avg_list = []
            for k in range(97):
                sum = 0
                for i in range(len(ROIlist[j])):
                    x = ROIlist[j][i][0]
                    y = ROIlist[j][i][1]
                    z = ROIlist[j][i][2]
                    sum += fmri[x][y][z][k]
                    avg = sum / len(ROIlist[j])
                avg_list.append(avg)  # 1*97 fmri time series
            time_avg_list.append(avg_list)  # 90*97 aal value별 time series

        correlation_matrix = np.zeros([90, 90], dtype=np.float32)

        for i in range(90):
            for j in range(90):
                correlation_matrix[i][j] = pearsonr(time_avg_list[i], time_avg_list[j])[0]  # 90*90 connectivity matrix
                correlation_matrix[i][j] = abs(correlation_matrix[i][j])

        plt.imshow(abs(correlation_matrix))
        plt.colorbar()
        plt.clim(0, 1)
        plt.show()  # 양의 상관관계그래프 띄우기

        binary_matrix = correlation_matrix
        for i in range(90):
            for j in range(90):
                if binary_matrix[i][j] != 0:
                    if abs(binary_matrix[i][j]) < 0.1:  # threshold 0.1
                        binary_matrix[i][j] = 0
                    else:
                        binary_matrix[i][j] = 1  # binarized

        coefficient_vector = bct.clustering_coef_bu(binary_matrix)
        sum1 = 0
        for i in range(len(coefficient_vector)):
            sum1 += coefficient_vector[i]
        coefficient = sum1 / len(coefficient_vector)  # coefficient vector 평균내어 Cf 계산
        print(coefficient)
        charlength = bct.charpath(binary_matrix)  # Lf , lambda_ 제외한 나머지값들 주석처리했음
        print(charlength)

        for i in range(20):
            ran_matrix = bct.randmio_und(binary_matrix, 10)[0]
            random_coefficient_vector = bct.clustering_coef_bu(ran_matrix)
            sum2 = 0
            for i in range(len(random_coefficient_vector)):
                sum2 += random_coefficient_vector[i]
            ran_coefficient = sum2 / len(random_coefficient_vector)
            random_list.append(ran_coefficient)  # random한 20개 coefficent list
            ran_length = bct.charpath(ran_matrix)
            random_list2.append(ran_length)  # random한 20개 pathlength list

        sum3 = 0

        for i in range(len(random_list)):
            sum3 += random_list[i]
        ran_coefficient_avg = sum3 / len(random_list)  # Cr
        print(ran_coefficient_avg)
        sum4 = 0
        for i in range(len(random_list2)):
            sum4 += random_list2[i]
        ran_length_avg = sum4 / len(random_list2)  # Lr
        print(ran_length_avg)

        small_worldness = (coefficient / ran_coefficient_avg) / (charlength / ran_length_avg)  # Sf= (Cf/Cr)/(Lf/Lr)
        print(small_worldness)

