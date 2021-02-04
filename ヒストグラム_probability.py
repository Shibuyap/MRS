#coding:utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import os

startTime = time.time()

concestor = {0: 'Homo', 1: 'Hominini', 2: 'Homininae', 3: 'Hominoidae', 4: 'Hominoidea', 5: 'Catarrhini', 6: 'Simiiformes', 7: 'Primate_Strepsirrhini', 8: 'Euarchonta', 9: 'Euarchontoglires', 10: 'Boreoeutheria', 11: 'Eutheria', 12: 'Theria_Metatheria', 13: 'Mammalia', 14: 'Amniota_Sauropsida', 15: 'Tetrapoda', 16: 'Sarcopterygii', 17: 'Euteleostomii', 18: 'Vertebrata'}

'''
    エンハンサーを器官ごとに分け、すべてのcolumnを各node毎(19個)に分割して、
    それらの確率をそのままヒストグラムにしてプロット
'''


new_dir_path_recursive = './histgram_probability'
os.makedirs(new_dir_path_recursive, exist_ok=True)


# まずはall
new_dir_path_recursive = './histgram_probability/all'
os.makedirs(new_dir_path_recursive, exist_ok=True)

lists = [[] for i in range(19)]

for loop in range(33):
    inputName = './probabilities_out/probabilities_out_'
    if loop <= 8:
        inputName += '0'
    inputName += str(loop+1)
    inputName += '.txt'

    fIn = open(inputName,'r')

    datalist = fIn.readlines()
    n = len(datalist)

    print('loop = ', loop, end=" ")
    print('n = ', n)

    i = 0
    enhancerCount = 0
    while i < n:
        i += 1
        m = int(datalist[i])
        i += 1

        for j in range(m):
            i += 1
        m = int(datalist[i])
        i += 1

        for j in range(m):
            temp = datalist[i].split()
            
            for k in range(19):
                lists[18-k].append(float(temp[k])) # rootとヒトの整合性を整えるため(あとでもっとリファクタする)
            
            i += 1

        enhancerCount += 1
        # print(enhancerCount)

    fIn.close()

for ite in range(19):
    plt.figure()
    plt.hist(lists[ite], bins = 20)
    figName = './histgram_probability/all/all_'
    if ite <= 9:
        figName += '0'
    figName += str(ite)
    figName += '_'
    figName += concestor[ite]
    plt.title(concestor[ite])
    plt.savefig(figName)
    print('OK ', ite)
    plt.close()





# 次に各器官ごと

'''
    setに器官名を入れる
'''
inputAllTissue = './all_tissues.txt'
fAllTisuues = open(inputAllTissue,'r')
keepTissue = fAllTisuues.readlines()
tissueCount = int(keepTissue[0])
tissues = set()
for i in range(tissueCount):
    tissues.add(''.join(keepTissue[i+1].splitlines())) # ここで改行コードを取り除いている
    

print(len(tissues)) # 41

dic = {}
iota = 0
for name in tissues:
    dic[name] = iota
    iota += 1
print(dic)


for name in tissues:
    new_dir_path_recursive = './histgram_probability/'
    new_dir_path_recursive += name
    os.makedirs(new_dir_path_recursive, exist_ok=True)

    print(name)

    lists = [[] for i in range(19)]

    for loop in range(33):
        inputName = './probabilities_out/probabilities_out_'
        if loop <= 8:
            inputName += '0'
        inputName += str(loop+1)
        inputName += '.txt'

        fIn = open(inputName,'r')

        datalist = fIn.readlines()
        n = len(datalist)

        print('loop = ', loop, end=" ")
        print('n = ', n)

        i = 0
        enhancerCount = 0
        while i < n:
            i += 1 # chr~の行

            mm = int(datalist[i])
            i += 1 # 関係ある器官の数(0多し)
            flag = 0
            for j in range(mm):
                if ''.join(datalist[i].splitlines()) == name:
                    flag = 1
                i += 1
            
            m = int(datalist[i]) # 塩基の数(ヒトがgapのcolumnは除去済み)
            i += 1

            for j in range(m):
                if flag == 1:
                    temp = datalist[i].split()
                    for k in range(19):
                        lists[18-k].append(float(temp[k]))
                i += 1

            enhancerCount += 1
            # print(enhancerCount)

        fIn.close()

    for ite in range(19):
        plt.figure()
        plt.hist(lists[ite], bins = 20)
        figName = './histgram_probability/' + name + '/' + name + '_'
        if ite <= 9:
            figName += '0'
        figName += str(ite)
        figName += '_'
        figName += concestor[ite]
        plt.title(concestor[ite])
        plt.savefig(figName)
        print('OK ', ite)
        plt.close()


elapsed_time = time.time() - startTime
print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")
