#coding:utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import pickle

startTime = time.time()


modeRootCut = 0


'''
    エンハンサーを関係する器官ごと(41個)に分け、
'''

humanPathSize = 19
concestor = {0: 'Homo', 1: 'Hominini', 2: 'Homininae', 3: 'Hominoidae', 4: 'Hominoidea', 5: 'Catarrhini', 6: 'Simiiformes', 7: 'Primate_Strepsirrhini', 8: 'Euarchonta', 9: 'Euarchontoglires', 10: 'Boreoeutheria', 11: 'Eutheria', 12: 'Theria_Metatheria', 13: 'Mammalia', 14: 'Amniota_Sauropsida', 15: 'Tetrapoda', 16: 'Sarcopterygii', 17: 'Euteleostomii', 18: 'Vertebrata'}
concestorLabel = ['Homo', 'Hominini', 'Homininae', 'Hominoidae', 'Hominoidea', 'Catarrhini', 'Simiiformes', 'Primate_Strepsirrhini', 'Euarchonta', 'Euarchontoglires', 'Boreoeutheria', 'Eutheria', 'Theria_Metatheria', 'Mammalia', 'Amniota_Sauropsida', 'Tetrapoda', 'Sarcopterygii', 'Euteleostomii', 'Vertebrata']
concestorJapanese = {0: 'ヒト', 1: 'Hominini', 2: 'Homininae', 3: 'Hominoidae', 4: 'Hominoidea', 5: 'Catarrhini', 6: 'Simiiformes', 7: 'Primate_Strepsirrhini', 8: 'Euarchonta', 9: 'Euarchontoglires', 10: 'Boreoeutheria', 11: 'Eutheria', 12: 'Theria_Metatheria', 13: 'Mammalia', 14: 'Amniota_Sauropsida', 15: 'Tetrapoda', 16: 'Sarcopterygii', 17: 'Euteleostomii', 18: 'Vertebrata'}
inputName = './lengthOfHumanPathEdges.txt'
fIn = open(inputName,'r')
lengthOfHumanPathEdges = [float(x.strip()) for x in fIn.readline().split()]


'''
    setに器官名を入れる
'''
inputAllTissue = './all_tissues.txt'
fAllTisuues = open(inputAllTissue,'r')
keepTissue = fAllTisuues.readlines()
tissueCount = int(keepTissue[0]) # 41
dic = {}
for i in range(tissueCount):
    dic[''.join(keepTissue[i+1].splitlines())] = i # ここで改行コードを取り除いている



inputName = './countAllArgMax.txt'
fIn = open(inputName,'r')
countAllArgMax = [int(x.strip()) for x in fIn.readline().split()]
fIn.close()

inputName = './countTissuesArgMax.txt'
fIn = open(inputName,'r')
countTissues = []
for i in range(tissueCount):
    temp = [int(x.strip()) for x in fIn.readline().split()]
    countTissues.append(temp)
fIn.close()



colorlist = ["r", "g", "b", "c", "m", "y", "k"]


'''
    横軸に器官、縦軸に比率のグラフをconcestor interval毎に描く
'''


label = ["adipose_tissue", "blood", "blood_vessel", "brain", "esophagus", "eye", "female_gonad", "gallbladder", "heart", "internal_male_genitalia", "kidney", "large_intestine", "liver", "lung", "lymph_node", "meninx", "olfactory_region", "pancreas", "parotid_gland", "penis", "placenta", "prostate_gland", "salivary_gland", "skeletal_muscle_tissue", "skin_of_body", "small_intestine", "smooth_muscle_tissue", "spinal_cord", "spleen", "stomach", "submandibular_gland", "testis", "throat", "thymus", "thyroid_gland", "tongue", "tonsil", "umbilical_cord", "urinary_bladder", "uterus", "vagina"]


ratioAll = [0.0] * humanPathSize
sumAll = sum(countAllArgMax)
for i in range(humanPathSize):
    ratioAll[i] = countAllArgMax[i] / sumAll


sumTissues = [0] * tissueCount
ratioTissues = [[0.0] * humanPathSize for i in range(tissueCount)]
for i in range(tissueCount):
    sumTissues[i] = sum(countTissues[i])
    for j in range(humanPathSize):
        ratioTissues[i][j] = countTissues[i][j] / sumTissues[i]



# 横軸Interbal、縦軸Countの組織ごと(Allも)のヒストグラム
if False:
    print('x axis:Interval, y axis: Count, for all tissues')
    new_dir_path_recursive = './intervalCount'
    os.makedirs(new_dir_path_recursive, exist_ok=True)
    # All
    left = list(range(humanPathSize))
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    plt.bar(left, countAllArgMax, color=colorlist, tick_label=concestorLabel, align="center")
    plt.xlabel('Concestor')
    plt.ylabel('counts')
    plt.xticks(rotation=90)
    figName = './intervalCount/intervalCount_all'
    pltTitle = 'all'
    plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('all OK')

    # 各Tissue
    for i in range(tissueCount):
        left = list(range(humanPathSize))
        plt.figure(figsize=(16.0,12.0))
        plt.bar(left, countTissues[i], color=colorlist, tick_label=concestorLabel, align="center")
        plt.xlabel('Concestor')
        plt.ylabel('counts')
        plt.xticks(rotation=90)
        figName = './intervalCount/intervalCount_'
        figName += label[i]
        pltTitle = label[i]
        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(label[i],' OK')

# 横軸Interbal、縦軸Countの組織ごと(Allも)のヒストグラム
# MRS入らないデータは除く
if False:
    print('x axis:Interval, y axis: Count, for all tissues 2')
    new_dir_path_recursive = './intervalCount2'
    os.makedirs(new_dir_path_recursive, exist_ok=True)
    # All
    left = list(range(humanPathSize-1))
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    countAllArgMax2 = []
    for i in range(len(countAllArgMax)-1):
        countAllArgMax2.append(countAllArgMax[i])
    concestorLabel2 = []
    for i in range(len(concestorLabel)-1):
        concestorLabel2.append(concestorLabel[i])
    plt.bar(left, countAllArgMax2, color=colorlist, tick_label=concestorLabel2, align="center")
    plt.xlabel('Concestor')
    plt.ylabel('counts')
    plt.xticks(rotation=90)
    figName = './intervalCount2/intervalCount2_all'
    pltTitle = 'all'
    plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('all OK')

    # 各Tissue
    for i in range(tissueCount):
        left = list(range(humanPathSize-1))
        plt.figure(figsize=(16.0,12.0))
        plt.rcParams["font.size"] = 18
        countTissues2 = []
        for j in range(len(countTissues[i])-1):
            countTissues2.append(countTissues[i][j])
        plt.bar(left, countTissues2, color=colorlist, tick_label=concestorLabel2, align="center")
        plt.xlabel('Concestor')
        plt.ylabel('counts')
        plt.xticks(rotation=90)
        figName = './intervalCount2/intervalCount2_'
        figName += label[i]
        pltTitle = label[i]
        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(label[i],' OK')




# 横軸Interbal、縦軸Countの組織ごと(Allも)のヒストグラムを枝長で割ったグラフ
# MRS入らないデータは0除算になるので除く
if False:
    print('x axis:Interval, y axis: Count/edgeLength, for all tissues')
    new_dir_path_recursive = './intervalCountPerLength'
    os.makedirs(new_dir_path_recursive, exist_ok=True)
    # All
    left = list(range(humanPathSize-1))
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    countAllPerLength = []
    for i in range(len(countAllArgMax)-1):
        countAllPerLength.append(countAllArgMax[i] / lengthOfHumanPathEdges[i])
    concestorLabel2 = []
    for i in range(len(concestorLabel)-1):
        concestorLabel2.append(concestorLabel[i])
    plt.bar(left, countAllPerLength, color=colorlist, tick_label=concestorLabel2, align="center")
    plt.xlabel('Concestor')
    plt.ylabel('counts per edge length')
    plt.xticks(rotation=90)
    figName = './intervalCountPerLength/intervalCountPerLength_all'
    pltTitle = 'all'
    plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('all OK')

    # 各Tissue
    for i in range(tissueCount):
        left = list(range(humanPathSize-1))
        plt.figure(figsize=(16.0,12.0))
        plt.rcParams["font.size"] = 18
        countTissuesPerLength = []
        for j in range(len(countTissues[i])-1):
            countTissuesPerLength.append(countTissues[i][j] / lengthOfHumanPathEdges[j])
        plt.bar(left, countTissuesPerLength, color=colorlist, tick_label=concestorLabel2, align="center")
        plt.xlabel('Concestor')
        plt.ylabel('counts per edge length')
        plt.xticks(rotation=90)
        figName = './intervalCountPerLength/intervalCountPerLength_'
        figName += label[i]
        pltTitle = label[i]
        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(label[i],' OK')



# 横軸Interbal、縦軸Countの組織ごと(All除く)の順位の折れ線グラフ
# 17はデータ数0のため除く
# rootは除く
if False:
    print('x axis:Interval, y axis: Rank, for all tissues')
    new_dir_path_recursive = './intervalRank'
    os.makedirs(new_dir_path_recursive, exist_ok=True)
    # 各Tissue
    for i in range(tissueCount):
        left = list(range(humanPathSize-2))
        rank = []
        for j in range(humanPathSize-2):
            keep = []
            for k in range(tissueCount):
                keep.append(ratioTissues[k][j])
            keep.sort(reverse=True)
            for k in range(tissueCount):
                if keep[k] == ratioTissues[i][j]:
                    rank.append(tissueCount - k)
                    break
        print(rank)
        
        plt.figure(figsize=(16.0,12.0))
        plt.rcParams["font.size"] = 18
        plt.plot(left, rank, label=str(label[i]))
        plt.xlabel('Concestor')
        plt.ylabel('rank')
        plt.xticks(rotation=90)

        ax = plt.gca()
        ax.axes.yaxis.set_ticks([])


        figName = './intervalRank/intervalRank_'
        figName += label[i]
        pltTitle = label[i]
        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(label[i],' OK')




# 折れ線グラフ用x軸ラベル
if True:
    print('x axis:Interval, y axis: Count, for all tissues 2')
    new_dir_path_recursive = './xLabel'
    os.makedirs(new_dir_path_recursive, exist_ok=True)
    # All
    left = list(range(humanPathSize-1))
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    countAllArgMax2 = []
    for i in range(len(countAllArgMax)-1):
        countAllArgMax2.append(countAllArgMax[i])
    concestorLabel2 = []
    for i in range(len(concestorLabel)-1):
        concestorLabel2.append(concestorLabel[i])
    plt.bar(left, countAllArgMax2, color=colorlist, tick_label=concestorLabel2, align="center")
    plt.xlabel('Concestor')
    plt.ylabel('counts')
    plt.xticks(rotation=90)

    #ax = plt.gca()
    #ax.axes.yaxis.set_ticks([])
    figName = './xLabel/xLabel18'
    pltTitle = 'all'
    plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('all OK')





# bar_ratio
if False:
    if modeRootCut == 0:
        new_dir_path_recursive = './bar_ratio_new'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    else:
        new_dir_path_recursive = './bar_ratio_root_cut_new'
        os.makedirs(new_dir_path_recursive, exist_ok=True)

    print('make ratio figures')
    for ite in range(humanPathSize):
        left = list(range(tissueCount))
        height = []
        labelTitle = []
        for j in range(tissueCount):
            height.append(ratioTissues[j][ite])
            labelString = label[j] + '\n' + str(countTissues[j][ite])
            labelTitle.append(labelString)

        plt.figure(figsize=(16.0,12.0))
        plt.bar(left, height, color=colorlist, tick_label=labelTitle, align="center")
        plt.xticks(rotation=90)
        figName = ''
        if modeRootCut == 0:
            figName = './bar_ratio_new/bar_ratio_'
        else:
            figName = './bar_ratio_root_cut_new/bar_ratio_root_cut_'
        if ite <= 9:
            figName += '0'
        figName += str(ite)
        figName += '_' + concestor[ite]

        pltTitle = ''
        if ite <= 9:
            pltTitle += '0'
        pltTitle += str(ite)
        pltTitle += '_' + concestor[ite]
        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(concestor[ite],' OK')




'''
    横軸に器官、縦軸に個数のグラフをconcestor interval毎に描く
'''
# bar_count
if False:
    if modeRootCut == 0:
        new_dir_path_recursive = './bar_count_new'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    else:
        new_dir_path_recursive = './bar_count_root_cut_new'
        os.makedirs(new_dir_path_recursive, exist_ok=True)

    print('make bar figures')
    for ite in range(humanPathSize):
        left = list(range(tissueCount))
        height = []
        for j in range(tissueCount):
            height.append(countTissues[j][ite])
        plt.figure(figsize=(16.0,12.0))
        plt.bar(left, height, color=colorlist, tick_label=label, align="center")
        plt.xticks(rotation=90)
        figName = ''
        if modeRootCut == 0:
            figName = './bar_count_new/bar_count_'
        else:
            figName = './bar_count_root_cut_new/bar_count_root_cut_'
        if ite <= 9:
            figName += '0'
        figName += str(ite)
        figName += '_' + concestor[ite]
        plt.title(concestor[ite])
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(concestor[ite],' OK')


# Top10のみのratioグラフ
# bar_ratio_top10
if False:
    sortedSum = sorted(sumTissues, reverse=True)

    if modeRootCut == 0:
        new_dir_path_recursive = './bar_ratio_top10'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    else:
        new_dir_path_recursive = './bar_ratio_root_cut_top10'
        os.makedirs(new_dir_path_recursive, exist_ok=True)

    print('make ratio top10 figures')
    for ite in range(humanPathSize):
        left10 = list(range(10))
        height10 = []
        label10 = []

        for j in range(tissueCount):
            if sumTissues[j] >= sortedSum[9]:
                height10.append(ratioTissues[j][ite])
                labelString = label[j]
                #labelString = label[j] + '\n' + str(countTissues[j][ite])
                label10.append(labelString)

        plt.figure(figsize=(16.0,12.0))
        plt.rcParams["font.size"] = 18
        plt.bar(left10, height10, color=colorlist, tick_label=label10, align="center")
        plt.xlabel('tissue')
        plt.ylabel('rate')
        plt.xticks(rotation=30)
        figName = ''
        if modeRootCut == 0:
            figName = './bar_ratio_top10/bar_ratio_top10_'
        else:
            figName = './bar_ratio_root_cut_top10/bar_ratio_root_cut_top10_'
        if ite <= 9:
            figName += '0'
        figName += str(ite)
        figName += '_' + concestor[ite]

        pltTitle = concestor[ite]
        '''
        if ite <= 9:
            pltTitle += '0'
        pltTitle += str(ite)
        pltTitle += '_' + concestor[ite]
        pltTitle += '_' + str(lengthOfHumanPathEdges[ite])
        '''

        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(concestor[ite],' OK')

# 横軸Human path, 縦軸ratio
# oresen_top10
if False:
    sortedSum = sorted(sumTissues, reverse=True)

    if modeRootCut == 0:
        new_dir_path_recursive = './oresen_top10'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    
    print('oresen top10 figures')
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    plt.xticks(rotation=30)

    for i in range(tissueCount):
        if sumTissues[i] >= sortedSum[9]:
            left = []
            height = []
            label10 = []
            for j in range(humanPathSize):
                if j == 18:
                    continue
                left.append(j)
                label10.append(concestor[j])
                height.append(ratioTissues[i][j])
            plt.plot(left, height, label=str(label[i]))
    figName = './oresen_top10/oresen_top10'
    pltTitle = 'Top_10'
    plt.legend(loc=0)
    #plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('OK')

# oresen_top10 その2(枝長で割った方)
if True:
    sortedSum = sorted(sumTissues, reverse=True)

    if modeRootCut == 0:
        new_dir_path_recursive = './oresen_top10'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    
    print('oresen top10 figures')
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    plt.xticks(rotation=30)

    for i in range(tissueCount):
        if sumTissues[i] >= sortedSum[9]:
            left = []
            height = []
            label10 = []
            for j in range(humanPathSize):
                if j == 18:
                    continue
                left.append(j)
                label10.append(concestor[j])
                height.append(ratioTissues[i][j] / lengthOfHumanPathEdges[j])
            plt.plot(left, height, label=str(label[i]))
    figName = './oresen_top10/oresen_top10_2'
    pltTitle = 'Top_10'
    plt.legend(loc=0)
    #plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('OK')




# oresen_top20
if False:
    sortedSum = sorted(sumTissues, reverse=True)

    if modeRootCut == 0:
        new_dir_path_recursive = './oresen_top20'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    
    print('oresen top20 figures')
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    plt.xticks(rotation=30)

    for i in range(tissueCount):
        if sumTissues[i] >= sortedSum[19]:
            left = []
            height = []
            label20 = []
            for j in range(humanPathSize):
                if j == 18:
                    continue
                left.append(j)
                label20.append(concestor[j])
                height.append(ratioTissues[i][j])
            plt.plot(left, height, label=str(label[i]))
    figName = './oresen_top10/oresen_top20'
    pltTitle = 'Top_20'
    plt.legend(loc=0)
    plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('OK')

# oresen_top20　その2(枝長で割った方)
if True:
    sortedSum = sorted(sumTissues, reverse=True)

    if modeRootCut == 0:
        new_dir_path_recursive = './oresen_top10'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    
    print('oresen top20 figures')
    plt.figure(figsize=(16.0,12.0))
    plt.rcParams["font.size"] = 18
    plt.xticks(rotation=30)

    for i in range(tissueCount):
        if sumTissues[i] >= sortedSum[19]:
            left = []
            height = []
            label20 = []
            for j in range(humanPathSize):
                if j == 18:
                    continue
                left.append(j)
                label20.append(concestor[j])
                height.append(ratioTissues[i][j] / lengthOfHumanPathEdges[j])
            plt.plot(left, height, label=str(label[i]))
    figName = './oresen_top10/oresen_top20_2'
    pltTitle = 'Top_20'
    plt.legend(loc=0)
    #plt.title(pltTitle)
    plt.tight_layout()
    plt.savefig(figName)
    plt.close()
    print('OK')


    

# Top10のbarグラフ
# bar_count_top10
if False:
    sortedSum = sorted(sumTissues, reverse=True)

    if modeRootCut == 0:
        new_dir_path_recursive = './bar_count_top10'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    else:
        new_dir_path_recursive = './bar_count_root_cut_top10'
        os.makedirs(new_dir_path_recursive, exist_ok=True)

    print('make bar top10 figures')
    for ite in range(humanPathSize):
        left10 = list(range(10))
        height10 = []
        label10 = []
        for j in range(tissueCount):
            if sumTissues[j] >= sortedSum[9]:
                height10.append(countTissues[j][ite])
                label10.append(label[j])
                
        plt.figure(figsize=(16.0,12.0))
        plt.bar(left10, height10, color=colorlist, tick_label=label10, align="center")
        plt.xticks(rotation=90)
        figName = ''
        if modeRootCut == 0:
            figName = './bar_count_top10/bar_count_top10_'
        else:
            figName = './bar_count_root_cut_top10/bar_count_root_cut_top10_'
        if ite <= 9:
            figName += '0'
        figName += str(ite)
        figName += '_' + concestor[ite]
        plt.title(concestor[ite])
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(concestor[ite],' OK')


# Top20のみのratioグラフ
#bar_ratio_top20
if False:
    sortedSum = sorted(sumTissues, reverse=True)

    if modeRootCut == 0:
        new_dir_path_recursive = './bar_ratio_top20'
        os.makedirs(new_dir_path_recursive, exist_ok=True)
    else:
        new_dir_path_recursive = './bar_ratio_root_cut_top20'
        os.makedirs(new_dir_path_recursive, exist_ok=True)

    print('make ratio top20 figures')
    for ite in range(humanPathSize):
        left20 = list(range(20))
        height20 = []
        label20 = []

        for j in range(tissueCount):
            if sumTissues[j] >= sortedSum[19]:
                height20.append(ratioTissues[j][ite])
                labelString = label[j] + '\n' + str(countTissues[j][ite])
                label20.append(labelString)

        plt.figure(figsize=(32.0,24.0))
        plt.rcParams["font.size"] = 18
        plt.bar(left20, height20, color=colorlist, tick_label=label20, align="center")
        plt.xticks(rotation=90)
        figName = ''
        if modeRootCut == 0:
            figName = './bar_ratio_top20/bar_ratio_top20_'
        else:
            figName = './bar_ratio_root_cut_top20/bar_ratio_root_cut_top20_'
        if ite <= 9:
            figName += '0'
        figName += str(ite)
        figName += '_' + concestor[ite]
        

        pltTitle = ''
        if ite <= 9:
            pltTitle += '0'
        pltTitle += str(ite)
        pltTitle += '_' + concestor[ite]
        pltTitle += '_' + str(lengthOfHumanPathEdges[ite])

        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print(concestor[ite],' OK')



# countPerLength.png
if False:
    print('make sum per edgeLength figure')

    for ite in range(1):
        left = list(range(humanPathSize))
        height = []
        labelConcestor = []

        for i in range(humanPathSize):
            if i == humanPathSize - 1:
                height.append(0.0)
            else:
                height.append(countAllArgMax[i] / lengthOfHumanPathEdges[i])
            labelConcestor.append(concestor[i])
            
        tempValue = height[0]
        for i in range(humanPathSize):
            height[i] = height[i] / tempValue

        plt.figure(figsize=(32.0,16.0))
        plt.rcParams["font.size"] = 12
        plt.bar(left, height, color=colorlist, tick_label=labelConcestor, align="center")
        plt.xticks(rotation=90)
        figName = ''
        if modeRootCut == 0:
            figName = './countPerLength'
        else:
            figName = './countPerLength'
        

        pltTitle = 'countPerLength'

        plt.title(pltTitle)
        plt.tight_layout()
        plt.savefig(figName)
        plt.close()
        print('OK')





elapsed_time = time.time() - startTime
print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")
