#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import numpy as np

os.chdir(r'/gss1/home/gaohao/Oryza_sativa_methylation_picture/A_1.methylation')#修改工作目录

inputfile1 = open("all.gff3")
inputfile2 = open("A_1.sortBYhao.CpG_report.txt") #修改要输入的文件名

def get_area_me_mean(all_area_num,gb_area_num,gud_area_length):#总共要划分为多少个区段（偶整数），基因内部划分多少个区段（偶整数），基因上下游取多少长度（bp）,基因上下游取的长度除以上下游的区段数目要是个整数
    ID = []
    gene_start = []
    gene_end = []
    strand = []
    x = 0
    for line in inputfile1:
        temp = line.split("\t")
        if (len(temp) <=3)&(x == 0):
            x = 2
            continue
        if (len(temp) <=3)&(x == 2):
            break
    
        if temp[2] == "gene":
            ID.append(int(temp[0][3:5]))
            gene_start.append(int(temp[3]))
            gene_end.append(int(temp[4]))
            strand.append(temp[6])

    for x in range(len(strand)):
        if strand[x] == "+":
            strand[x] = 1
        elif strand[x] == "-":
            strand[x] = -1

    data = {"chromosome":ID,"gene_start":gene_start,"gene_end":gene_end,"strand":strand} #构建DataFrame
    Tigr7_genes_start_end = pd.DataFrame(data)

#Tigr7_genes_start_end.to_csv(path_or_buf = r"/gss1/home/gaohao/Oryza_sativa_methylation_picture/A_1.methylation/Tigr7_genes_start_end.csv")

    _TSS_TES_ = pd.DataFrame()
    for x in range(81): #基因上游分20段，基因内部分40段，基因下游分20段
        gene_distance = list(np.array(Tigr7_genes_start_end["gene_end"]) - np.array(Tigr7_genes_start_end["gene_start"]))
        gene_bin = map(lambda y : y/gb_area_num,gene_distance)
        gup_area_num = (all_area_num - gb_area_num)/2
        gdown_area_num = all_area_num - gup_area_num
        if x <= gup_area_num:
            _TSS_TES_[x] = Tigr7_genes_start_end["gene_start"].map(lambda y : y - (gud_area_length/gup_area_num)*(gup_area_num-x))
        if (x > gup_area_num)&(x <= gdown_area_num):
            _TSS_TES_[x] = list(map(lambda y: int(y),list(np.array(Tigr7_genes_start_end["gene_start"]) + np.array(list(map(lambda y : y*(x-gup_area_num),gene_bin))))))
        #解释一下：先将应该加的基因区间数值算出来，放在列表中，再转化为array，gene_start也转化为array，相加之后再转化为列表，并且使列表中的每一个数都转化为整数
        if x > gdown_area_num:
            _TSS_TES_[x] = Tigr7_genes_start_end["gene_end"].map(lambda y : y + (gud_area_length/gup_area_num)*(x-gdown_area_num))

    _TSS_TES_["chromosome"] = Tigr7_genes_start_end["chromosome"]

    _TSS_TES_.to_csv(path_or_buf  = r"/gss1/home/gaohao/Oryza_sativa_methylation_picture/A_1.methylation/_TSS_TES_.csv")

    chromosome = []
    position = []
    #strand = []
    count_methylated = []
    count_unmethylated = []
    #C_context = []
    #trinucleotide_context = []
    for line in inputfile2:
        temp = line.split("\t")
        if len(temp[0]) >= 6:
            continue
        else:
            chromosome.append(temp[0])
            position.append(int(temp[1]))
            #strand.append(temp[2])
            count_methylated.append(int(temp[3]))
            count_unmethylated.append(int(temp[4]))
            #C_context.append(temp[5])
            #trinucleotide_context.append(str(temp[6]).strip('\n'))_1.sortBYhao.CpG_report.txt

    data1 = {"chromosome":chromosome,"position":position,"count_methylated":count_methylated,"count_unmethylated":count_unmethylated}

    A_1_CpG = pd.DataFrame(data1)

    A_1_CpG.to_csv(path_or_buf  = r"/gss1/home/gaohao/Oryza_sativa_methylation_picture/A_1.methylation/CpG.csv")

    inputfile3 = open("CpG.csv")#每列依次为：行号,chromosome,position,count_methylated,count_unmethylated
    inputfile4 = open("_TSS_TES_.csv") #每列依次为：行号，0-80也就是各区段的的端点，chromosome
    inputfile5 = open("_TSS_TES_.csv")
    testfile = open("test.txt","w")

#将A_1_CpG.csv里的数据转化为{染色体：{胞嘧啶位置：[甲基化的数目，未甲基化的数目]}}这样的格式
    gdict = {}
    for line in inputfile3:
        temp = line.split(',')
        if temp[1] == "chromosome":
            continue
        else:
            chr = int(temp[1][3:5])
            c_pos = int(temp[2])
            me_num = [int(temp[3]),int(temp[4])]
            if chr not in gdict:
                gdict[chr] = {}
                gdict[chr][c_pos]= me_num
            else:
                if c_pos not in gdict[chr]:
                    gdict[chr][c_pos] = me_num
                else:
                    gdict[chr][c_pos].append(me_num)
    testfile.write("gdict finished\n")

#这样就可以通过染色体号以及该胞嘧啶的位置以及它的甲基化状态



#将_TSS_TES_.csv文件中的数据转化为{染色体：{"区段号-基因号"：[区段前端位置，区段后端位置]}}
    chr_gene_num = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0} #把每条染色体上的基因数目记录一下，下面遍历要用
    col = all_area_num + 2
    for line in inputfile4:
        temp = line.split(',')
        if temp[col] == "chromosome\n":
            continue
        else:
            chr = int(temp[col])
            chr_gene_num[chr] = chr_gene_num[chr] + 1

    def chr_gene_num_acc(x):
        y = 0
        for x in range(x):
            y = y + chr_gene_num[x+1]
        return y


    hdict = {}
    row = 0
    for line in inputfile5:
        temp = line.split(',')
        if temp[col] == "chromosome\n":
            continue
        else:
            for area_point in range(all_area_num):
                chr  = int(temp[col])
                pos_name = (area_point + 1,row +1 - chr_gene_num_acc(chr-1))
                pos = [float(temp[area_point+1]),float(temp[area_point + 2])]
                if chr not in hdict:
                    hdict[chr] = {}
                    hdict[chr][pos_name] = pos
                else:
                    if pos_name not in hdict[chr]:
                        hdict[chr][pos_name] = pos
                    else:
                        hdict[chr][pos_name].append(pos)
            row = row + 1
    testfile.write("hdict finished\n")
#这样可以通过染色体号以及区段号-基因号(每个染色体都要从1开始编号)找到某个区段的具体位置

    area_id_me_level = {}
    for x in range(1,all_area_num+1):
        area_id_me_level[x] = []

    area_id = 1                                  #将所有胞嘧啶分配到一号染色体一号区域上，并算出这个区域每个区段的甲基化水平，存到第一个列表中（测试时候的想法）
    for area_id in range(1,all_area_num+1):

        for chr_id in range(1,13):

            me = 0
            all_me = 0
            gene_id = 1

            for c_pos in gdict[chr_id]:
    
                 if gene_id == chr_gene_num[chr_id] + 1:
                     break
    
                 if c_pos < hdict[chr_id][(area_id,gene_id)][0]:
                     continue
                 elif c_pos <= hdict[chr_id][(area_id,gene_id)][1]:
                     me = me + gdict[chr_id][c_pos][0]
                     all_me  = all_me + gdict[chr_id][c_pos][0] + gdict[chr_id][c_pos][1]
                 elif (c_pos > hdict[chr_id][(area_id,gene_id)][1])&(me > 0):
                     this_area_me_level = me/all_me
                     area_id_me_level[area_id].append(this_area_me_level)
                     me = 0
                     all_me = 0
                     gene_id = gene_id + 1
                 else :
                     this_area_me_level = 0
                     area_id_me_level[area_id].append(this_area_me_level)
                     all_me = 0
                     gene_id = gene_id + 1

            for append_zero in range(len(area_id_me_level[area_id]),chr_gene_num_acc(chr_id)):
                     area_id_me_level[area_id].append(0)

              
    area_id_me_level_reverse = pd.DataFrame(area_id_me_level)
    #area_id_me_level_reverse.to_csv(path_or_buf = r"/gss1/home/gaohao/Oryza_sativa_methylation_picture/A_1.methylation/area_id_me_level_noreverse_new.csv")


    me_mean_area = []
    all_gene_num = chr_gene_num[1]+chr_gene_num[2]+chr_gene_num[3]+chr_gene_num[4]+chr_gene_num[5]+chr_gene_num[6]+chr_gene_num[7]+chr_gene_num[8]+chr_gene_num[9]+chr_gene_num[10]+chr_gene_num[11]+chr_gene_num[12]
    for x in range(all_gene_num):
        if Tigr7_genes_start_end.loc[x,"strand"] < 0:
            area_id_me_level_reverse.loc[x] = list(area_id_me_level_reverse.loc[x])[::-1]

    #data2 = area_id_me_level_reverse
    #data2.to_csv(path_or_buf = r"/gss1/home/gaohao/Oryza_sativa_methylation_picture/A_1.methylation/area_id_me_level_reverse_new.csv")

    for area_id in range(1,all_area_num + 1):
        me_mean = np.mean(list(area_id_me_level_reverse[area_id]))
        me_mean_area.append(me_mean)
    
    print(me_mean_area)
#把负链基因的区段甲基化程度颠倒之后，求各区段的平均值,输出在log文件里
    
    for x in range(all_area_num):
        testfile.write(str(len(area_id_me_level[x+1]))+"\n")
    
    inputfile1.close()
    inputfile2.close()
    inputfile3.close()
    inputfile4.close()
    inputfile5.close()
    testfile.close()

get_area_me_mean(80,40,1000)
    

