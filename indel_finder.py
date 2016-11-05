#!/usr/bin/env python
#fn: bindel.py
import os, getopt,sys 
import gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def usage():
    print("\nUsage: python -r ref.fa(5' to 3') -t target.fa(5' to 3' with PAM) -1 1_S1_L001_R1_001.fastq(or fastq.gz) -2 2_S1_L001_R2_001.fastq -o out_dir\n")

out_dir = "./"
opts, args = getopt.getopt(sys.argv[1:], "hr:t:1:2:o:")

if len(opts) == 0:
    usage()
    sys.exit()
for op, value in opts:
    if op == "-r":
        ref_file = value
    elif op == "-t":
        target_file = value
    elif op == "-1":
        read1 = value
    elif op == "-2":
        read2 = value
    elif op == "-o":
        if out_dir[-1] != '/':
            out_dir = value+'/'
        else:
            out_dir = value
if not os.path.exists(ref_file):
    print("The path [{}] do not exit!".format(ref_file))
    sys.exit()
if not os.path.exists(target_file):
    print("The path [{}] do not exit!".format(target_file))
    sys.exit()
if not os.path.exists(read1):
    print("The path [{}] do not exit!".format(read1))
    sys.exit()
if not os.path.exists(read2):
    print("The path [{}] do not exit!".format(read2))
    sys.exit()

# judge fastq.gz or fastq
def isFastq(input_file):
    #gzext = ("fq.gz", "fastq.gz")
    #for ext in gzext:
    if input_file.endswith("q.gz"):
        return gzip.open(input_file, 'rb')
    else:
        return open(input_file, 'r')


# get indel from '*.clear.out' file        
def getindel(clear_file, start_val, rlt_dir, fn):
    ra = 36 # half range 
    rang = [start_val - ra, start_val + ra, start_val] # star_val is the third base upstream of PAM
    FILE = open(clear_file, 'r')
    #OUTFILE = open(rlt_dir+fn[:-7]+'.bigap.out', 'w')
    #OUT2 = open(rlt_dir+fn[:-7]+'.loc2.out','w') # hit two location and gap < 50
    SUMMARY = open(rlt_dir+fn[:-3]+'.info','w')
    H1 = open(rlt_dir+fn[:-3]+'.noindel.out','w')
    #H2 = open(rlt_dir+fn[:-3]+'.h2.out','w')
    DELET = open(rlt_dir+fn[:-3]+'.delet.out','w')
    INSER = open(rlt_dir+fn[:-3]+'.inser.out','w')
    #LOC_3 = open('3_loc.txt','w')
    #LOC_3.write('>'+fn+'\n')
    LOC_DEL = open(rlt_dir+fn[:-3]+'.pos_del.txt','w')
    LOC_INS = open(rlt_dir+fn[:-3]+'.pos_ins.txt','w')
    LOC_IDL = open(rlt_dir+fn[:-3]+'.pos_idl.txt','w')
    NUM_DEL = open(rlt_dir+fn[:-3]+'.size_del.txt','w')
    NUM_INS = open(rlt_dir+fn[:-3]+'.size_ins.txt','w')
    NUM_IDL = open(rlt_dir+fn[:-3]+'.size_idl.txt','w')

    LOC_DEL.write(fn[:-3]+'\n')
    LOC_INS.write(fn[:-3]+'\n')
    LOC_IDL.write(fn[:-3]+'\n')
    NUM_DEL.write(fn[:-3]+'\n')
    NUM_INS.write(fn[:-3]+'\n')
    NUM_IDL.write(fn[:-3]+'\n')


    delet_ls=np.zeros(60)
    inser_ls=np.zeros(60)
    loc_inser = np.zeros(2*ra)
    loc_delet = np.zeros(2*ra)
    rsum = 0

    l_4 = 0
    l_3 = 0
    l_2 = 0
    l_2_r = 0
    l_1 = 0
    bigap = 0
    r = 0
    lnum = 0
    for line in FILE:
        line = line.strip('\n')
        ls = line.split('\t')
        loc_query = ls[3].split('|')
        loc_sbjct = ls[4].split('|')
        if len(loc_sbjct) == 2:
            l_2 += 1
            loc_1 = loc_sbjct[0].split(';')#[0].split(':')[0]
            loc_2 = loc_sbjct[1].split(';')#[0].split(':')[1]
            loc_1_1 = loc_1[0].split(':')
            loc_2_1 = loc_2[0].split(':')
            # -
            if int(loc_1_1[0]) - int(loc_1_1[1]) > 0 and int(loc_2_1[0]) - int(loc_1_1[1]) > 0:
                #gap > 50
                if int(loc_1[-1].split(':')[1]) - int(loc_2[0].split(':')[0]) > 50:
                    bigap += 1
                    #OUTFILE.write(line+'\t-\n')
                #elif int(loc_1[-1].split(':')[1]) - int(loc_2[0].split(':')[0]) < 50:
                    #OUT2.write(line+'\n')
                    #print(line)
            # +
            elif int(loc_1_1[0]) - int(loc_1_1[1]) < 0 and int(loc_2_1[0]) - int(loc_1_1[1]) < 0:
                if int(loc_1[0].split(':')[0]) - int(loc_2[-1].split(':')[1]) > 50:
                    bigap += 1
                    #OUTFILE.write(line+'\t+\n')
                #elif int(loc_1[0].split(':')[0]) - int(loc_2[-1].split(':')[1]) < 50:
                    #OUT2.write(line+'\n')
                    #print(line)
            # +:- or -:+
            else:
                l_2_r += 1

        # hit 1 location
        if len(loc_sbjct) == 1:
            l_1 += 1
            query = ls[1]
            sbjct = ls[2]
            sbjct_ls = loc_sbjct[0].split(':')
            start = int(sbjct_ls[0])
            end = int(sbjct_ls[-1])
            #+
            if start - end < 0:
                if start < rang[0] + 30 and end > rang[1] - 20:
                    #if '-' in ls[1] or '-' in ls[2]:
                        #H1.write(line+'\n')
                    rsum += 1
                    count = 0
                    delet = 0
                    inser = 0
                    #for delet
                    for i in query:
                        count += 1
                        if i == '-' and query[count] != '-':
                            delet += 1
                            if count +start - delet -1 <= rang[2] + 5 and count +start - 1 >= rang[2] - 5: #- 30
                                #if delet > 49:
                                    #print(line)
                                loc_delet[count+start-delet-rang[0]] += 1
                                delet_ls[delet] += 1
                                DELET.write(line+'\n')
                                break
                            else:
                                if '-' not in ls[2]:
                                    H1.write(line+'\n')
                                delet = 0
                                #continue
                        if i == '-' and query[count] == '-':
                            delet += 1
                    #for inser
                    count = 0
                    for i in sbjct:
                        count += 1
                        if i == '-' and sbjct[count] != '-':
                            inser += 1
                            if count + start - inser -1 < rang[2] + 5 and count + start -1 >= rang[2] - 5:
                                loc_inser[count+start-inser-rang[0]] += 1
                                inser_ls[inser] += 1
                                INSER.write(line+'\n')
                                break
                            else:
                                H1.write(line+'\n')
                                inser = 0
                                #continue
                        if i == '-' and sbjct[count] == '-':
                            inser += 1
            #-
                #else:
                    #H2.write(line+'\n')
            if start - end > 0:
                query = query[::-1]
                a = start
                start = end
                end = a
                if start < rang[0] + 30 and end > rang[1] - 20:
                    #if '-' in ls[1] or '-' in ls[2]:
                        #H1.write(line+'\n')
                    rsum += 1
                    count = 0
                    delet = 0
                    inser = 0
                    #for delet
                    for i in query:
                        count += 1
                        if i == '-' and query[count] != '-':
                            delet += 1
                            if count +start - delet -1 < rang[2] + 5 and count +start - 1 > rang[2] - 5:#- 30:
                                loc_delet[count+start-delet-rang[0]] += 1
                                delet_ls[delet] += 1
                                DELET.write(line+'\n')
                                break
                            else:
                                if '-' not in ls[2]:
                                    H1.write(line+'\n')
                                delet = 0
                                #continue
                        if i == '-' and query[count] == '-':
                            delet += 1
                    #for inser
                    count = 0
                    for i in sbjct:
                        count += 1
                        if i == '-' and sbjct[count] != '-':
                            inser += 1
                            if count + start - inser - 1 < rang[2] + 5 and count + start - 1 > rang[2] - 5:
                                loc_inser[count+start-inser-rang[0]] += 1
                                inser_ls[inser] += 1
                                INSER.write(line+'\n')
                                break
                            else:
                                H1.write(line+'\n')
                                inser = 0
                                #continue
                        if i == '-' and sbjct[count] == '-':
                            inser += 1
            #print(line)
                #else:
                    #H2.write(line+'\n')

        if len(loc_sbjct) == 3:
            #print(line)
            #LOC_3.write(line+'\n')
            l_3 += 1
        if len(loc_sbjct) == 4:
            #print(line)
            l_4 += 1
    #print(delet_ls)
    #print(inser_ls)

    print('l_1: '+str(l_1))
    print('l_2: '+str(l_2))
    #print('l_3: '+str(l_3))
    #print('l_4: '+str(l_4))
    print('for l_2:\nr2_+/-: '+str(l_2_r))
    print('bigap: '+str(bigap))
    
    indel = sum(delet_ls) + sum(inser_ls)
    indel_ls = delet_ls + inser_ls
    loc_indel = loc_delet + loc_inser
    ratio = indel*1.0/rsum
    print('delet: '+str(sum(delet_ls)))
    print('inser: '+str(sum(inser_ls)))
    print('covered reads: '+str(rsum))
    print('indel: '+str(indel))
    print('ratio: '+str(ratio))
    inser_str = '1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\t25\n'
    delet_str = '1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\t25\n'
#initialize lists for plotting
    num_ins_x = []
    num_ins_y = []
    num_del_x = []
    num_del_y = []
    num_idl_x = []
    num_idl_y = []

    loc_ins_x = []
    loc_ins_y = []
    loc_del_x = []
    loc_del_y = []
    loc_idl_x = []
    loc_idl_y = []

    for i in range(len(inser_ls)-1):
        NUM_INS.write(str(i+1)+'\t'+str(inser_ls[i+1]*100/sum(inser_ls))+'\n')
        num_ins_x.append(i+1)
        num_ins_y.append(inser_ls[i+1]*100/sum(inser_ls))
    for i in range(len(delet_ls)-1):
        NUM_DEL.write(str(i+1)+'\t'+str(delet_ls[i+1]*100/sum(delet_ls))+'\n')
        num_del_x.append(i+1)
        num_del_y.append(delet_ls[i+1]*100/sum(delet_ls))
    for i in range(len(indel_ls)-1):
        NUM_IDL.write(str(i+1)+'\t'+str(indel_ls[i+1]*100/sum(indel_ls))+'\n')
        num_idl_x.append((i+1))
        num_idl_y.append(indel_ls[i+1]*100/sum(indel_ls))
    for i in range(len(loc_delet)-1):
        LOC_DEL.write(str(i+1-ra)+'\t'+str(loc_delet[i+1]*100/sum(loc_delet))+'\n')
        loc_del_x.append(i+1-ra)
        loc_del_y.append(loc_delet[i+1]*100/sum(loc_delet))
    for i in range(len(loc_inser)-1):
        LOC_INS.write(str(i+1-ra)+'\t'+str(loc_inser[i+1]*100/sum(loc_inser))+'\n')
        loc_ins_x.append((i+1-ra))
        loc_ins_y.append(loc_inser[i+1]*100/sum(loc_inser))
    for i in range(len(loc_indel)-1):
        LOC_IDL.write(str(i+1-ra)+'\t'+str(loc_indel[i+1]*100/sum(loc_indel))+'\n')
        loc_idl_x.append((i+1-ra))
        loc_idl_y.append(loc_indel[i+1]*100/sum(loc_indel))
    getFigure(num_ins_x, num_ins_y, fn[:-3]+' Insertion Size', 'Size', 'Ratio(%)', fn[:-3]+'_size_ins.jpg',fn)
    getFigure(num_del_x, num_del_y, fn[:-3]+' Deletion Size', 'Size', 'Ratio(%)', fn[:-3]+'_size_del.jpg', fn)
    getFigure(num_idl_x, num_idl_y, fn[:-3]+' Indel Size', 'Size', 'Ratio(%)', fn[:-3]+'_size_idl.jpg',fn)
    getFigure(loc_ins_x, loc_ins_y, fn[:-3]+' Insertion Position', 'Position', 'Ratio(%)', fn[:-3]+'_pos_ins.jpg',fn)
    getFigure(loc_del_x, loc_del_y, fn[:-3]+' Deletion Position', 'Position', 'Ratio(%)', fn[:-3]+'_pos_del.jpg',fn)
    getFigure(loc_idl_x, loc_idl_y, fn[:-3]+' Indel Position', 'Position', 'Ratio(%)', fn[:-3]+'_pos_idl.jpg',fn)

    for i in range(len(inser_ls)-1):
        if i > 50:
            break
        inser_str = inser_str+str(inser_ls[i+1])+'\t'
        if i == 24:
            inser_str = inser_str+'\n26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t41\t42\t43\t44\t45\t46\t47\t48\t49\t50\n'
    for i in range(len(delet_ls)-1):
        if i > 50:
            break
        delet_str = delet_str+str(delet_ls[i+1])+'\t'
        if i == 24:
            delet_str = delet_str+'\n26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t41\t42\t43\t44\t45\t46\t47\t48\t49\t50\n'
    SUMMARY.write('hit_1: '+str(l_1)+'\n')
    SUMMARY.write('hit_2: '+str(l_2)+'\n')
    #SUMMARY.write('hit_3: '+str(l_3)+'\n')
    #SUMMARY.write('hit_4: '+str(l_4)+'\n')
    SUMMARY.write('for hit_2:\n\thit2_+/-: '+str(l_2_r)+'\n')
    SUMMARY.write('\tbigap(gap > 50): '+str(bigap)+'\n')
    SUMMARY.write('Covered reads: '+str(rsum)+'\n')
    SUMMARY.write('Indel: '+str(indel)+'\n')
    SUMMARY.write('Ratio(Indel/Covered reads): '+str(ratio)+'\n')
    SUMMARY.write('Insertion: '+str(sum(inser_ls))+'\n')
    SUMMARY.write('Deletion: '+str(sum(delet_ls))+'\n')
    SUMMARY.write('\n')
    SUMMARY.write('Insertion: '+str(sum(inser_ls))+'\n')
    SUMMARY.write(inser_str+'\n')
    SUMMARY.write('Deletion: '+str(sum(delet_ls))+'\n')
    SUMMARY.write(delet_str+'\n')
    #SUMMARY.write('---------------------------------------------------------------------------\n')
    SUMMARY.close()



# from blastn output result to another format that is easy to be processed
def clear(fl,fn):
    rlt_dir = out_dir + fn[:-3]+'_result/'
    if not os.path.exists(rlt_dir):
        os.makedirs(rlt_dir)
    FILE = open(fl,'r')
    #print(fn[-5])
    if int(fn[-1]) == 1:
        OUTFILE = open(fl[:-7]+'.clear.out','w')
# output blastn result to a format that is easy to process 
    if int(fn[-1]) == 2:
        OUTFILE = open(fl[:-7]+'.clear.out','a') 
# output basic information of the blastn result
    INFOFILE = open(rlt_dir+fn[:-3]+'.info','w')
    INFOFILE.write(fn[:-3]+'\n')
    lnum = 0
    mark = 0
    query = ''
    sbjct = ''
    loc_q = ''
    loc_s = ''
    rn = ''
    nohit = 0
    hit = 0
    for line in FILE:
        lnum += 1
        if lnum < 15:
            continue
        line = line.strip('\n')
        ls = line.split('\t')
        if len(line) == 0:
            continue
        if '*****' in line:
            nohit += 1
        if line[:6] == 'Query=':
            #print(line)
            newrn = line[7:] 
        if line[:8] == ' Strand=':
            query = query + ';'
            sbjct = sbjct + ';'
            loc_q = loc_q[:-1] + '|'
            loc_s = loc_s[:-1] + '|'

        if line[0] == '>':
            hit += 1
            if hit == 1:
                rn = newrn
                continue
            OUTFILE.write(rn+'\t'+query[1:]+'\t'+sbjct[1:]+'\t'+loc_q[1:-1]+'\t'+loc_s[1:-1]+'\n')
            rn = newrn
            query = ''
            sbjct = ''
            loc_q = ''
            loc_s = ''
            
            #OUTFILE.write(line)
        elif line[:6] == 'Query ':
            ls = line.split(' ')
            #print(ls[6])
            if len(ls) == 8:
                query = query+ls[5]
                loc_q = loc_q+ls[2]+':'+ls[7]+';'
            if len(ls) == 7:
                query = query+ls[4]
                loc_q = loc_q+ls[2]+':'+ls[6]+';'
            if len(ls) == 9:
                query = query+ls[6]
                loc_q = loc_q+ls[2]+':'+ls[8]+';'

        elif line[:6] == 'Sbjct ':
            ls = line.split(' ')
            if len(ls) == 9:
                sbjct = sbjct+ls[6]
                loc_s = loc_s+ls[2]+':'+ls[8]+';'
            if len(ls) == 8:
                sbjct = sbjct+ls[5]
                loc_s = loc_s+ls[2]+':'+ls[7]+';'
            if len(ls) == 7:
                sbjct = sbjct+ls[4]
                loc_s = loc_s+ls[2]+':'+ls[6]+';'

    OUTFILE.write(newrn+'\t'+query[1:-1]+'\t'+sbjct[1:-1]+'\t'+loc_q[1:-1]+'\t'+loc_s[1:-1]+'\n')
    INFOFILE.write(fn[:-4]+':\nhit: '+str(hit)+'\n')
    INFOFILE.write(fn[:-4]+':\nnohit: '+str(nohit)+'\n')
    OUTFILE.close()
    FILE.close()
    INFOFILE.close()
    #print(fn[:12]+':\nhit: '+str(hit))
    #print('nohit: '+str(nohit))
    return rlt_dir

# generate fasta format file from fastq file and run blast
def blastn(read):
    read_full_name = read.split('/')[-1]
    read_name = read_full_name.split('_R')[0]
    rname_len = len(read_name)
    r_name = read_full_name[:rname_len+3]
    fqFILE = isFastq(read)
    faFILE = open(temp_dir+r_name+'.fa','w')

    lnum = 0
    #count = 0
    #count2 = 0
    for line in fqFILE:
        lnum += 1
        line = line.strip('\n')
        #print(line)
        if lnum % 4 == 1:
            ls = line.split(' ')
            faFILE.write('>'+ls[0][1:]+'\n')
        if lnum % 4 == 2:
            faFILE.write(line + '\n')
            '''
            if len(line) < 60:
                count += 1
            if len(line) < 100:
                count2 += 2
                '''
    fqFILE.close()
    faFILE.close()
    print('blastn:')
    os.system('nohup blastn -task dc-megablast -db '+ref_file+' -query '+temp_dir+r_name+'.fa -outfmt 0 -out '+temp_dir+r_name+'.out ')

    return r_name

# get the position of the third base close to PAM
def getpos(ref_file, target_file):
    lib = {}
    ref = ''
    target = ''
    RFILE = open(ref_file, 'r')
    TFILE = open(target_file, 'r')

    for line in TFILE:
        line = line.strip('\n')
        if line[0] == '>':
            continue
        target = target + line.upper()
    
    word_len = len(target)

    for line in RFILE:
        line = line.strip('\n')
        if line[0] == '>':
            continue
        ref = ref+line.upper()
    ref_len = len(ref)
    i = ref_len
    j = 0

    while i >= word_len:
        word = ref[j:j+word_len]
        if word in lib.keys():
            lib[word] += ';' + str(j+1)
        else:
            lib[word] = str(j+1)
        j += 1
        i -= 1
# judge whether target in reference. If not, throw error and break!
    if target not in lib.keys():
        print('The target region  not in reference!!!')
        sys.exit()
# target position in reference
    pos = int(lib[target].split(';')[0])
    print('pos: '+str(pos))
# if there two target regions, throw error and break!
    if len(lib[target].split(';')) >= 2:
        print('There 2 target regions in reference!!!')
        sys.exit()
    if target[0:2] == 'GG' or target[0:2] == 'CC':
        return pos + 5
    if target[-2:] == 'GG' or target[-2:] == 'CC':
        return pos + word_len - 5
    else:
        print('The target is without PAM. Please add PAM to the target sequence: '+target_file)
        sys.exit()


# draw figure
def getFigure(x, y, title, x_label, y_label, name, fn):
    fig_dir = out_dir + fn[:-3]+'_figures/'
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig, ax = plt.subplots(figsize = (10, 10))
    ax.yaxis.grid(True, linestyle = '-', which = 'major', color = 'lightgrey', alpha = 0.5)
    ax.xaxis.grid(True, linestyle = '-', which = 'major', color = 'lightgrey', alpha = 0.5)
    ax.set_axisbelow('True')
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
#Bar Plot
    plt.bar(x, y, facecolor = 'red', edgecolor = 'black')
    plt.savefig(fig_dir+name, format = 'jpg')


# make temp files directory
temp_dir = ''
# make figure files directory
fig_dir = ''
if os.path.exists(out_dir):
    temp_dir = out_dir+'cas_temp/'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
else:
    os.makedirs(out_dir)
    temp_dir = out_dir+'cas_temp/'
    os.makedirs(temp_dir)
# make blast database
os.system('makeblastdb -in '+ref_file+ ' -dbtype nucl -parse_seqids')
# get the position of the third base close to PAM
rang = getpos(ref_file, target_file)
print('rang: '+str(rang))
print('temp: '+temp_dir)
# invoke blastn function for read1
r_name = blastn(read1)
print('r_name: '+r_name)
result_dir = clear(temp_dir+r_name+'.out', r_name)
print('r_dir: '+result_dir)
#getindel(temp_dir+r_name[:-3]+'.clear.out', rang, result_dir, r_name)
# invoke blastn function for read2
r_name = blastn(read2)
result_dir = clear(temp_dir+r_name+'.out', r_name)
getindel(temp_dir+r_name[:-3]+'.clear.out', rang, result_dir, r_name)
#r_ls = [257, 442, 508, 782, 543, 603, 685]
