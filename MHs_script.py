# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 17:35:22 2022

@authors: Rodrigues, M.P. & Rodrigues, M.L.B.
"""
from datetime import datetime, date
import os
from sys import exit
from collections import Counter  # in the FUTURE, observe if there is a change in this method

t_ini = datetime.now()

SEP = '\n' + '='*115 + '\n'

## DEFINITIONS
dir_base = os.getcwd()
dir_SNPs = dir_base + '\\SNPs_in_MH_vcfs'
dir_MH = dir_base + '\\MHs_vcfs'

# DEFINE THE DIMENSION OF THE MICROHAPLOTYPE
dim_mH = 300

# Define the output files prefixes
prefix_SNPs_in_MH = 'SNPs_in_MH_' #files with selected positions
prefix_MH = 'MH_' #files with MHS

print(SEP)
print('Working directory:')
print(dir_base)



"""
###################################################################################################
##                                       DEFINE FUNCTIONS                                        ##
###################################################################################################
"""

def starting(scope):
    print(SEP)
    print('PROCESSING - ' + scope + ' ...............................................\n')



def record_line(file_2, mode, line_rec):
    with open(file_2, mode) as arc_VCF:
        arc_VCF.write(line_rec)



def sel_mH(file_analyzing):
    # Function that selects SNPs in MH and assembles VCF with every respective positions

    file_VCF_SNP = open(dir_base + '\\' + file_1, 'r') #specifies the input file being processed
    vcf_sel_mH = open(dir_SNPs + '\\' + prefix_SNPs_in_MH+file_1, 'w')  #creates the output file

    nr_POS = 0
    POS_base_ant = 0
    POS_current = 0
    POS_base = 0
    qtd_POS = 0
    mH = 0

    global line_SNP

    for line_SNP in file_VCF_SNP: # SELECT AND WRITE THE LINES OF THE MICROHAPLOTYPES
        if (line_SNP[0:1]) != '#':  # verification to avoid an "out of range" error, first lines
            line_current = line_SNP.split('\t')
            POS_current = int(line_current[1])
            CHR = line_current[0]

            if POS_current-POS_base > dim_mH:
                POS_base = POS_current
                line_base = line_current
            else:
                if POS_base != POS_base_ant:
                    mH = mH + 1
                    qtd_POS = qtd_POS + 1
                    #note: 'line_base' was defined by previous cycle
                    line_base[2] = 'MH' + f"{mH:06d}" + '_CHR' + CHR
                    line_base[5] = '.'
                    line_base[6] = '.'
                    line_base[7] = '.'
                    line_base[8] = 'GT'
                    line_current[2] = 'MH' + f"{mH:06d}" + '_CHR' + CHR
                    line_current[5] = '.'
                    line_current[6] = '.'
                    line_current[7] = '.'
                    line_current[8] = 'GT'
                    vcf_sel_mH.write('\t'.join(line_base[0:]))  # writes first line of SNPs in MH
                    vcf_sel_mH.write('\t'.join(line_current[0:]))  # writes the second line
                    POS_base_ant = POS_base
                else:
                    line_current[2] = 'MH'+f"{mH:06d}" + '_CHR' + CHR
                    line_current[5] = '.'
                    line_current[6] = '.'
                    line_current[7] = '.'
                    line_current[8] = 'GT'
                    qtd_POS = qtd_POS + 1
                    vcf_sel_mH.write('\t'.join(line_current[0:])) #writes third and following lines

            nr_POS = nr_POS + 1

            """ LIMITS PROCESSING TO A FEW LINES """
            #if nr_POS > 50:
            #    break

        else:  # SETS AND WRITES THE HEADING LINES - SNPs VCF
            if (line_SNP[0:2]) == '##':
                continue
            else:
                vcf_sel_mH.write('##fileformat=VCFv4.2\n')
                vcf_sel_mH.write('##fileDate=' + date.today().strftime('%Y%m%d')+ '\n')
                vcf_sel_mH.write('##source=SNPs_in_MHs_script\n')
                vcf_sel_mH.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                vcf_sel_mH.write(line_SNP)

    print('processed lines:\t\t\t\t\t' + str(nr_POS))
    print('Number of MHs formed:\t\t\t\t' + str(mH))
    print('Quantity of SNPs composing MHs:\t\t' + str(qtd_POS + mH))

    # ASSEMBLES TABLE OF RESULTS
    arc_results1 = open(dir_SNPs + '\\MHs_REPORT.txt', 'a')
    arc_results1.write(str(CHR) + '\t' + str(nr_POS) + '\t\t\t' + str(mH) + '\t' + str(qtd_POS+mH) + '\n')

    # CLOSE OPENED FILES
    file_VCF_SNP.close()
    vcf_sel_mH.close()
    arc_results1.close()



def group_mH(MH_file_analyzed):
    ##  FUNCTION that groups by MH, using the MH.VCF file (with the MH positions)
     #  And assemble VCF with MH grouped by line

    count = 0
    nr_mH = 0
    mH_base = ''
    Ref_mH = ''
    Alt_mH = ''
    allele_dict_MH = {}
    mHz = 0

    global out_file
    out_file = dir_MH + '\\' + file_2.replace(prefix_SNPs_in_MH, prefix_MH)
    vcf_MH = open(out_file, 'w')  #creates the output file
    file_VCF_MH = open(dir_SNPs + '\\' + file_2, 'r') #opens the input file

    global line_MH
    for line_MH in file_VCF_MH:

        if (line_MH[0:1]) != '#':
            current_line_MH = line_MH.split()
            global CHR
            CHR = current_line_MH[0]
            global POS
            POS = current_line_MH[1]
            global MH_line
            MH_line = current_line_MH[2]
            REF = current_line_MH[3]
            ALT = current_line_MH[4]
            allele_dict = {}

            ## ASSEMBLY THE ALLELE DICTIONARY, EXCHANGE 0 or 1 with the respective letter REF or ALT
            ## Here is the treatment of non-existent ./.
            x = 0
            for allele in current_line_MH[9:]: #creates allele_dict by line, by individuals
                if allele[0] == '0':
                    alle_ref = current_line_MH[3]
                elif allele[0] == '1':
                    alle_ref = current_line_MH[4]
                else:
                    alle_ref = '.'

                if allele[2] == '0':
                    alle_alt = current_line_MH[3]
                elif allele[2] == '1':
                    alle_alt = current_line_MH[4]
                else:
                    alle_alt = '.'

                x = x + 1
                allele_dict[x] = [alle_ref,alle_alt]  # square brackets create dictionary with lists, not tuples
            ##  END OF ASSEMBLY OF THE ALLELE DICTIONARY

            ##  HAPLOTYPE PROCESSING
            if MH_line != mH_base:
                nr_mH = nr_mH + 1

                ##  FINALIZATION OF THE MICROHAPLOTYPES, ALLELE COUNTING AND GENERATION OF REF AND ALT
                if count > 0:
                    mHz = mHz + (assemble_MH(allele_dict_MH)) #call FUNCTION like this to get RETURN

                allele_dict_MH = allele_dict.copy()

            else:
                Ref_mH = Ref_mH + REF
                Alt_mH = Alt_mH + ALT

                global last_line
                last_line = line_MH

                for key in allele_dict_MH.keys():
                    allele_dict_MH[key][0] = allele_dict_MH[key][0] + allele_dict[key][0]
                    allele_dict_MH[key][1] = allele_dict_MH[key][1] + allele_dict[key][1]

                Ref_mH = REF
                Alt_mH = ALT

            mH_base = MH_line

            ##  END OF HAPLOTYPE PROCESSING

            count = count + 1

        else:              # SETS AND WRITES THE HEADING LINES - MHs VCF
            if line_MH[0:7] == '##sourc':
                record_line(out_file, 'a', '##source=MHs_script\n')
            else:
                record_line(out_file, 'a', line_MH)

    vcf_MH.close()

    mHz = mHz + (assemble_MH(allele_dict_MH))

    file_report = dir_MH + '\\MHs_REPORT.txt'
    lin_report = str(CHR) + '\t' + str(count) + '\t\t' + str(nr_mH) + '\t' + str(mHz) + '\t' + str(nr_mH - mHz)
    record_line(file_report, 'a', lin_report + '\n')
    print('processed lines:\t\t',count)
    print('mHs descartados:\t\t', mHz)
    print('number of MHs formed:\t', nr_mH)

    file_VCF_MH.close()



def assemble_MH(allele_dict_MH):
    ##  FUNCTION that analyzes the microhaplotype formed AND assembles the VCF line

    out_line = ''

    allele_MH = [xx for key in allele_dict_MH for xx in allele_dict_MH[key]] # DICT COMPREHENSION
    allele_MH = [yy for yy in allele_MH if '.' not in yy]  # DELETION OF NON-EXISTENT DATA
    alleles_qty = Counter(allele_MH).most_common()  #COUNTS THE ALLELES, CLASSIFYING FROM MOST TO LESS FREQUENT
    alleles_list = [x for x,y in alleles_qty]

    if len(alleles_qty) == 1:
        ctg_z = 1 #monomorphic MHs count
    else:
        ctg_z = 0

    ##  GENERATING THE ALLELE INDEX
    y = 0
    Ref_Alt = {}
    for nr_alelo in alleles_qty:
        Ref_Alt[y] = alleles_qty[y][0]
        y = y + 1

    Ref_Alt_Inv = {v: k for k, v in Ref_Alt.items()}
    #https://acervolima.com/python-maneiras-de-inverter-o-mapeamento-do-dicionario/
    Reference = alleles_list[0]
    Alternate = ','.join([allele for allele in alleles_list if allele!= Reference])

    out_line_1 = last_line.split()[0:9]
    out_line_2 = ''

    ##  CONSTRUCTION OF THE OUTPUT LINE
    for key2 in allele_dict_MH:
        allele_MH_1 = allele_dict_MH[key2][0]
        allele_MH_2 = allele_dict_MH[key2][1]

        if '.' in allele_MH_1 or '.' in allele_MH_2:  #assembling the non-existent allele
            allele_local = '.|.'
        else:
            allele_local = str(Ref_Alt_Inv[allele_MH_1])+'|'+str(Ref_Alt_Inv[allele_MH_2])

        out_line_2 = out_line_2 + '\t' + allele_local

    ##  ASSEMBLY OF THE OUTPUT LINE
    out_line_1[3] = Reference
    out_line_1[4] = Alternate
    out_line = out_line_1 + out_line_2.split()
    line_VCF = '\t'.join(out_line)
    if out_line_1[4] != '':   # RECORDS only those with alternates <> ""
        record_line(out_file, 'a', line_VCF + '\n')
    ##  END OF CONSTRUCTION OF THE OUTPUT LINE

    return ctg_z  # if there is a RETURN command, it must be on the last line of the function



#%%
"""
###################################################################################################
##                           CREATING DIRECTORIES and RESULT FILES                               ##
##          "SNPs_in_MH_vcfs" and "MHs_vcfs" will be created if they do not exist                ##
###################################################################################################
"""
starting('Verifying and preparing the directories')

#  CREATES DIRECTORY FOR SNPS IN MH VCF FILES, IF IT DOES NOT ALREADY EXIST
if os.path.exists(dir_SNPs):
    print('Directory', dir_SNPs, 'already exists: ', os.path.isdir(dir_SNPs))
else:
    print('Directory "SNPs_in_MH_vcfs" does not exist and will be created\n...')
    os.mkdir(dir_SNPs)
    print('Directory', dir_SNPs, 'was created:', os.path.isdir(dir_SNPs),'\n')

#  CREATE DIRECTORY FOR MHS VCF FILES, IF IT DOES NOT ALREADY EXIST
if os.path.exists(dir_MH):
    print('Directory', dir_MH, 'already exists: ', os.path.isdir(dir_MH))
else:
    print('Directory "MHs_vcfs" does not exist and will be created\n...')
    os.mkdir(dir_MH)
    print('Directory', dir_MH, 'was created:', os.path.isdir(dir_MH), '\n')
print(SEP)

print('files present in the directory where this script is:\n')
for arc in os.listdir(dir_base):
    print(arc)

# ASSEMBLY OF THE TABLE OF RESULTS HEADING - SNPS IN MHS
arc_results1=open(dir_SNPs+'\\MHs_REPORT.txt','w') # use 'a' for append
arc_results1.write('Dimension of the Microhaplotypes: '+str(dim_mH)+'\n\n')
arc_results1.write('CHR\tTotal_Positions\t\tMH\tPositions_in_MH\n')
arc_results1.close()

## ASSEMBLY OF THE TABLE OF RESULTS HEADING - MHS
arc_results2 = open(dir_MH+'\\MHs_REPORT.txt','w')
arc_results2.write('MICROHAPLOTYPE FORMATION REPORT \n')
arc_results2.write('\nCHR\tSNPs_in_MH\tMH\tMH_Mono\tMH_Poli\n')
arc_results2.close()



#%%
"""
###################################################################################################
##                                   SELECTING SNPS IN MHS                                       ##
###################################################################################################
"""
starting('Selecting SNPS in MHS')

# CREATES A LIST WITH ONLY THE EXISTING VCF FILES IN THE CURRENT DIRECTORY
vcfs1 = []
for x in os.listdir(dir_base):
    if x.endswith('.vcf'):
        vcfs1 = vcfs1 + [x]

print('Origin  ',dir_base)
print('Destiny ',dir_SNPs)

for file_1 in vcfs1:
    print('\n\nfile in process:',file_1)
    sel_mH(file_1)

arc_results1.close()



#%%
""" 
###################################################################################################
##                                         CREATING MHS                                          ##
###################################################################################################
"""
starting('Assembling the MHS')

# CREATES A LIST WITH ONLY THE EXISTING VCF FILES IN THE DIRECTORY "SNPs_in_MH_vcfs"
vcfs2 = []
for x in os.listdir(dir_SNPs):
    if x.endswith('.vcf'):
        vcfs2 = vcfs2 + [x]

print('\nOrigin  ',dir_SNPs)
print('Destiny ',dir_MH)

for file_2 in vcfs2:
    print('\nfile in process: ',file_2)
    group_mH(file_2)

arc_results2.close()



#%%
starting('Ending')
t_end = datetime.now()

print('all processing time: ', t_end - t_ini)
exit
