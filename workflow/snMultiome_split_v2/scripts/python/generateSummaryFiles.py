#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Jan 30, 2024
usage: python generateSummaryFiles.py metric_summary
Modified from the generateSummaryFiles.py from SF by Ying

"""


import glob, xlsxwriter, csv, sys, ntpath, os
from shutil import copyfile
from sys import argv

rawdataPath = argv[1]
project = argv[2]
metricsPath = 'cellranger_finalreport/'
summaryPath = 'cellranger_finalreport/summaries/'
#print (metricsPath, summaryPath)

def main(arg1='metric_summary'):
    createMetricsSummary(arg1)
    copyWebSummary()

def is_number(s):
    try:
        float(s)
        if '_' in s:
            return False
        return True
    except ValueError:
        return False

def createMetricsSummary(arg1):
    try:
        os.makedirs(metricsPath)
    except OSError:
        if not os.path.isdir(metricsPath):
            raise
    files = glob.glob(f'{rawdataPath}/*/outs/summary.csv')
    workbook = xlsxwriter.Workbook(metricsPath + project+'.xlsx')
    worksheet = workbook.add_worksheet("metrics_summary")
    worksheet.set_column(0, 12, 10.1)
    worksheet.set_column(13, 16, 12.2)
    worksheet.set_column(17,20, 10)

    formatNum = workbook.add_format({'num_format': '#,###'})
    formatDec = workbook.add_format({'num_format': '#.####'})
    formatPer = workbook.add_format({'num_format': '0.00%'})
    formatHead = workbook.add_format({'bold': True, 'italic': True, 'text_wrap': True, 'align': 'center'})
    row = 1
    samples = list()
    for filename in files:
        with open(filename, 'r') as csvfile:
            f = csv.reader(csvfile, delimiter=',', quotechar='"')
            header = next(f)
            line = next(f)
            #worksheet.write(row, 0, filename.split('/')[1])
            samples.append(filename.split('/')[-3])
            col = 0
            for i in line:
                i = i.strip('"')
                if '%' in i:
                    worksheet.write(row, col, float(i.strip('%'))/100, formatPer)
                elif is_number(i):
                    if i.isdigit():
                        worksheet.write(row, col, int(i.replace(',','')), formatNum)
                    else:
                        worksheet.write(row, col, float(i.replace(',','')), formatDec)
                else:
                    worksheet.write(row, col, i)
                col += 1
            row += 1

    col = 0
    row = 0
    print(samples)
    for i in header:
        worksheet.write(row, col, i, formatHead)
        col += 1

    #for i in samples:
    #    worksheet = workbook.add_worksheet(i)

    workbook.close()

def copyWebSummary():
    try:
        os.makedirs(summaryPath)
    except OSError:
        if not os.path.isdir(summaryPath):
            raise
    files = glob.glob(f'{rawdataPath}/*/outs/web_summary.html')
 #   print("now let's copy files")
 #   print(files)
    for filename in files:
        print(filename.split('/')[-3])
        copyfile(filename, '%s/%s_web_summary.html' % (summaryPath, filename.split('/')[-3]))
        sample_name = filename.split('/')[-3]
  #      print(sample_name)
        if sample_name != "AggregatedDatasets":
            gex_summary = glob.glob(f"{rawdataPath}/{sample_name}/SC_ATAC_GEX_COUNTER_CS/SC_ATAC_GEX_COUNTER/GEX_SUMMARIZE_REPORTS/*/*/files/web_summary.html")
            copyfile(gex_summary[0], '%s/%s_GEX_web_summary.html' % (summaryPath, sample_name))
            atac_summary = glob.glob(f"{rawdataPath}/{sample_name}/SC_ATAC_GEX_COUNTER_CS/SC_ATAC_GEX_COUNTER/_SC_ATAC_REPORTER/CREATE_WEBSUMMARY/*/files/web_summary.html") 
            copyfile(atac_summary[0], '%s/%s_ATAC_web_summary.html' % (summaryPath, sample_name))
    
    
if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        main(sys.argv[1])
