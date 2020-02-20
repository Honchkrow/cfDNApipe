# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 10:38:07 2019

@author: zhang
"""

from cfDNApipe import *

# set global configure
pipeConfigure(threads=100,
              genome="hg19",
              refdir="/home/zhangwei/Genome/hg19_bismark",
              outdir="/home/zhangwei/pipeline-WGBS-comp",
              data="WGBS",
              type="paired",
              case="case",
              ctrl="ctrl",
              build=True)

# case processing
switchConfigure("case")
res_case1 = inputprocess(inputFolder = r"/home/zhangwei/pipeline-WGBS-cc/raw/case")
res_case2 = fastqc(upstream = res_case1)
res_case3 = identifyAdapter(upstream = res_case1)
res_case4 = adapterremoval(upstream = res_case3)
res_case5 = bismark(upstream = res_case4)
res_case6 = bismark_deduplicate(upstream = res_case5)
res_case7 = bismark_methylation_extractor(upstream = res_case6)
res_case8 = compress_methyl(upstream = res_case7)
res_case9 = calculate_methyl(upstream = res_case8)
res_case10 = bamsort(upstream = res_case6)
res_case11 = bam2bed(upstream = res_case10)
res_case12 = fraglenplot(upstream = res_case11)

# ctrl processing
switchConfigure("ctrl")
res_ctrl1 = inputprocess(inputFolder = r"/home/zhangwei/pipeline-WGBS-cc/raw/control")
res_ctrl2 = fastqc(upstream = res_ctrl1)
res_ctrl3 = identifyAdapter(upstream = res_ctrl1)
res_ctrl4 = adapterremoval(upstream = res_ctrl3)
res_ctrl5 = bismark(upstream = res_ctrl4)
res_ctrl6 = bismark_deduplicate(upstream = res_ctrl5)
res_ctrl7 = bismark_methylation_extractor(upstream = res_ctrl6)
res_ctrl8 = compress_methyl(upstream = res_ctrl7)
res_ctrl9 = calculate_methyl(upstream = res_ctrl8)
res_ctrl10 = bamsort(upstream = res_ctrl6)
res_ctrl11 = bam2bed(upstream = res_ctrl10)
res_ctrl12 = fraglenplot(upstream = res_ctrl11)







###########
from cfDNApipe import *
Configure.setData('WGBS')
Configure.setThreads(100)
Configure.setGenome("hg19")
Configure.setRefDir(r"/home/zhangwei/Genome/hg19_bismark")
Configure.setOutDir(r'/home/zhangwei/pipeline-WGBS')
Configure.pipeFolderInit()
Configure.refCheck(build = True)

res_case1 = inputprocess(inputFolder = r"/home/zhangwei/pipeline-WGBS/raw")
res_case2 = fastqc(upstream = res_case1)
res_case3 = identifyAdapter(upstream = res_case1)
res_case4 = adapterremoval(upstream = res_case3)
res_case5 = bismark(upstream = res_case4)
res_case6 = bismark_deduplicate(upstream = res_case5)
res_case7 = bismark_methylation_extractor(upstream = res_case6)
res_case8 = compress_methyl(upstream = res_case7)
res_case9 = calculate_methyl(upstream = res_case8)
res_case10 = bamsort(upstream = res_case6)
res_case11 = bam2bed(upstream = res_case10)
res_case12 = fraglenplot(upstream = res_case11)

