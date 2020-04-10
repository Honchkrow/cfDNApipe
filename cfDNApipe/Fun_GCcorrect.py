# -*- coding: utf-8 -*-
"""
Created on Wed Apr 8 14:25:33 2020

@author: Jiaqi Huang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError, wig2df, correctReadCount
import pandas as pd
import numpy as np
import os
from .Configure import Configure

__metaclass__ = type


class GCCorrect(StepBase):
    def __init__(
            self,
            readInput=None, # list
            gcwigInput=None, # list
            readtype=None, # int
            corrkey=None, # char
            outputdir=None,  # str
            stepNum=None,
            readupstream=None,
            gcupstream=None,
            **kwargs
    ):
        """
        This function is used for processing GC correction on read count data in wig or csv/txt files.

        GCcorrect(readInput=None, gcwigInput=None, outputdir=None, stepNum=None, readupstream=None, gcupstream=None)
        {P}arameters:
            readInput: list, paths of input files of read counts.
            gcwigInput: list, paths of wig files of gc contents.
            readtype: int, file type of readInput, 1 for .wig, 2 for .txt/.csv.; 1 is set by default.
            corrkey: char, type of GC correction, "-" for minus, "/" for divide, None for process without GC correction.
            outputdir: str, output result folder, None means the same folder as input files.
            stepNum: Step number for folder name.
            readupstream: Not used parameter, do not set this parameter.
            gcupstream: Not used parameter, do not set this parameter.
        """
            
        super(GCCorrect, self).__init__(stepNum, readupstream)

        if readupstream is None or gcupstream is None:
            self.setInput("readInput", readInput)
            self.setInput("gcwigInput", gcwigInput)
            self.checkInputFilePath()
            
            if outputdir is None:
                self.setOutput("outputdir", os.path.dirname(
                    os.path.abspath(self.getInput("readInput")[0])))
            else:
                self.setOutput("outputdir", outputdir)
        
        else:
            Configure.configureCheck()
            readupstream.checkFilePath()
            gcupstream.checkFilePath()

            if readupstream.__class__.__name__ == "runCounter":
                self.setInput("readInput", readupstream.getOutput("wigOutput"))
            elif readupstream.__class__.__name__ == "fraglenCounter":
                self.setInput("readInput", readupstream.getOutput("txtOutput"))
            else:
                raise commonError("Parameter upstream must from runCounter or fraglenCounter.")
            if gcupstream.__class__.__name__ == "runCounter":
                self.setInput("gcwigInput", gcupstream.getOutput("wigOutput"))
            else:
                raise commonError("Parameter upstream must from runCounter.")
            self.checkInputFilePath()

            self.setOutput("outputdir", self.getStepFolderPath())

        self.setOutput("txtOutput", [os.path.join(self.getOutput(
            "outputdir"), self.getMaxFileNamePrefixV2(x)) + ".txt" for x in self.getInput("readInput")])
            
        self.setOutput("plotOutput", [os.path.join(self.getOutput(
            "outputdir"), self.getMaxFileNamePrefixV2(x)) + ".png" for x in self.getInput("readInput")])

        finishFlag = self.stepInit(readupstream)
        
        multi_run_len = len(self.getInput("readInput"))
        
        if not finishFlag:
            gc_df = wig2df(self.getInput("gcwigInput")[0])
            if 
            for i in range(multi_run_len):
                print("Now, processing", self.getMaxFileNamePrefixV2(
                    self.getInput("readInput")[i]), "...")
                read_df = wig2df(self.getInput("readInput")[i])
                correctReadCount(
                    read_df, gc_df, self.getOutput("txtOutput")[i], self.getOutput("plotOutput")[i])

        self.stepInfoRec(cmds=[], finishFlag=finishFlag)
