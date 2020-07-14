# -*- coding: utf-8 -*-
# Created on Mar 28, 2020
# Author: Maxwell4444
# All right reserved

import numpy as np
import sys
import os

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QLabel, QPushButton, QMessageBox, QComboBox
from PyQSM.QSMCalAll import QSMCalAll, QSMBet, getBasepath
import nibabel as nib

def ComboboxFastCreate(strlist,curind=0):
    combobox = QComboBox()
    for key in strlist:
        combobox.addItem(key)
    combobox.setCurrentIndex(curind)
    return combobox   

def Save1Nii(img,affine, filepath):
    nib.save(nib.Nifti1Image(img, affine),filepath)
         
def Load1Nii(filepath):
    fimg  = nib.load(filepath)
    data = np.asarray(fimg.dataobj).astype(float)
    return data    
class QSMClass(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self,parent)
        self.createAppUI()
        self._data, self._fimg, self._dinfo = [], [], []
        
    def getAppName(self):
        return 'QSM'

    def createAppUI(self):
        #######################################################################
        # T1/T2/T2* mapping block
        appUI = self
        self._appUI = appUI
        button_layout = QtWidgets.QGridLayout()
        button_layout.setContentsMargins(0,0,0,0)
        appUI.setLayout(button_layout)

        index = 0
        button_layout.addWidget(QtWidgets.QLabel(' Import:'),index,1)
        
        button_importnii = QtWidgets.QPushButton(' ImportNii')
        button_layout.addWidget(button_importnii,index+1,1)
        button_importnii.clicked.connect(self.ImportData)
        
        ####################################################################
        ## create QSM setting UI module
        index +=2
        appUI, QSMConfigUI = self.CreateQSMConfigUI(appUI)
        appUI._QSMConfigUI = QSMConfigUI
        button_layout.addWidget(QSMConfigUI,index,1,1,3)  

        return appUI        
    def CreateQSMConfigUI(self, appUI):
        UI = QtWidgets.QWidget()
        button_layout = QtWidgets.QGridLayout()
        button_layout.setContentsMargins(0,0,0,0)
        UI.setLayout(button_layout)  # change       
        
        ########################################################################
        ## general param
        index = 0
        button_layout.addWidget(QLabel('QSMParam:'),index,1,1,1)
    
        combobox_dataorder = ComboboxFastCreate(['mm..pp..', 'mpmp..'],curind=0)
        button_calall = QPushButton('QSMCalAll')
        button_layout.addWidget(QLabel('   DataOrder:'), index+1,1)
        button_layout.addWidget(combobox_dataorder,index+1,2)
        button_layout.addWidget(button_calall, index+1,3)
        
        button_batchqsmcal = QPushButton('BatchQSMCal')
        button_Loadqsmresults = QPushButton('LoadQSMResult')
        button_layout.addWidget(button_batchqsmcal, index+2,1)
        button_layout.addWidget(button_Loadqsmresults, index+2,2)
    
        button_calall.clicked.connect(self.QSMCalAll)
        button_batchqsmcal.clicked.connect(self.BatchQSMCal)
        button_Loadqsmresults.clicked.connect(self.LoadQSMResults)
        UI._combobox_dataorder = combobox_dataorder
        
        ########################################################################
        ## detail step param
    #         index += 2
    #         combobox_unwarpAlgo = ComboboxFastCreate(['laplacian'],curind=0)
    #         combobox_BFRemovalAlgo = ComboboxFastCreate(['PDF'], curind=0)
    #         combobox_DipoleInversionAlgo = ComboboxFastCreate(['MEDI'], curind=0)
    #         
    #         button_layout.addWidget(combobox_unwarpAlgo,index+1,1)
    #         button_layout.addWidget(combobox_BFRemovalAlgo,index+2,1)
    #         button_layout.addWidget(combobox_DipoleInversionAlgo,index+3,1)
    #         
    #         UI._combobox_unwarpAlgo = combobox_unwarpAlgo
    #         UI._combobox_BFRemovalAlgo = combobox_BFRemovalAlgo 
    #         UI._combobox_DipoleInversionAlgo = combobox_DipoleInversionAlgo
    #         
    #         button_calfitting = QPushButton('CalFitting')
    #         button_calunwrapping = QPushButton('CalUnwrapping')
    #         button_calBFRemoval = QPushButton('CalBFRemoval')
    #         button_calBFRemoval.setToolTip('Background Field Removal')
    #         button_DipoleInv = QPushButton('DipoleInversion')
    #         
    #         button_layout.addWidget(button_calfitting, index, 2)
    #         button_layout.addWidget(button_calunwrapping, index+1, 2)
    #         button_layout.addWidget(button_calBFRemoval, index+2, 2)
    #         button_layout.addWidget(button_DipoleInv, index+3, 2)
    # 
    #         button_calfitting.clicked.connect(self.QSMCalFitting)
    #         button_calunwrapping.clicked.connect(self.QSMUnwrapping)
    #         button_calBFRemoval.clicked.connect(self.QSMBFRemoval)
    #         button_DipoleInv.clicked.connect(self.QSMDipoleInv)
        return appUI, UI 

    def BatchQSMCal(self):
        pass
    def LoadQSMResults(self):
        pass    
    def ImportData(self):
        dialog = QtWidgets.QFileDialog()
        startpath = 'D:/eclipse-workspace/OncoImageAnalysis/src/ExampleData/QSM/batch_test/case2/010_QSM_3D_high_res_2mm_merge.nii.gz'
#         startpath = r'.' 
        filelist, ext = dialog.getOpenFileNames(parent=None, caption='Import Multi Nii',directory=startpath,filter='Nii(*.nii.gz);;Nii(*.nii);;All(*.*)')
        if filelist == []: return
        self._datapath = filelist[0]

        def getDInfo(dinfopath):
            dinfo = []
            if os.path.isfile(dinfopath):
                # how to load dinfo  
                arch = np.load(dinfopath)
                if 'dinfo' in arch.keys():
                    dinfo = arch['dinfo'].tolist() 
            return dinfo
        def LoadNiiData(fdwi):
            fimg  = nib.load(fdwi)
            data = np.asarray(fimg.dataobj).astype(float)
            dinfo = getDInfo(getBasepath(fdwi)+'.npz')
            return data, fimg, dinfo
        
        self._rawdata, self._fimg, self._dinfo = LoadNiiData(self._datapath)
        self._pixel_spacing = self._fimg.header['pixdim'][1:4]        
        self._tr = self._dinfo['tr']
        self._te = self._dinfo['te']
        self._fa = self._dinfo['fa']
    def QSMCalAll(self):
        try:
            self.doQSMCalAll()
            self.SaveNii()
        except:
            print('QSMCallAll error')
    def SaveNii(self):
        for key in self._outputParam:
            Save1Nii(self._outputParam[key], self._fimg.affine, getBasepath(self._datapath)+'_'+str(key)+'.nii.gz')
    def doQSMCalAll(self):
        dataorder = self._appUI._QSMConfigUI._combobox_dataorder.currentText()
        
        # TE = np.array([0.0036,0.0095,0.0154,0.0213,0.0273,0.0332,0.0391,0.0450,   0.0036,0.0095,0.0154,0.0213,0.0273,0.0332,0.0391,0.0450])
        TE = np.array(self._te)
        volnum = int(self._rawdata.shape[3])
        if volnum%2 != 0: 
            print('data should have same number of mag and pha data')
            return
        def ReScalePhase(pha): return np.round(pha) / 4095 * np.pi
        if dataorder == 'mm..pp..':
            pha = self._rawdata[...,int(volnum/2):]
            #pha = pha / (np.abs(pha).max()) * np.pi;
            pha = ReScalePhase(pha)
            mag = self._rawdata[...,:int(volnum/2)]
            TE = TE[:int(volnum/2)]
        elif dataorder == 'mpmp..':
            pha = self._rawdata[...,1::2]        
            mag = self._rawdata[..., ::2]
            #pha = pha / (np.abs(pha).max()) * np.pi;
            pha = ReScalePhase(pha)
            TE = TE[..., ::2]
        
#         refdatareader = copy.copy(self._image2Dframe._imgaxes[0]._datareader)
#         qsm_filepath = refdatareader._filepath
#         voxel_size = refdatareader._pixel_spacing #np.array([0.9375,0.9375,2.0])
        qsm_filepath, voxel_size = self._datapath, self._pixel_spacing
        if 'B0_dir' in self._dinfo.keys():
            B0_dir = self._dinfo['B0_dir'] 
            CF = self._dinfo['CF']
            PrecessionIsClockwise = self._dinfo['PrecessionIsClockwise']
        else: 
            B0_dir = np.array([0,0,1])
            CF = 123225803
            PrecessionIsClockwise = -1
        
#         if PrecessionIsClockwise == -1:
#             pha = -pha
             
        # param['mag']
        qsm_maskpath = getBasepath(qsm_filepath) +'_mask.nii.gz'
        qsm_magfilepath = getBasepath(qsm_filepath) +'_magsum.nii.gz'
        magsum = np.sqrt(np.sum(mag**2,axis=3))
#         datareader = DataReader()
#         datareader.SaveNiftyFilebyData(magsum, qsm_magfilepath, qsm_filepath)
        Save1Nii(magsum.astype(int), self._fimg.affine, qsm_magfilepath)
        QSMBet(qsm_magfilepath, qsm_maskpath,f=0.5)
        
#         datareader.LoadNiftyFile(qsm_maskpath)
#         Mask = datareader._data.astype(int)
        Mask = Load1Nii(qsm_maskpath)
        try:
            self._outputParam = QSMCalAll(mag,pha,voxel_size, B0_dir, TE, CF, Mask=Mask,qsm_filepath='')
        except:
            return    

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    qsmapp = QSMClass()
    qsmapp.show()
    sys.exit(app.exec_())