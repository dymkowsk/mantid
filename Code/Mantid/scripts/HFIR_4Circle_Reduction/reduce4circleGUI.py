#pylint: disable=invalid-name
################################################################################
#
# MainWindow application for reducing HFIR 4-circle 
#
################################################################################
import sys
import os
import numpy

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import reduce4circleControl as r4c


try:
    import mantid
except ImportError:
    sys.path.append('/home/wzz/Mantid/Code/debug/bin/')
    import mantid
finally:
    import mantid.simpleapi as api
    import mantid.kernel
    from mantid.simpleapi import AnalysisDataService
    from mantid.kernel import ConfigService

import fourcircle_utility as fcutil


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

from Ui_MainWindow import Ui_MainWindow #import line for the UI python class


class MainWindow(QtGui.QMainWindow):
    """ Class of Main Window (top)
    """
    def __init__(self, parent=None):
        """ Initialization and set up
        """
        # Base class
        QtGui.QMainWindow.__init__(self,parent)

        # UI Window (from Qt Designer)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Mantid configuration
        config = ConfigService.Instance()
        self._instrument = config["default.instrument"]

        # Event handling definitions
        # Tab 'Data Access'
        self.connect(self.ui.pushButton_browseLocalDataDir, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_spice_data)
        self.connect(self.ui.pushButton_testURLs, QtCore.SIGNAL('clicked()'),
                     self.do_test_url)
        self.connect(self.ui.pushButton_ListScans, QtCore.SIGNAL('clicked()'),
                     self.do_list_scans)
        self.connect(self.ui.pushButton_downloadExpData, QtCore.SIGNAL('clicked()'),
                     self.do_download_spice_data)
        self.connect(self.ui.comboBox_mode, QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.change_data_access_mode)

        # Tab 'Advanced'
        self.connect(self.ui.pushButton_browseLocalCache, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_cache_dir)

        # Tab ...
        self.connect(self.ui.pushButton_load, QtCore.SIGNAL('clicked()'), self.doLoad)


        self.connect(self.ui.pushButton_plotScan, QtCore.SIGNAL('clicked()'),
                self.doPlotScanPt)

        self.connect(self.ui.pushButton_prevScan, QtCore.SIGNAL('clicked()'),
                self.doPlotPrevScanPt)

        self.connect(self.ui.pushButton_nextScan, QtCore.SIGNAL('clicked()'),
                self.doPlotNextScanPt)

        # Event handling for tab 'calculate ub matrix'
        self.connect(self.ui.pushButton_findPeak, QtCore.SIGNAL('clicked()'),
                self.doFindPeak)

        self.connect(self.ui.pushButton_calUB, QtCore.SIGNAL('clicked()'),
                self.doCalUBMatrix)

        self.connect(self.ui.pushButton_acceptUB, QtCore.SIGNAL('clicked()'),
                self.doAcceptCalUB)

        self.connect(self.ui.pushButton_resetCalUB, QtCore.SIGNAL('clicked()'),
                self.doResetCalUB) 

        # Event handling for tab 'refine ub matrix'
        self.connect(self.ui.pushButton_addToRefine, QtCore.SIGNAL('clicked()'),
                self.doAddScanPtToRefineUB)

        # Validator


        # Declaration of class variable
        self._runID = None
        self._expID = None
        self._currPt = None
        self._xmlwsbasename = None
        
        # Some configuration
        self._homeSrcDir = os.getcwd()
        self._homeSaveDir = os.getcwd()

        # Control
        self._myControl = r4c.CWSCDReductionControl()

        # Initial setup
        # Tab 'Access'
        self.ui.lineEdit_url.setText('http://neutron.ornl.gov/user_data/hb3a/')
        self.ui.comboBox_mode.setCurrentIndex(0)
        self.ui.lineEdit_localSpiceDir.setEnabled(False)
        self.ui.pushButton_browseLocalDataDir.setEnabled(False)

        return

    #---------------------------------------------------------------------------
    # Event handling methods
    #---------------------------------------------------------------------------
    def change_data_access_mode(self):
        """ Change data access mode between downloading from server and local
        :return:
        """
        new_mode = str(self.ui.comboBox_mode.currentText())
        print '[DB] New Mode = ', new_mode
        if new_mode.startswith('Local') is True:
            self.ui.lineEdit_localSpiceDir.setEnabled(True)
            self.ui.pushButton_browseLocalDataDir.setEnabled(True)
            self.ui.lineEdit_url.setEnabled(False)
        else:
            self.ui.lineEdit_localSpiceDir.setEnabled(False)
            self.ui.pushButton_browseLocalDataDir.setEnabled(False)
            self.ui.lineEdit_url.setEnabled(True)

        return

    def doAcceptCalUB(self):
        """ Accept the calculated UB matrix
        """

        return

    def doAddScanPtToRefineUB(self):
        """ Add scan/pt numbers to the list of data points for refining ub matrix

        And the added scan number and pt numbers will be reflected in the (left sidebar)

        """
        raise NotImplementedError("ASAP")

        return

    def do_browse_local_spice_data(self):
        """ Browse local source SPICE data directory
        """
        src_spice_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Directory',
                                                                   self._homeSrcDir))
        self._homeSrcDir = src_spice_dir
        self.ui.lineEdit_localSpiceDir.setText(src_spice_dir)

        return


    def doBrowseSaveDir(self):
        """ Browse the local directory to save the data
        """
        targetdatadir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Directory', self._homeSaveDir))
        self._homeSaveDir = targetdatadir

        self.ui.lineEdit_dirSave.setText(targetdatadir)

        return 
    
    
    def doCalUBMatrix(self):
        """ Calculate UB matrix by 2 or 3 reflections
        """

        return

    def do_download_spice_data(self):
        """ Download SPICE data
        :return:
        """
        # Check scans to download
        scan_list_str = str(self.ui.lineEdit_downloadScans.text())
        if len(scan_list_str) > 0:
            # user specifies scans to download
            valid, scan_list = fcutil.parse_int_array(scan_list_str)
            if valid is False:
                error_message = scan_list
                self.pop_one_button_dialog(scan_list)
        else:
            # Get all scans
            scan_list = fcutil.get_scans_list(server_url, exp_no, return_list=True)
        self.pop_one_button_dialog('Going to download scans %s.' % str(scan_list))

        # Check location
        destination_dir = str(self.ui.lineEdit_localSrcDir.text())
        if os.path.exists(destination_dir) is False:
            self.pop_one_button_dialog('Destination directory %s cannot be found.' % destination_dir)
            return
        else:
            if os.access(destination_dir, os.W_OK) is False:
                self.pop_one_button_dialog('Destination directory %s is not writable.' % destination_dir)
                return
            else:
                self.pop_one_button_dialog('Spice files will be downloaded to %s.' % destination_dir)

        self._myControl.downloadSelectedDataSet(scan_list)

        return
    
    def doFindPeak(self):
        """ Find peak in a given scan/pt
        """
        scanNo = self._getInt(self.ui.lineEdit_scanNumber)
        ptNo = self._getInt(self.ui.lineEdit_ptNumber)

        self._myProject.findPeak(scanNo, ptNo)


        return

    def do_list_scans(self):
        """ List all scans available
        :return:
        """
        # Experiment number
        exp_no = int(self.ui.lineEdit_exp.text())

        access_mode = str(self.ui.comboBox_mode.currentText())
        if access_mode == 'Local':
            spice_dir = str(self.ui.lineEdit_localSpiceDir.text())
            message = fcutil.get_scans_list_local_disk(spice_dir, exp_no)
        else:
            url = str(self.ui.lineEdit_url.text())
            message = fcutil.get_scans_list(url, exp_no)

        self.pop_one_button_dialog(message)

        return

    def doLoad(self):
        """ Download and optinally load the data
        """
        # get experiment and run 
        expid = int(self.ui.lineEdit_exp.text())
        runid = int(self.ui.lineEdit_run.text())

        workdir = str(self.ui.lineEdit_dirSave.text())

        # load mode
        uselocalfile = self.ui.checkBox_dataLocal.status()
        if uselocalfile == 0:
            uselocalfile = False
        else:
            uselocalfile = True

        # determine operation mode
        if uselocalfile is True:
            source = str(self.ui.lineEdit_localSrcDir.text())
            mode = ['Copy', 'Reduce']
        else:
            source = str(self.ui.lineEdit_url.text())
            modestr = str(self.ui.comboBox_mode.currentText())
            mode = ['Download']
            if modestr.count('Reduce') == 1:
                mode.append('Reduce')

        self._loadData(source, workdir, expid, runid, mode)

        return


    def doPlotScanPt(self):
        """ Plot the Pt. 
        """
        # get measurement pt and the file number
        pt = int(self.ui.lineEdit_ptPlot.text())

        xmlwsname = self._xmlwsbasename + "_%d" % (pt)
        if self._xmlwkspdict.has_key(xmlwsname) is False:
            self._logError('Pt %d does not does not exist.' % (pt))
        xmlws = self._xmlwkspdict[xmlwsname]
        self._currPt = xmlws

        self._plotRawXMLWksp(self._currPt)

        return

        
    def doPlotPrevScanPt(self):
        """ Plot the Pt. 
        """
        # get measurement pt and the file number
        curindex = self._wkspNameList.index(self._currPt)
        prevwsname = self._wkspNameList[curindex-1]
        self._currPt = prevwsname

        self._plotRawXMLWksp(self._currPt)

        return


    def doPlotNextScanPt(self):
        """ Plot the Pt. 
        """
        # get measurement pt and the file number
        nextindex = self._wkspNameList.index(self._currPt) + 1
        if nextindex == len(self._wkspNameList):
            nextindex = 0
        nextwsname = self._wkspNameList[nextindex]
        self._currPt = nextwsname

        self._plotRawXMLWksp(self._currPt)

        return

    def doResetCalUB(self):
        """ Reset/reject the UB matrix calculation
        """

        return

    def do_test_url(self):
        """ Test whether the root URL provided specified is good
        """
        url = str(self.ui.lineEdit_url.text())

        url_is_good = fcutil.check_url(url)
        if url_is_good is True:
            self.pop_one_button_dialog("URL %s is valid." % url)
        else:
            self.pop_one_button_dialog("Unable to access %s.  Check internet access." % url)

        return


    def pop_one_button_dialog(self, message):
        """ Pop up a one-button dialog
        :param message:
        :return:
        """
        QtGui.QMessageBox.information(self, '4-circle Data Reduction', message)

        return

    #---------------------------------------------------------------------------
    # Private event handling methods
    #---------------------------------------------------------------------------
    def _loadData(self, source, targetdir, expid, runid, mode):
        """ Copy/download data to a directory and reduce them as an option
        Arguments:
         - source
         - targetdir
         - mode: 
        """
        basefilename =  "HB3A_exp%d_scan%0d.txt" % (expid, runid)
        localfilename = os.path.join(targetdir, basefilename)

        # load SPICE's run file
        if 'Download' in mode:
            # download from internet
            # generate the URL from 
            if source.endswith('/') is False:
                source = source+'/'
            spicerunfileurl = source + "HB3A_exp%d_scan%0d.txt" % (expid, runid)

            # download
            try:
                api.DownloadFile(Address=spicerunfileurl, Filename=localfilename)
            except Exception as e:
                return (False, str(e))

            # check file exist?
            if os.path.doesExist(localfilename) is False:
                return (False, "NO DOWNLOADED FILE")
            
        else:
            # copy from local disk
            # check whether the source and target directory are same
            source = os.path.absolutePath(source)
            targetdir = os.path.abosolutePath(targetdir)

            # copy file
            if source != targetdir:
                sourcefilename = os.path.join(source, basefilename)
                os.copyFile(sourcefilename, localfilename)

            # check file exist?
            if os.path.doesExist(localfilename) is False:
                return (False, "NO COPIED FILE")

        # ENDIFELSE

        # process SPICE's scan data file
        if 'Reduce' in mode:
            # load scan/run spice file
            spicetablews, infows = api.LoadSpiceAscii(Filename=localfilename, OutputWorkspace=spicetablewsname, 
                    RunInfoWorkspace=infowsname)

            # get Pt. data 
            ptlist = self._getPtList(spicetablews)

            self._xmlwkspdict = {} 
            for pt in ptlist:
                # generate xml file name
                basename = 'HB3A_exp%d_scan%04d_%04d.xml' % (expid, runid, pt)
                xmlfilename = os.path.join(targetdir, basename)
                if os.path.doesExist(xmlfilename) is False:
                    self._logError("File %s does not exist for exp %d scan %d pt %d" % (xmlfilename, expid, runid, pt))

                # load
                xmlwkspname = 'HB3A_e%d_s%d_m%d_raw' % (expid, runid, pt)
                xmlwksp = api.LoadSpiceXMLData(Filename=xmlfilename, OutputWorkspace=xmlwkspname)
                # FIXME - emit an signal?: for tree structure and log

                self._xmlwkspdict[pt] = xmlwksp
            # ENDFOR
        # ENDIF

        return
        

    def _plotRawXMLWksp(self, xmlws):
        """ Plot raw workspace from XML file for a measurement/pt.
        """
        # get data
        numspec = xmlws.getNumberHistograms()
        vecylist = []
        for iws in xrange(len(numspec)): 
            vecy = xmlws.readY(0)
            vecylist.append(vecy)
        # ENDFOR(iws)

        # plot 
        self._plot2DData(vecylist)

        return


