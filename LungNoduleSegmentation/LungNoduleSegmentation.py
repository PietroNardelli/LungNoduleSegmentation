import os
import unittest
from __main__ import vtk, qt, ctk, slicer
import numpy
import timeit

#
# LungNoduleSegmentation
#

class LungNoduleSegmentation:
  def __init__(self, parent):
    parent.title = "Lung Nodule Segmentation" # TODO make this more human readable by adding spaces
    parent.categories = ["Segmentation"]
    parent.dependencies = []
    parent.contributors = ["Pietro Nardelli (University College Cork)"] 
    parent.helpText = """
    Scripted module bundled in an extension for Lung Segmentation.
    """
    parent.acknowledgementText = """
    This file was originally developed by Pietro Nardelli, University College Cork.
""" 
    self.parent = parent

#
# qLungNoduleSegmentationWidget
#
class LungNoduleSegmentationWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent

    self.layout = self.parent.layout()
    #self.updateGUI()

    if not parent:
      self.setup()
      self.parent.show()
      #self.updateGUI()

  def cleanup(self):
    pass

  def setup(self):
    #
    # Reload and Test area
    #
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Reload Section"
    self.layout.addWidget(reloadCollapsibleButton)
    self.layout.setSpacing(6)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "LungNoduleSegmentation Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Lung Segmentation"
    self.layout.addWidget(parametersCollapsibleButton)
    self.layout.setSpacing(20)
    # Layout within the dummy collapsible button
    IOFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    ###################################################################################
    #############################  Input Data Selector  ###############################
    ###################################################################################
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = True
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the chest CT data to segment." )
    IOFormLayout.addRow("Input CT Data: ", self.inputSelector)

    ###################################################################################
    ###############################  Label Selector  ##################################
    ###################################################################################
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 1 )
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = False
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output (label) of the algorithm." )
    IOFormLayout.addRow("Lung Label: ", self.outputSelector)

    ###################################################################################
    ###########################  Fiducial Point Selector  #############################
    ################################################################################### 

    self.fiducialListSelector = slicer.qMRMLNodeComboBox()
    self.fiducialListSelector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.fiducialListSelector.selectNodeUponCreation = True
    self.fiducialListSelector.addEnabled = False
    self.fiducialListSelector.removeEnabled = True
    self.fiducialListSelector.noneEnabled = True
    self.fiducialListSelector.showHidden = False
    self.fiducialListSelector.showChildNodeTypes = False
    self.fiducialListSelector.setMRMLScene( slicer.mrmlScene )
    self.fiducialListSelector.setToolTip( "Pick the list containing the ROI fiducial" )
    IOFormLayout.addRow("Seed Point: ", self.fiducialListSelector)

    ###################################################################################
    #######################  Button to trigger segmentation  ##########################
    ###################################################################################    

    self.LungNoduleSegmentationButton = qt.QPushButton("Start Lung Segmentation")
    self.LungNoduleSegmentationButton.toolTip = "Run the algorithm to segment the lung"
    self.LungNoduleSegmentationButton.setFixedSize(300,50)
    self.LungNoduleSegmentationButton.enabled = False

    boxLayout = qt.QVBoxLayout()

    IOFormLayout.addRow(boxLayout)
    boxLayout.addWidget(self.LungNoduleSegmentationButton,0,4)

    ########################################################################################
    ################################ Create Connections ####################################
    ########################################################################################
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    self.LungNoduleSegmentationButton.connect('clicked(bool)', self.onStartSegmentationButton)
  
    #
    # Add Vertical Spacer
    #
    self.layout.addStretch(1)

  def onSelect(self):

    if self.inputSelector.currentNode() and self.outputSelector.currentNode():
      self.LungNoduleSegmentationButton.enabled = True

  def onStartSegmentationButton(self):
    inputVolume = self.inputSelector.currentNode()
    outputVolume = self.outputSelector.currentNode()
    seedPoint = self.fiducialListSelector.currentNode()

    inputNodeArray = slicer.util.array(inputVolume.GetID())
    maskedImage = inputNodeArray[inputNodeArray >= -1024]

    edges=[] 
    I_max = maskedImage.max()
    I_min = maskedImage.min()

    for val in range(I_min,I_max+2): 
      edges.append(val)

    I_max =float(I_max)
    I_min =float(I_min)

    maskedImage = maskedImage.astype(float)
   
    histogram,bins = numpy.histogram(maskedImage,bins=edges)

    mx, mn = self.peakdetect(histogram)

    threshold = (-1024+bins[mx[0][0]]) / 2
    print threshold

    LungNoduleSegmentationModule = slicer.modules.lungnodulesegmentationcli

    parameters = {
       "InputVolume": inputVolume,
       "OutputVolume": outputVolume,
       "seed": seedPoint,
       "Lower": -1024,
       "Upper": threshold,
       "noduleColor": 6,
    }

    slicer.cli.run( LungNoduleSegmentationModule,None,parameters,wait_for_completion=False )

  def datacheck_peakdetect(self, x_axis, y_axis):
    if x_axis is None:
        x_axis = range(len(y_axis))
    
    if len(y_axis) != len(x_axis):
        raise (ValueError, 
                'Input vectors y_axis and x_axis must have same length')
    
    #needs to be a numpy array
    y_axis = numpy.array(y_axis)
    x_axis = numpy.array(x_axis)
    return x_axis, y_axis

  def peakdetect(self, y_axis, x_axis = None, lookahead = 300, delta=0):
    """
    Converted from/based on a MATLAB script at: 
    http://billauer.co.il/peakdet.html
    
    function for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- (optional) A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function
    
    return -- two lists [max_peaks, min_peaks] containing the positive and
        negative peaks respectively. Each cell of the lists contains a tupple
        of: (position, peak_value) 
        to get the average peak value do: np.mean(max_peaks, 0)[1] on the
        results to unpack one of the lists into x, y coordinates do: 
        x, y = zip(*tab)
    """
    max_peaks = []
    min_peaks = []
    dump = []   #Used to pop the first hit which almost always is false
       
    # check input data
    x_axis, y_axis = self.datacheck_peakdetect(x_axis, y_axis)
    # store data length for later use
    length = len(y_axis)
    
    
    #perform some checks
    if lookahead < 1:
        raise ValueError, "Lookahead must be '1' or above in value"
    if not (numpy.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"
    
    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = numpy.Inf, -numpy.Inf
    
    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], 
                                        y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x
        
        ####look for max####
        if y < mx-delta and mx != numpy.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                mx = numpy.Inf
                mn = numpy.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue
            #else:  #slows shit down this does
            #    mx = ahead
            #    mxpos = x_axis[np.where(y_axis[index:index+lookahead]==mx)]
        
        ####look for min####
        if y > mn+delta and mn != -numpy.Inf:
            #Minima peak candidate found 
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                min_peaks.append([mnpos, mn])
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -numpy.Inf
                mx = -numpy.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
            #else:  #slows shit down this does
            #    mn = ahead
            #    mnpos = x_axis[np.where(y_axis[index:index+lookahead]==mn)]
    
    
    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            max_peaks.pop(0)
        else:
            min_peaks.pop(0)
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass
        
    return [max_peaks, min_peaks]


  def onReload(self,moduleName="LungNoduleSegmentation"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    import os
    import unittest
    from __main__ import vtk, qt, ctk, slicer, numpy
    import imp, sys

    widgetName = moduleName + "Widget"

    # reload the source code
    # - set source file path
    # - load the module to the global space
    filePath = eval('slicer.modules.%s.path' % moduleName.lower())
    p = os.path.dirname(filePath)
    if not sys.path.__contains__(p):
      sys.path.insert(0,p)
    fp = open(filePath, "r")
    globals()[moduleName] = imp.load_module(
        moduleName, fp, filePath, ('.py', 'r', imp.PY_SOURCE))
    fp.close()

    # rebuild the widget
    # - find and hide the existing widget
    # - create a new widget in the existing parent
    parent = slicer.util.findChildren(name='%s Reload' % moduleName)[0].parent().parent()
    for child in parent.children():
      try:
        child.hide()
      except AttributeError:
        pass
    # Remove spacer items
    item = parent.layout().itemAt(0)
    while item:
      parent.layout().removeItem(item)
      item = parent.layout().itemAt(0)

    # delete the old widget instance
    if hasattr(globals()['slicer'].modules, widgetName):
      getattr(globals()['slicer'].modules, widgetName).cleanup()

    # create new widget inside existing parent
    globals()[widgetName.lower()] = eval(
        'globals()["%s"].%s(parent)' % (moduleName, widgetName))
    globals()[widgetName.lower()].setup()
    setattr(globals()['slicer'].modules, widgetName, globals()[widgetName.lower()])

