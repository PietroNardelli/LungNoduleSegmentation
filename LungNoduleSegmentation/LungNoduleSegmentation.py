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
    Scripted module bundled in an extension for Lung Nodule Segmentation.
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
    self.fiducialListSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    self.LungNoduleSegmentationButton.connect('clicked(bool)', self.onStartSegmentationButton)
  
    #
    # Add Vertical Spacer
    #
    self.layout.addStretch(1)

  def onSelect(self):

    if self.inputSelector.currentNode() and self.fiducialListSelector.currentNode():
      self.LungNoduleSegmentationButton.enabled = True

  def onStartSegmentationButton(self):
    # Analyze histogram of the image to work out the second peak intensity value
    start = timeit.default_timer()

    self.LungNoduleSegmentationButton.enabled = False   
    inputVolume = self.inputSelector.currentNode()
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
    if threshold < -800:
      maskedImage = inputNodeArray[inputNodeArray >= -1000]

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
    
    stop = timeit.default_timer()
    print "Histogram anaylis: ", stop - start, "sec" 

    # Segment the lung surface

    start = timeit.default_timer()

    nodeType = 'vtkMRMLScalarVolumeNode'
    self.lungLabelNode = slicer.mrmlScene.CreateNodeByClass(nodeType)
    self.lungLabelNode.SetLabelMap(1)
    self.lungLabelNode.SetScene(slicer.mrmlScene)
    self.lungLabelNode.SetName(slicer.mrmlScene.GetUniqueNameByString('LungLabel'))
    slicer.mrmlScene.AddNode(self.lungLabelNode)

    LungSegmentationModule = slicer.modules.lungsegmentationcli
    lungLabelColor = 1

    parameters = {
       "InputVolume": inputVolume,
       "OutputVolume": self.lungLabelNode,
       "Lower": -1024,
       "Upper": threshold,
       "lungColor": lungLabelColor,
    }

    slicer.cli.run( LungSegmentationModule, None, parameters, wait_for_completion = True )

    stop = timeit.default_timer()
    print "Lung segmentation: ", stop - start, "sec" 

    # Segment the selected nodule
    start = timeit.default_timer()

    nodeType = 'vtkMRMLScalarVolumeNode'
    self.noduleLabelNode = slicer.mrmlScene.CreateNodeByClass(nodeType)
    self.noduleLabelNode.SetLabelMap(1)
    self.noduleLabelNode.SetScene(slicer.mrmlScene)
    self.noduleLabelNode.SetName(slicer.mrmlScene.GetUniqueNameByString('NoduleLabel'))
    slicer.mrmlScene.AddNode(self.noduleLabelNode)

    LungNoduleSegmentationModule = slicer.modules.lungnodulesegmentationcli
    labelColor = 6

    parameters = {
       "InputVolume": inputVolume,
       "InputLungLabel": self.lungLabelNode, 
       "OutputVolume": self.noduleLabelNode,
       "seed": seedPoint,
       "noduleColor": labelColor,
    }

    slicer.cli.run( LungNoduleSegmentationModule, None, parameters, wait_for_completion = True )

    stop = timeit.default_timer()
    print "Nodule segmentation: ", stop - start, "sec" 

    self.LobesNode = slicer.mrmlScene.CreateNodeByClass(nodeType)
    self.LobesNode.SetLabelMap(1)
    self.LobesNode.SetScene(slicer.mrmlScene)
    self.LobesNode.SetName(slicer.mrmlScene.GetUniqueNameByString('LobesVolume'))
    slicer.mrmlScene.AddNode(self.LobesNode)

    LobesSegmentationModule = slicer.modules.lobesegmentationcli

    params = {
       "InputVolume": inputVolume,
       "InputLungLabel": self.lungLabelNode, 
       "OutputVolume": self.LobesNode,
    }

    slicer.cli.run( LobesSegmentationModule, None, params, wait_for_completion = True )

    # Create 3D model of the segment nodule 
    start = timeit.default_timer()

    modelHierarchyCollection = slicer.mrmlScene.GetNodesByName('NoduleModelHierarchy')
    if( modelHierarchyCollection.GetNumberOfItems() >= 1 ):
      modelHierarchy = modelHierarchyCollection.GetItemAsObject(0)
    else:
      modelHierarchy = slicer.vtkMRMLModelHierarchyNode()
      modelHierarchy.SetName(slicer.mrmlScene.GetUniqueNameByString('NoduleModelHierarchy'))
      slicer.mrmlScene.AddNode(modelHierarchy)
    
    parameters = {}
    parameters["InputVolume"] = self.noduleLabelNode.GetID()
    parameters["ModelSceneFile"] = modelHierarchy.GetID()
    parameters["Name"] = 'NoduleModel'
    parameters["Smooth"] = 20
    parameters["Decimate"] = 0.10
    
    modelMaker = slicer.modules.modelmaker
    slicer.cli.run( modelMaker, None, parameters, wait_for_completion = True )
  
    lm = slicer.app.layoutManager()
    threeDView = lm.threeDWidget( 0 ).threeDView()
    threeDView.resetFocalPoint()
    threeDView.lookFromViewAxis(ctk.ctkAxesWidget().Anterior)

    stop = timeit.default_timer()
    print "Model maker: ", stop - start, "sec" 

    self.LungNoduleSegmentationButton.enabled = True

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

    return True

