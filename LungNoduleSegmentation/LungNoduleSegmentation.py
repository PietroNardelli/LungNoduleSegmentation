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
    parametersCollapsibleButton.text = "IO Images"
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
    #########################  Prior Prob Of Malignancy  ##############################
    ################################################################################### 

    priorProbCollapsibleButton = ctk.ctkCollapsibleButton()
    priorProbCollapsibleButton.text = "Prior Probability of Malignancy"
    self.layout.addWidget(priorProbCollapsibleButton)

    # Layout within the dummy collapsible button
    priorFormLayout = qt.QFormLayout(priorProbCollapsibleButton)

    PriorProbBox = qt.QHBoxLayout()
    priorFormLayout.addRow(PriorProbBox)

    self.priorComboBox = qt.QComboBox()
    self.priorComboBox.setLayoutDirection(1)
    self.priorComboBox.setFixedWidth(50)

    priorLabel = qt.QLabel()
    priorLabel.setText('Select Prior Probability (%): ')
    priorLabel.setFixedWidth(150)

    PriorSelectionBox = qt.QFrame()
    PriorSelectionBox.setFixedWidth(210)
    priorHBox = qt.QHBoxLayout()
    PriorSelectionBox.setLayout(priorHBox)

    PriorSelectionBox.layout().addWidget(priorLabel,0,4)
    PriorSelectionBox.layout().addWidget(self.priorComboBox)

    PriorProbBox.addWidget(PriorSelectionBox,0,4)

    ###################################################################################
    ##############################  Patient Info Area  ################################
    ################################################################################### 

    inputInfoCollapsibleButton = ctk.ctkCollapsibleButton()
    inputInfoCollapsibleButton.text = "Patient Information"
    self.layout.addWidget(inputInfoCollapsibleButton)

    # Layout within the dummy collapsible button
    infoFormLayout = qt.QFormLayout(inputInfoCollapsibleButton)

    PatientInfoSelectionBoxOne = qt.QHBoxLayout()
    infoFormLayout.addRow(PatientInfoSelectionBoxOne)

    PatientInfoSelectionBoxTwo = qt.QHBoxLayout()
    infoFormLayout.addRow(PatientInfoSelectionBoxTwo)

    ###################################################################################
    ################################  Age Selector  ###################################
    ################################################################################### 

    self.age = qt.QLineEdit()
    self.age.setReadOnly(1)
    self.age.setFixedWidth(40)

    ageLabel = qt.QLabel()
    ageLabel.setText('Patient Age: ')
    ageLabel.setFixedWidth(60)

    AgeSelectionBox = qt.QFrame()
    AgeSelectionBox.setFixedWidth(120)
    ageHBox = qt.QHBoxLayout()
    AgeSelectionBox.setLayout(ageHBox)

    AgeSelectionBox.layout().addWidget(ageLabel,0,4)
    AgeSelectionBox.layout().addWidget(self.age)

    PatientInfoSelectionBoxOne.addWidget(AgeSelectionBox,0,4)

    ###################################################################################
    ##########################  Smoking Status Selector  ##############################
    ################################################################################### 

    self.smokingComboBox = qt.QComboBox()
    self.smokingComboBox.setLayoutDirection(1)
    self.smokingComboBox.setFixedWidth(90)

    smokeLabel = qt.QLabel()
    smokeLabel.setText('Smoker (Pk-Yrs): ')
    smokeLabel.setFixedWidth(80)

    SmokeSelectionBox = qt.QFrame()
    SmokeSelectionBox.setFixedWidth(190)
    smokeHBox = qt.QHBoxLayout()
    SmokeSelectionBox.setLayout(smokeHBox)

    SmokeSelectionBox.layout().addWidget(smokeLabel,0,4)
    SmokeSelectionBox.layout().addWidget(self.smokingComboBox)

    PatientInfoSelectionBoxOne.addWidget(SmokeSelectionBox,0,4)

    ###################################################################################
    ############################  Hemoptysis Selector  ################################
    ################################################################################### 

    self.hemoptysisComboBox = qt.QComboBox()
    self.hemoptysisComboBox.setLayoutDirection(1)
    self.hemoptysisComboBox.setFixedWidth(80)

    hemoptysisLabel = qt.QLabel()
    hemoptysisLabel.setText('Hemoptysis: ')
    hemoptysisLabel.setFixedWidth(70)

    HemoptysisSelectionBox = qt.QFrame()
    HemoptysisSelectionBox.setFixedWidth(170)
    hemoptysisHBox = qt.QHBoxLayout()
    HemoptysisSelectionBox.setLayout(hemoptysisHBox)

    HemoptysisSelectionBox.layout().addWidget(hemoptysisLabel,0,4)
    HemoptysisSelectionBox.layout().addWidget(self.hemoptysisComboBox)

    PatientInfoSelectionBoxTwo.addWidget(HemoptysisSelectionBox,0,4)

    ###################################################################################
    ############################  Hx Prev Selector  ################################
    ################################################################################### 

    self.prevMalComboBox = qt.QComboBox()
    self.prevMalComboBox.setLayoutDirection(1)
    self.prevMalComboBox.setFixedWidth(80)

    prevMalLabel = qt.QLabel()
    prevMalLabel.setText('Hx Prev Malig: ')
    prevMalLabel.setFixedWidth(70)

    PrevMalSelectionBox = qt.QFrame()
    PrevMalSelectionBox.setFixedWidth(170)
    prevMalHBox = qt.QHBoxLayout()
    PrevMalSelectionBox.setLayout(prevMalHBox)

    PrevMalSelectionBox.layout().addWidget(prevMalLabel,0,4)
    PrevMalSelectionBox.layout().addWidget(self.prevMalComboBox)

    PatientInfoSelectionBoxTwo.addWidget(PrevMalSelectionBox,0,4)

    ###################################################################################
    #######################  Button to trigger segmentation  ##########################
    ###################################################################################    

    self.LungNoduleSegmentationButton = qt.QPushButton("Start Nodule Evaluation")
    self.LungNoduleSegmentationButton.toolTip = "Run the algorithm to segment the lung"
    self.LungNoduleSegmentationButton.setFixedSize(200,50)
    self.LungNoduleSegmentationButton.enabled = False

    boxLayout = qt.QVBoxLayout()

    #IOFormLayout.addRow(boxLayout)
    #boxLayout.addWidget(self.LungNoduleSegmentationButton,0,4)
    self.layout.addWidget(self.LungNoduleSegmentationButton,0,4)

    ###################################################################################
    #################### Metric Results + Perc. of Malignancy #########################
    ###################################################################################

    resultCollapsibleButton = ctk.ctkCollapsibleButton()
    resultCollapsibleButton.text = "Results Area"
    self.layout.addWidget(resultCollapsibleButton)
    
    # Layout within the dummy collapsible button
    resultsFormLayout = qt.QFormLayout(resultCollapsibleButton)

    FirstRowMetricHBox = qt.QHBoxLayout()
    resultsFormLayout.addRow(FirstRowMetricHBox)
    
    # Nodule Size
    self.size = qt.QLineEdit()
    self.size.setReadOnly(1)
    self.size.setFixedWidth(80)

    sizeLabel = qt.QLabel()
    sizeLabel.setText('Nodule Diameter: ')
    sizeLabel.setFixedWidth(90)

    SizeSelectionBox = qt.QFrame()
    SizeSelectionBox.setFixedWidth(180)
    sHBox = qt.QHBoxLayout()
    SizeSelectionBox.setLayout(sHBox)

    SizeSelectionBox.layout().addWidget(sizeLabel,0,4)
    SizeSelectionBox.layout().addWidget(self.size)

    FirstRowMetricHBox.addWidget(SizeSelectionBox,0,4)

    # Nodule Location
    self.location = qt.QLineEdit()
    self.location.setReadOnly(1)
    self.location.setFixedWidth(180)

    locationLabel = qt.QLabel()
    locationLabel.setText('Nodule Location: ')
    locationLabel.setFixedWidth(80)

    LocationSelectionBox = qt.QFrame()
    LocationSelectionBox.setFixedWidth(270)
    lHBox = qt.QHBoxLayout()
    LocationSelectionBox.setLayout(lHBox)

    LocationSelectionBox.layout().addWidget(locationLabel,0,4)
    LocationSelectionBox.layout().addWidget(self.location)

    FirstRowMetricHBox.addWidget(LocationSelectionBox,0,4)

    # Nodule Cavity Wall Thickness
    self.thickness = qt.QLineEdit()
    self.thickness.setReadOnly(1)
    self.thickness.setFixedWidth(80)

    thickLabel = qt.QLabel()
    thickLabel.setText('Cavity Wall Thickness: ')
    thickLabel.setFixedWidth(110)

    ThicknessSelectionBox = qt.QFrame()
    ThicknessSelectionBox.setFixedWidth(200)
    tHBox = qt.QHBoxLayout()
    ThicknessSelectionBox.setLayout(tHBox)

    ThicknessSelectionBox.layout().addWidget(thickLabel,0,4)
    ThicknessSelectionBox.layout().addWidget(self.thickness)

    FirstRowMetricHBox.addWidget(ThicknessSelectionBox,0,4)

    # Spiculation
    SecondRowMetricHBox = qt.QHBoxLayout()
    resultsFormLayout.addRow(SecondRowMetricHBox)
    
    self.spiculation = qt.QLineEdit()
    self.spiculation.setReadOnly(1)
    self.spiculation.setFixedWidth(130)

    spicLabel = qt.QLabel()
    spicLabel.setText('Nodule Smoothness: ')
    spicLabel.setFixedWidth(100)

    SpiculationSelectionBox = qt.QFrame()
    SpiculationSelectionBox.setFixedWidth(240)
    spicHBox = qt.QHBoxLayout()
    SpiculationSelectionBox.setLayout(spicHBox)

    SpiculationSelectionBox.layout().addWidget(spicLabel,0,4)
    SpiculationSelectionBox.layout().addWidget(self.spiculation)
    
    SecondRowMetricHBox.addWidget(SpiculationSelectionBox,0,4)

    # Calcification
    self.calcification = qt.QLineEdit()
    self.calcification.setReadOnly(1)
    self.calcification.setFixedWidth(150)

    calcLabel = qt.QLabel()
    calcLabel.setText('Calcification Pattern: ')
    calcLabel.setFixedWidth(100)

    CalcificationSelectionBox = qt.QFrame()
    CalcificationSelectionBox.setFixedWidth(260)
    calcHBox = qt.QHBoxLayout()
    CalcificationSelectionBox.setLayout(calcHBox)

    CalcificationSelectionBox.layout().addWidget(calcLabel,0,4)
    CalcificationSelectionBox.layout().addWidget(self.calcification)

    SecondRowMetricHBox.addWidget(CalcificationSelectionBox,0,4)
    
    # Percentage of malignancy
    self.percentageOfMalignancy = qt.QLineEdit()
    self.percentageOfMalignancy.setReadOnly(1)
    self.percentageOfMalignancy.setFixedWidth(50)

    malignancyLabel = qt.QLabel()
    malignancyLabel.setText('Percentage Of Malignancy: ')
    malignancyLabel.setStyleSheet("background-color: rgb(255,255,102)")
    malignancyLabel.setFixedWidth(130)

    MalignancySelectionBox = qt.QFrame()
    MalignancySelectionBox.setFixedWidth(210)
    mHBox = qt.QHBoxLayout()
    MalignancySelectionBox.setLayout(mHBox)

    MalignancySelectionBox.layout().addWidget(malignancyLabel,0,4)
    MalignancySelectionBox.layout().addWidget(self.percentageOfMalignancy)

    resultsHBoxLayout = qt.QHBoxLayout()
    resultsHBoxLayout.addWidget(MalignancySelectionBox,0,4)

    resultsFormLayout.addRow(resultsHBoxLayout)


    self.RecomputeButton = qt.QPushButton("Recalculate Probability")
    self.RecomputeButton.toolTip = "Recalculate probabilty of malignancy with new parameters"
    self.RecomputeButton.setFixedSize(200,50)
    
    #recomputeBoxLayout = qt.QVBoxLayout()

    #IOFormLayout.addRow(boxLayout)
    #boxLayout.addWidget(self.LungNoduleSegmentationButton,0,4)
    self.layout.addWidget(self.RecomputeButton,0,4)
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
    self.updateGUI()

  def updateGUI(self):
    self.smokingComboBox.addItem('Not Known')
    self.smokingComboBox.addItem('Non-Smoker')
    self.smokingComboBox.addItem('< 30 Pk-Yrs ')
    self.smokingComboBox.addItem('30-39 Pk-Yrs')
    self.smokingComboBox.addItem('> 40 Pk-Yrs ')

    self.smokingComboBox.setCurrentIndex(0)
    self.smokingComboBox.setEditable(0)

    self.hemoptysisComboBox.addItem('Absent')
    self.hemoptysisComboBox.addItem('Present')
    
    self.hemoptysisComboBox.setCurrentIndex(0)
    self.hemoptysisComboBox.setEditable(0)

    self.prevMalComboBox.addItem('Absent')
    self.prevMalComboBox.addItem('Present')
    
    self.prevMalComboBox.setCurrentIndex(0)
    self.prevMalComboBox.setEditable(0)

    for i in xrange(1,100):
      self.priorComboBox.addItem(i)
      self.priorComboBox.setCurrentIndex(49)
    
  def onSelect(self):
    if self.inputSelector.currentNode():
      # Read patient's age from DICOM metadata
      inputVolume = self.inputSelector.currentNode()
      inputName = inputVolume.GetName()

      inputNode = slicer.util.getNode(inputName)
      instUIDs = inputNode.GetAttribute('DICOM.instanceUIDs').split()
      fileName = slicer.dicomDatabase.fileForInstance(instUIDs[0])
      patientAge = slicer.dicomDatabase.fileValue(fileName,'0010,1010')
      patientAge = patientAge[1:3]
      age_str = str(patientAge)
      self.age.setText(age_str)
      
    if self.inputSelector.currentNode() and self.fiducialListSelector.currentNode():
      self.LungNoduleSegmentationButton.enabled = True

  def onStartSegmentationButton(self):
    inputVolume = self.inputSelector.currentNode()
    inputName = inputVolume.GetName()
    
    inputNode = slicer.util.getNode(inputName)
    instUIDs = inputNode.GetAttribute('DICOM.instanceUIDs').split()
    fileName = slicer.dicomDatabase.fileForInstance(instUIDs[0])
    patientAge = slicer.dicomDatabase.fileValue(fileName,'0010,1010')
    patientAge=patientAge[1:3]
    age_str = str(patientAge)
    self.age.setText(age_str)
    
    # Analyze histogram of the image to work out the second peak intensity value
    start = timeit.default_timer()

    self.LungNoduleSegmentationButton.enabled = False   
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

    NoduleSegmentationNode = slicer.cli.run( LungNoduleSegmentationModule, None, parameters, wait_for_completion = True )

    stop = timeit.default_timer()
    print "Nodule segmentation: ", stop - start, "sec"     

    self.FissuresNode = slicer.mrmlScene.CreateNodeByClass(nodeType)
    self.FissuresNode.SetLabelMap(1)
    self.FissuresNode.SetScene(slicer.mrmlScene)
    self.FissuresNode.SetName(slicer.mrmlScene.GetUniqueNameByString('FissuresVolume'))
    slicer.mrmlScene.AddNode(self.FissuresNode)

    FissuresSegmentationModule = slicer.modules.lungfissuresegmentationcli

    params = {
       "InputVolume": inputVolume,
       "InputLungLabel": self.lungLabelNode, 
       "OutputVolume": self.FissuresNode,
    }

    slicer.cli.run( FissuresSegmentationModule, None, params, wait_for_completion = True )

    self.LobesNode = slicer.mrmlScene.CreateNodeByClass(nodeType)
    self.LobesNode.SetLabelMap(1)
    self.LobesNode.SetScene(slicer.mrmlScene)
    self.LobesNode.SetName(slicer.mrmlScene.GetUniqueNameByString('LobesVolume'))
    slicer.mrmlScene.AddNode(self.LobesNode)

    LobesSegmentationModule = slicer.modules.lunglobesegmentationcli

    params = {
       "InputVolume": inputVolume,
       "LabelMapVolume": self.lungLabelNode, 
       "FissuresVolume": self.FissuresNode,
       "OutputVolume": self.LobesNode
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

    # Evaluate malignancy of the nodule

    nodulePosition = NoduleSegmentationNode.GetParameterDefault(2,0)
    noduleSize = NoduleSegmentationNode.GetParameterDefault(2,1)
    noduleRoundness = NoduleSegmentationNode.GetParameterDefault(2,2)
    noduleCavityWallThickness = NoduleSegmentationNode.GetParameterDefault(2,3) 
    noduleCalcificationPattern = NoduleSegmentationNode.GetParameterDefault(2,4)

    size_str = str(noduleSize) + ' mm'
    self.size.setText(size_str)

    smoothness_str = 'Smooth'
    noduleRoundness = float(noduleRoundness)
    if noduleRoundness < 0.7:
      smoothness_str = 'Spiculated'

    self.spiculation.setText(smoothness_str)

    thickness_str = 'Not Cavitated'
    noduleCavityWallThickness = float(noduleCavityWallThickness)
    if (noduleCavityWallThickness >= 0.1 and noduleCavityWallThickness < 100.0):
      print noduleCavityWallThickness
      thickness_str = str(noduleCavityWallThickness) + ' mm'
    self.thickness.setText(thickness_str)

    lobesImageData = self.LobesNode.GetImageData()

    nodulePosition = eval(nodulePosition)

    lobeValue = lobesImageData.GetScalarComponentAsDouble(int(nodulePosition[0]),int(nodulePosition[1]),int(nodulePosition[2]), 0)
    position_str = 'Right Upper Lobe'
    if lobeValue == 5:
      position_str = 'Right Middle Lobe'
    elif lobeValue == 6:
      position_str = 'Right Lower Lobe'
    elif lobeValue == 7:
      position_str = 'Left Upper Lobe'
    elif lobeValue == 8:
      position_str = 'Left Lower Lobe'

    self.location.setText(position_str)

    calc_str = 'Benign Pattern'
    if noduleCalcificationPattern == 'false':
      calc_str = 'No Calcification'

    self.calcification.setText(calc_str)
    
    string_percentage = str(lobeValue) + '%'    
       
    self.percentageOfMalignancy.setText(string_percentage)
 
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

