/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Modules/CLI/Threshold.cxx $
  Language:  C++
  Date:      $Date: 2006-12-20 16:00:24 -0500 (Wed, 20 Dec 2006) $
  Version:   $Revision: 1892 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// ITK includes
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMedianImageFilter.h"

#include "itkImageDuplicator.h"

#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryDilateImageFilter.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"

#include "itkRegionOfInterestImageFilter.h"

#include "itkConnectedThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkPasteImageFilter.h"

#include "itkNotImageFilter.h"
#include "itkAndImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "itkImageMomentsCalculator.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "LungSegmentationCLICLP.h"


const unsigned int Dim = 3;

typedef signed short 	InputPixelType;
typedef unsigned short 	OutputPixelType;
typedef unsigned int 	LabelType;
typedef float		MSOutputPixelType;

typedef itk::Image< InputPixelType,  Dim > InputImageType;
typedef itk::Image< OutputPixelType, Dim > OutputImageType;

OutputImageType::Pointer AirwaySegmentation( InputImageType::Pointer image, InputImageType::IndexType indexFiducialSlice, 
                                             InputImageType::IndexType index, int labelColor = 1 )
{
    OutputImageType::Pointer airways 		= OutputImageType::New(); 
    OutputImageType::Pointer airwaysPrev 	= OutputImageType::New();
	
    /*airways->SetRegions( image->GetRequestedRegion() );                                   
    airways->SetBufferedRegion( image->GetBufferedRegion() );
    airways->SetLargestPossibleRegion( image->GetLargestPossibleRegion() );
    airways->CopyInformation( image );
    airways->Allocate();*/
  		
    /** AIRWAYS SEGMENTATION PIPELINE */
    typedef itk::ConnectedThresholdImageFilter< InputImageType, InputImageType > ConnectedFilterType; 
    ConnectedFilterType::Pointer thresholdConnected = ConnectedFilterType::New(); 
	
    thresholdConnected->SetInput( image );	
    thresholdConnected->SetReplaceValue( labelColor ); 
    thresholdConnected->AddSeed( index );   
   	
    // Starting upper threshold value
    InputPixelType UpperThreshold = -900;
    //thresholdConnected->SetLower( -2000 );
    thresholdConnected->SetUpper( UpperThreshold );          

    typedef itk::CastImageFilter<InputImageType, OutputImageType> CastingFilterType;  
    CastingFilterType::Pointer  caster = CastingFilterType::New();		
	
    caster->SetInput( thresholdConnected->GetOutput() );  
    caster->Update();
    airways = caster->GetOutput();
    
    /** COMPUTING THE LABEL SIZES */ 	                  
    OutputImageType::Pointer airwaysAxialCopy = OutputImageType::New(); 
    OutputImageType::Pointer airwaysCoronalCopy = OutputImageType::New(); 

    typedef itk::ImageDuplicator<OutputImageType> DuplicatorFilterType;

    DuplicatorFilterType::Pointer duplicatorFilter = DuplicatorFilterType::New();
    duplicatorFilter->SetInputImage(airways);
    duplicatorFilter->Update();

    // Extracting the axial slice containing the airways fiducial point
    OutputImageType::SizeType  oneAxialSliceSize;
    InputImageType::IndexType  indexAxialSlice = indexFiducialSlice;
  	
    oneAxialSliceSize[0] = airways->GetLargestPossibleRegion().GetSize(0);
    oneAxialSliceSize[1] = airways->GetLargestPossibleRegion().GetSize(1);
    unsigned int diff = airways->GetLargestPossibleRegion().GetSize(2)-indexAxialSlice[2];
    if( airways->GetLargestPossibleRegion().GetSize(2) > 40 &&
        indexAxialSlice[2] >= 20 &&
        diff >= 20 )
    {
        oneAxialSliceSize[2] = 40;
        indexAxialSlice[2]  -= 20;
    }
    else if( airways->GetLargestPossibleRegion().GetSize(2) > 40 &&
             indexAxialSlice[2] >= 20 &&
             diff < 20 )
    {
        oneAxialSliceSize[2] = 40;
        indexAxialSlice[2]   = airways->GetLargestPossibleRegion().GetSize(2) - 40;
    }
    else if( airways->GetLargestPossibleRegion().GetSize(2) > 40 && indexAxialSlice[2] < 20 )
    {
        oneAxialSliceSize[2] = 40;
        indexAxialSlice  [2] = 0;
    }
    else if( airways->GetLargestPossibleRegion().GetSize(2) <= 40 )
    {
        oneAxialSliceSize[2] = airways->GetLargestPossibleRegion().GetSize(2);
        indexAxialSlice  [2] = 0;
    }
	
    OutputImageType::RegionType axialSlice;
 
    typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > ROIFilterType;
    ROIFilterType::Pointer axialAirwaysFilter = ROIFilterType::New();

    typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > ShapeLabelType;	
    ShapeLabelType::Pointer axialLabelSizeFilter = ShapeLabelType::New();

    axialLabelSizeFilter->SetInputForegroundValue( labelColor );
    axialLabelSizeFilter->SetFullyConnected(1);
	
    // Extracting the coronal slice containing the airways fiducial point
    OutputImageType::SizeType oneCoronalSliceSize;
    oneCoronalSliceSize[0] = airways->GetLargestPossibleRegion().GetSize(0);
    oneCoronalSliceSize[1] = 1;
    oneCoronalSliceSize[2] = 6;
  
    InputImageType::IndexType indexCoronalSlice;
    indexCoronalSlice.Fill(0);
    indexCoronalSlice[1] = index[1];
    if( indexFiducialSlice[2] >= 3 )
    {
        indexCoronalSlice[2] = indexFiducialSlice[2] - 3;
    }
    else
    {
        indexCoronalSlice[2] = indexFiducialSlice[2];
    }
    OutputImageType::RegionType coronalSlice;

    ROIFilterType::Pointer coronalAirwaysFilter = ROIFilterType::New();

    ShapeLabelType::Pointer coronalLabelSizeFilter = ShapeLabelType::New();
    coronalLabelSizeFilter->SetInputForegroundValue( labelColor );
    coronalLabelSizeFilter->SetFullyConnected(1);
  
    // Computing the sizes 
    double xSize      = 0;
    double ySize      = 0;
    bool   firstCheck = 0;
    bool   check      = 0;
    bool   decrease   = 0;
       
    double airwaysYSize  = airways->GetLargestPossibleRegion().GetSize(1) * 0.25;
    double airwaysXSize  = airways->GetLargestPossibleRegion().GetSize(0) * 2/3;

    do{
        axialSlice.SetSize( oneAxialSliceSize );
        axialSlice.SetIndex( indexAxialSlice );

        duplicatorFilter->Update();	
        airwaysAxialCopy = duplicatorFilter->GetOutput();	
            
        axialAirwaysFilter->SetInput( airwaysAxialCopy );
        axialAirwaysFilter->SetRegionOfInterest( axialSlice );
        axialAirwaysFilter->Update();

        axialLabelSizeFilter->SetInput( axialAirwaysFilter->GetOutput() );
        axialLabelSizeFilter->Update();

        if( axialLabelSizeFilter->GetOutput()->GetNumberOfLabelObjects() > 0 )
        {
            bool labelOverSize = 0; 
            unsigned int numberOfObjects = axialLabelSizeFilter->GetOutput()->GetNumberOfLabelObjects();

            for( unsigned int i = 0; i < numberOfObjects; ++i )
            {
                ySize += axialLabelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(1);
                if( ySize > airwaysYSize)
                {
                    UpperThreshold = UpperThreshold - 20;

                    thresholdConnected->SetUpper( UpperThreshold ); 
                    caster->SetInput( thresholdConnected->GetOutput() ); 
                    caster->Update();
                    airways = caster->GetOutput();
                    duplicatorFilter->Update();	
                    airwaysAxialCopy = duplicatorFilter->GetOutput();	
                    axialAirwaysFilter->SetInput( airwaysAxialCopy );
                    axialAirwaysFilter->SetRegionOfInterest( axialSlice );
                    axialAirwaysFilter->Update();
                    axialLabelSizeFilter->SetInput( axialAirwaysFilter->GetOutput() );
                    axialLabelSizeFilter->Update();
                    decrease = 1;
                    xSize = 0;
                    ySize = 0;
                    labelOverSize = 1;
                    i = numberOfObjects - 1;    
                }
            }

            if( !labelOverSize )
            {
                coronalSlice.SetIndex( indexCoronalSlice );
                coronalSlice.SetSize( oneCoronalSliceSize );

                duplicatorFilter->Update();	
                airwaysCoronalCopy = duplicatorFilter->GetOutput();	

                coronalAirwaysFilter->SetInput( airwaysCoronalCopy );
                coronalAirwaysFilter->SetRegionOfInterest( coronalSlice );
                coronalAirwaysFilter->Update();

                coronalLabelSizeFilter->SetInput( coronalAirwaysFilter->GetOutput() );
                coronalLabelSizeFilter->Update();
		 
                unsigned int numberOfObjects = coronalLabelSizeFilter->GetOutput()->GetNumberOfLabelObjects();

                for( unsigned int i = 0; i < numberOfObjects; i++ )
                {
                    xSize += coronalLabelSizeFilter->GetOutput()->GetNthLabelObject( i )->GetBoundingBox().GetSize(0);

                    if( xSize > airwaysXSize )
                    {
                        UpperThreshold = UpperThreshold - 20;

                        thresholdConnected->SetUpper( UpperThreshold ); 
                        caster->SetInput( thresholdConnected->GetOutput() );  
                        caster->Update();
                        airways = caster->GetOutput();
                        duplicatorFilter->Update();	
                        airwaysCoronalCopy = duplicatorFilter->GetOutput();	
                        coronalAirwaysFilter->SetInput( airwaysCoronalCopy );
                        coronalAirwaysFilter->SetRegionOfInterest( coronalSlice );
                        coronalAirwaysFilter->Update();
                        i = numberOfObjects - 1;	 
                        coronalLabelSizeFilter->SetInput( axialAirwaysFilter->GetOutput() );
                        coronalLabelSizeFilter->Update();
                        decrease = 1;
                        xSize = 0;
                        ySize = 0;
                    }
                }
            }
            if( xSize != 0 && ySize != 0 )
            {
                xSize = xSize + xSize * 30 / 100;
                ySize = ySize + ySize * 30 / 100;
                firstCheck = 1;
            }
        }
        else
        {
            bool isMinor = 0;
            InputImageType::SizeType    radius,regionSize;
            InputImageType::IndexType   regionIndex;
            InputImageType::RegionType  region;	  

            regionSize.Fill(3);
            regionIndex = index;
            radius.Fill(3);

            region.SetSize(regionSize);
            region.SetIndex(regionIndex);

            typedef itk::ConstNeighborhoodIterator< InputImageType > NeighborhoodIterator;
            NeighborhoodIterator iterator(radius, image, region);
	
            unsigned int counter = 0;

            while( counter < iterator.Size() && !isMinor )
            {
                if( iterator.GetPixel(counter) < UpperThreshold )
                {
                    index = iterator.GetIndex( counter );
            
                    indexCoronalSlice[1] = index[1];
                    indexCoronalSlice[2] = index[2] - 3;

                    thresholdConnected->ClearSeeds();
                    thresholdConnected->AddSeed( index );
                    thresholdConnected->Update();                
    
                    caster->SetInput( thresholdConnected->GetOutput() );
                    caster->Update();	

                    airways = caster->GetOutput();	
                    isMinor = 1;
                }
                counter++;   
            }

            if ( !isMinor && !decrease)
            {
                if( UpperThreshold < -800 )
                {
                    UpperThreshold = UpperThreshold + 50;

                    thresholdConnected->SetUpper( UpperThreshold ); 
                    thresholdConnected->Update();                
    
                    caster->SetInput( thresholdConnected->GetOutput() );
                    caster->Update();

                    airways = caster->GetOutput();
                }
                else
                {
                    std::cout<<"Please move the seed point in a different position."<<std::endl;
                    return airways;
                }
            }
            else if( !isMinor && decrease )
            {
                if( UpperThreshold < -800 )
                {
                    UpperThreshold = UpperThreshold + 1;

                    thresholdConnected->SetUpper( UpperThreshold ); 
                    thresholdConnected->Update();

                    caster->SetInput( thresholdConnected->GetOutput() ); 
                    caster->Update();

                    airways = caster->GetOutput();	
                }
                else
                {
                    std::cout<<"Please move the seed point to a different location."<<std::endl;
                    return airways;
                }
            }
			
            axialSlice.SetSize( oneAxialSliceSize );
            axialSlice.SetIndex( indexAxialSlice );

            duplicatorFilter->Update();	
            airwaysAxialCopy = duplicatorFilter->GetOutput();	
            
            axialAirwaysFilter->SetInput( airwaysAxialCopy );
            axialAirwaysFilter->SetRegionOfInterest( axialSlice );
            axialAirwaysFilter->Update();

            axialLabelSizeFilter->SetInput( axialAirwaysFilter->GetOutput() );
            axialLabelSizeFilter->Update();

            coronalSlice.SetIndex( indexCoronalSlice );
            coronalSlice.SetSize( oneCoronalSliceSize );

            duplicatorFilter->Update();	
            airwaysCoronalCopy = duplicatorFilter->GetOutput();	

            coronalAirwaysFilter->SetInput( airwaysCoronalCopy );
            coronalAirwaysFilter->SetRegionOfInterest( coronalSlice );
            coronalAirwaysFilter->Update();

            coronalLabelSizeFilter->SetInput( coronalAirwaysFilter->GetOutput() );
            coronalLabelSizeFilter->Update();
	  		
            xSize = 0;
            ySize = 0;
        }
    }
    while( !firstCheck && UpperThreshold > -1100 );
    
    duplicatorFilter->SetInputImage(airways);
    duplicatorFilter->Update();
    airwaysPrev = duplicatorFilter->GetOutput();
        
    /** INCREASING THE THRESHOLD ITERATIVELY UNTIL LEAKAGE OCCURS */

    typedef itk::SubtractImageFilter< OutputImageType,OutputImageType,OutputImageType > SubtractLabelImageType; 
    SubtractLabelImageType::Pointer addedLabel = SubtractLabelImageType::New();
  
    bool overThreshold = 0;

    ShapeLabelType::Pointer labelSizeFilter = ShapeLabelType::New();

    do{
        addedLabel->SetInput1( airways );
        addedLabel->SetInput2( airwaysPrev );
        addedLabel->Update();

        labelSizeFilter->SetInput( addedLabel->GetOutput() );
        labelSizeFilter->SetInputForegroundValue( labelColor );
        labelSizeFilter->Update();
        unsigned int numberOfObjects = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects(); 
        double       xSz             = 0;
        double       ySz             = 0;
        double       zSz             = 0;
        
        if( numberOfObjects > 0 )
        {
            for( unsigned int i = 0; i < numberOfObjects; i++ )
            {
                xSz = labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(0);
                ySz = labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(1);
                zSz = labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(2);
                if( xSz > xSize || ySz > ySize || zSz > airways->GetLargestPossibleRegion().GetSize(2) / 3 )
                {
                    if( decrease )
                    {
                        UpperThreshold = UpperThreshold - 10;
                        thresholdConnected->SetUpper( UpperThreshold ); 
                        thresholdConnected->Update();
                        caster->SetInput( thresholdConnected->GetOutput() );
                        caster->Update();
                        airways = caster->GetOutput();
                    }
                    check = 1;
                    i = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects() - 1;
                }
            }
        }
        if( !check )
        {
            if( UpperThreshold < -800 )
            {
                duplicatorFilter->SetInputImage(airways);
                duplicatorFilter->Update();
                airwaysPrev = duplicatorFilter->GetOutput();
                if( !decrease )
                {
                    UpperThreshold = UpperThreshold + 50;
                }
                else
                {
                    UpperThreshold = UpperThreshold + 10;
                }
                thresholdConnected->SetUpper( UpperThreshold ); 
                thresholdConnected->Update();
                caster->SetInput( thresholdConnected->GetOutput() );
                caster->Update();
                airways = caster->GetOutput();	
            }
            else
            {
                check = 1;
                overThreshold = 1;
            }
        }
    }
    while( !check );
  
    // Decreasing the threshold to find a better segmentation
    if( !overThreshold && !decrease )
    {
        while( check && UpperThreshold > -1100 )
        {
            UpperThreshold = UpperThreshold - 10;

            thresholdConnected->SetUpper( UpperThreshold ); 
            caster->SetInput( thresholdConnected->GetOutput() );  
            airways = caster->GetOutput();			
            caster->Update();

            addedLabel->SetInput1( airways );
            addedLabel->SetInput2( airwaysPrev );
            addedLabel->Update();

            labelSizeFilter->SetInput( addedLabel->GetOutput() );
            labelSizeFilter->Update();

            unsigned int count = 0;
            unsigned int numberOfObjects = labelSizeFilter->GetOutput()->GetNumberOfLabelObjects();

            if( numberOfObjects > 0 )
            {
                for( unsigned int i = 0; i < numberOfObjects; i++ )
                {
                    if( labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(0) < xSize && 
                        labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(1) < ySize &&
                        labelSizeFilter->GetOutput()->GetNthLabelObject(i)->GetBoundingBox().GetSize(2) < airways->GetLargestPossibleRegion().GetSize(2) / 3)
                    {
                        count++;
                    }
                }
                if( count == numberOfObjects )
                {
                    check = 0;
                }
            }
            else
            {
                check = 0;
            }
        }
    }

    std::cout<<"Threshold: "<<UpperThreshold<<std::endl;
    return airways;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef itk::Image< LabelType,       Dim > 				LabelImageType;
  typedef itk::Image< MSOutputPixelType, Dim > 				MSOutputImageType;

  typedef itk::ImageFileReader<InputImageType>  			ReaderType;

  typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > 	I2LType;
  typedef I2LType::OutputImageType 					LabelMapType;
  typedef LabelMapType::LabelObjectType					ShapeLabelObjectType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( InputVolume.c_str() );
  reader->Update();
  InputImageType::Pointer inputImage = reader->GetOutput();  

  typedef itk::ImageDuplicator< InputImageType >InputDuplicatorType;
  InputDuplicatorType::Pointer inputDuplicator = InputDuplicatorType::New();
  inputDuplicator->SetInputImage(inputImage);
  inputDuplicator->Update();
  InputImageType::Pointer clonedInput = inputDuplicator->GetOutput();

  // Median filter
  typedef itk::MedianImageFilter<InputImageType, InputImageType > MedianFilterType;
  MedianFilterType::InputSizeType medianRadius;
  medianRadius.Fill(1);
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  medianFilter->SetRadius( medianRadius );
  medianFilter->SetInput( inputImage );

  // Binarize the image based on the specified threshold
  typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > BinaryThresholdFilterType;
  BinaryThresholdFilterType::Pointer binaryThresholdFilter = BinaryThresholdFilterType::New();
  binaryThresholdFilter->SetInput( 0, medianFilter->GetOutput() );
  binaryThresholdFilter->SetInsideValue( lungsColor );
  binaryThresholdFilter->SetOutsideValue( 0 );
  binaryThresholdFilter->SetLowerThreshold(Lower);
  binaryThresholdFilter->SetUpperThreshold(Upper);
  
  // Remove holes not connected to the boundary of the image.
  typedef itk::BinaryFillholeImageFilter< OutputImageType > FillHolesFilterType;
  FillHolesFilterType::Pointer holeFillingFilter = FillHolesFilterType::New();
  holeFillingFilter->SetInput( binaryThresholdFilter->GetOutput() );
  holeFillingFilter->SetForegroundValue( lungsColor );
  holeFillingFilter->SetFullyConnected( 1 );

  // Convert the binary image into a label map to valuate the shape attributes 
  I2LType::Pointer i2l = I2LType::New();
  i2l->SetInput( holeFillingFilter->GetOutput() );
  i2l->SetInputForegroundValue( lungsColor );
  i2l->SetFullyConnected( 1 );
  i2l->Update();

  LabelMapType::Pointer labelMap = i2l->GetOutput();
  //std::cout << "Number of objects beforehand: " << labelMap->GetNumberOfLabelObjects() << std::endl;
  std::vector<unsigned long> labelsToRemove;

  for(unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
  {
     ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(n);
     if(labelObject->GetNumberOfPixels() <= 2000 || labelObject->GetPerimeter() <= 200 )
     {
        labelsToRemove.push_back(labelObject->GetLabel());
     }
     else 
     {
        if(labelObject->GetBoundingBox().GetIndex()[0] == 0 || labelObject->GetBoundingBox().GetIndex()[1] == 0 || labelObject->GetBoundingBox().GetIndex()[2] == 0)
        {
          labelsToRemove.push_back(labelObject->GetLabel());
        }
     }
  }
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i)
  {
      labelMap->RemoveLabel(labelsToRemove[i]);
  }
  
  //std::cout << "Number of objects afterwards: " << labelMap->GetNumberOfLabelObjects() << std::endl;

  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, OutputImageType > LabelMapToBinaryImageType; 
  LabelMapToBinaryImageType::Pointer labelToBinaryFilter = LabelMapToBinaryImageType::New();
  labelToBinaryFilter->SetInput( labelMap );
  labelToBinaryFilter->SetBackgroundValue( 0 );
  labelToBinaryFilter->SetForegroundValue( lungsColor );  
  labelToBinaryFilter->Update();

  typedef itk::ImageDuplicator<OutputImageType> LungDuplicatorFilterType;

  LungDuplicatorFilterType::Pointer lungDuplicatorFilter = LungDuplicatorFilterType::New();
  lungDuplicatorFilter->SetInputImage(labelToBinaryFilter->GetOutput());
  lungDuplicatorFilter->Update();
  
  unsigned int Zposition = 0;
  if(clonedInput->GetLargestPossibleRegion().GetSize(2) > 150)
  {
    Zposition = int(clonedInput->GetLargestPossibleRegion().GetSize(2) - (clonedInput->GetLargestPossibleRegion().GetSize(2)*0.2));
  }
  else
  {
    Zposition = int(clonedInput->GetLargestPossibleRegion().GetSize(2) - (clonedInput->GetLargestPossibleRegion().GetSize(2)*0.1));
  }
  //std::cout<<"Zposition: "<<Zposition<<std::endl;

  // Extract a slice to find the trachea location
  OutputImageType::SizeType sliceSize = clonedInput->GetLargestPossibleRegion().GetSize();
  sliceSize[2] = 1;
  std::cout<<"sliceSize: "<<sliceSize<<std::endl;

  OutputImageType::IndexType sliceIndex;
  sliceIndex.Fill(0);
  sliceIndex[2] = Zposition;

  OutputImageType::RegionType sliceRegion( sliceIndex, sliceSize );

  typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > ROIFilterType;

  ROIFilterType::Pointer sliceExtractor = ROIFilterType::New();
  sliceExtractor->SetRegionOfInterest( sliceRegion );
  sliceExtractor->SetInput( labelToBinaryFilter->GetOutput() );

  I2LType::Pointer slice2labels = I2LType::New();
  slice2labels->SetInput( sliceExtractor->GetOutput() );
  slice2labels->SetInputForegroundValue( lungsColor );
  slice2labels->SetFullyConnected( 1 );
  slice2labels->SetComputeFeretDiameter( 1 );
  slice2labels->Update();

  std::cout << "Number of labels on slice: " << slice2labels->GetOutput()->GetNumberOfLabelObjects() << std::endl;

  float minDiameter = itk::NumericTraits< float >::max();
  unsigned long minLabel = 0;
  float minRoundness = 1.0;
  for(unsigned int i = 0; i < slice2labels->GetOutput()->GetNumberOfLabelObjects(); ++i)
  {
    ShapeLabelObjectType::Pointer labelObject = slice2labels->GetOutput()->GetNthLabelObject(i);
    float currentDiam = labelObject->GetFeretDiameter();
    float currentRoundness = labelObject->GetRoundness();
    //std::cout<<"currentDiam: "<<currentDiam<<std::endl;
    if( currentDiam > 10 && currentDiam < 40 && currentDiam < minDiameter && currentRoundness < minRoundness )
    {
      minDiameter = currentDiam;
      minLabel = i;
      minRoundness = currentRoundness;
    }   
  }

  itk::Point<double, Dim> centroid = slice2labels->GetOutput()->GetNthLabelObject(minLabel)->GetCentroid();
  std::cout<<"centroid: "<<centroid<<std::endl;

  clonedInput = inputDuplicator->GetOutput();
  
  /** AIRWAYS SEGMENTATION */
  InputImageType::IndexType airwayIndex;
  clonedInput->TransformPhysicalPointToIndex(centroid, airwayIndex);

  std::cout<<"index: "<<airwayIndex<<std::endl;

  InputImageType::SizeType  cropSize;	                    
  InputImageType::IndexType cropIndex;              

  cropSize[0] = 100; 
  cropIndex[0] = airwayIndex[0] - 50;   
  if( clonedInput->GetLargestPossibleRegion().GetSize(2) > 200 )
  {
      cropSize[2] = clonedInput->GetLargestPossibleRegion().GetSize(2) - 100;
      cropIndex[2] = 100;
  }
  else
  {
      cropSize[2] = clonedInput->GetLargestPossibleRegion().GetSize(2) - (clonedInput->GetLargestPossibleRegion().GetSize(2)/4); 
      cropIndex[2] = clonedInput->GetLargestPossibleRegion().GetSize(2)/4;
  }

  cropSize[1] = clonedInput->GetLargestPossibleRegion().GetSize(1);        
  cropIndex[1] = 0;		     
    
  InputImageType::RegionType cropRegion;                                                  
  cropRegion.SetSize(  cropSize  );                                                                
  cropRegion.SetIndex( cropIndex );
  
  std::cout<<"cropIndex: "<<cropIndex<<std::endl;
  std::cout<<"cropSize: "<<cropSize<<std::endl;
  // Cropping the trachea 
  typedef itk::RegionOfInterestImageFilter< InputImageType, InputImageType > inputROIFilterType;  
  inputROIFilterType::Pointer airwayROIFilter = inputROIFilterType::New();	                 
  
  airwayROIFilter->SetInput( clonedInput );						          
  airwayROIFilter->SetRegionOfInterest( cropRegion );					 
  airwayROIFilter->Update();

  InputImageType::IndexType fiducialSlice;
  fiducialSlice.Fill(0);
  fiducialSlice[2] = cropSize[2] - ( clonedInput->GetLargestPossibleRegion().GetSize(2) - airwayIndex[2] ); 
  
  double fidPs  = double( fiducialSlice[2] );
  double trSz   = double( cropSize[2] );
  double ratio  = fidPs/trSz;

  if(ratio >= 0.85)
  {
      fiducialSlice[2] = cropSize[2]*0.8;
  }

  InputImageType::IndexType croppedAirwayIndex;
  airwayROIFilter->GetOutput()->TransformPhysicalPointToIndex(centroid, croppedAirwayIndex);

  std::cout<<"index: "<<airwayIndex<<std::endl;

  OutputImageType::Pointer airwayImage = AirwaySegmentation(airwayROIFilter->GetOutput(), fiducialSlice, croppedAirwayIndex, lungsColor);

  // Remove small (i.e., smaller than the structuring element) holes and tube like structures in the interior or at the boundaries of the image.
  typedef itk::BinaryBallStructuringElement< OutputPixelType, Dim > StructuringElementType;  	
  StructuringElementType structElement;
  StructuringElementType::SizeType radius;
  radius.Fill( 10 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
	
  typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > ClosingFilterType;
  ClosingFilterType::Pointer airwayClosing = ClosingFilterType::New();

  airwayClosing->SetInput( airwayImage );
  airwayClosing->SetKernel( structElement );
  airwayClosing->SetForegroundValue( lungsColor );
  airwayClosing->SetSafeBorder( 1 );
  airwayClosing->Update();
  airwayImage = airwayClosing->GetOutput();

  /* FINDING THE CARINA 
  clonedInput = inputDuplicator->GetOutput();
  int          value;
  int          xDist;
  int          yDist; 

  unsigned int xCarinaPosition = 0;
  unsigned int yCarinaPosition = 0;
  unsigned int carinaIdx       = 0;

  InputImageType::IndexType carinaIdxPrevSlice;
  carinaIdxPrevSlice.Fill(0);
  InputImageType::IndexType tmpIdx;
  tmpIdx.Fill(0);

  typedef itk::ImageSliceIteratorWithIndex<OutputImageType> SliceIterator;
  SliceIterator sIt(airwayImage, airwayImage->GetLargestPossibleRegion());

  sIt.SetFirstDirection(0);
  sIt.SetSecondDirection(1);

  sIt.GoToBegin();
        
  unsigned int limit = fiducialSlice[2] - 5;
 
  while( sIt.GetIndex()[2] < limit )
  {
    unsigned int yMaxPos     = 0;
    unsigned int yMinPos     = airwayImage->GetLargestPossibleRegion().GetSize(1); 

    unsigned int yMidPos     = 0; 

    unsigned int yCurrentMax = 0;
    unsigned int yCurrentMin = 0;

    unsigned int prevPos     = 0;
    unsigned int prevLine    = 0;
    
    unsigned int prevDiff    = 0;

    while( !sIt.IsAtEndOfSlice() )
    {
        while( !sIt.IsAtEndOfLine() )
        {
            value = sIt.Get();
            if( value == lungsColor )
            {
	      double d = sIt.GetIndex()[1] - prevLine;
              yDist = abs(d);
              if( prevLine != 0 && yDist <= 5 )
                {
                    if( sIt.GetIndex()[1] > yMaxPos )
                    {
                        yMaxPos = sIt.GetIndex()[1];
                    }
                    if( sIt.GetIndex()[1] < yMinPos )
                    {
                        yMinPos = sIt.GetIndex()[1];
                    }
                    prevLine = sIt.GetIndex()[1];

                    unsigned int diff = yMaxPos - yMinPos;
                    if( diff > prevDiff)
                    {
                        yCurrentMax = yMaxPos;
                        yCurrentMin = yMinPos;
                     }
                }
                else if( prevLine != 0 && yDist >= 5 )
                {
                    prevDiff = yCurrentMax - yCurrentMin;

                    yMinPos = airwayImage->GetLargestPossibleRegion().GetSize(1); 
                    yMaxPos = 0; 
                    if( sIt.GetIndex()[1] > yMaxPos )
                    {
                        yMaxPos = sIt.GetIndex()[1];
                    }
                    if( sIt.GetIndex()[1] < yMinPos )
                    {
                        yMinPos = sIt.GetIndex()[1];
                    }
                }
                prevLine = sIt.GetIndex()[1];
            }
            ++sIt;
        }
        sIt.NextLine();
    }

    if( yCurrentMax > yCurrentMin )
    {
        yMidPos = yCurrentMin + (yCurrentMax - yCurrentMin) / 2;
    }

    sIt.GoToBeginOfSlice();
    while( !sIt.IsAtEndOfSlice() )
    {
        while( !sIt.IsAtEndOfLine() )
        {
            value = sIt.Get();
            if( value == lungsColor )
            {
                xDist = sIt.GetIndex()[0] - prevPos;
                if( prevPos != 0 && sIt.GetIndex()[1] == yMidPos 
                    && xDist >= 10 && xDist < 20
                    && sIt.GetIndex()[0] > int(cropSize[0]/3) )
                {
                    carinaIdx       = sIt.GetIndex()[2];
                    xCarinaPosition = prevPos + (xDist/2);
                    yCarinaPosition = yMidPos;
                }
                prevPos = sIt.GetIndex()[0];
            }
            ++sIt;
        }
        sIt.NextLine();
    }        
      carinaIdxPrevSlice[0] = xCarinaPosition; 
      carinaIdxPrevSlice[1] = yCarinaPosition;
      carinaIdxPrevSlice[2] = carinaIdx;

      sIt.NextSlice();
  }

  carinaIdx -= 10;
  cropSize[2] -= carinaIdx;
  cropIndex[2] += carinaIdx;
  fiducialSlice[2] = cropSize[2] - ( clonedInput->GetLargestPossibleRegion().GetSize(2) - airwayIndex[2] ); 

  fidPs  = double( fiducialSlice[2] );
  trSz   = double( cropSize[2] );
  ratio  = fidPs/trSz;

  if(ratio >= 0.85)
  {
      fiducialSlice[2] = cropSize[2]*0.8;
  }

  bool exit = 0;
  sIt.GoToBegin();
        
  while( !exit && !sIt.IsAtEnd() )
  {
      while( !exit && !sIt.IsAtEndOfSlice() )
      {
          while( !exit && !sIt.IsAtEndOfLine() )
          {
              value = sIt.Get();
              if( value == lungsColor && sIt.GetIndex()[2] >= (carinaIdx + 10 ) && sIt.GetIndex()[2] <= (cropSize[2] - 5) )
              {
                  if( sIt.GetIndex()[0] == 0 || sIt.GetIndex()[0] == cropSize[0] )
                  {
                      cropSize[0] += 10; 
                      cropIndex[0] -= 5;
                      xCarinaPosition += 5;
                      exit = 1;
                  }
              }
              ++sIt;
          }
          sIt.NextLine();
      }
      sIt.NextSlice();
  }

  cropRegion.SetSize(  cropSize  );                                                                
  cropRegion.SetIndex( cropIndex );                                                          

  airwayROIFilter->SetInput( clonedInput );						          
  airwayROIFilter->SetRegionOfInterest( cropRegion );					 
  airwayROIFilter->Update();

  airwayROIFilter->GetOutput()->TransformPhysicalPointToIndex(centroid, croppedAirwayIndex);
  airwayImage = AirwaySegmentation(airwayROIFilter->GetOutput(), fiducialSlice, croppedAirwayIndex, lungsColor);*/

  typedef itk::BinaryDilateImageFilter <OutputImageType, OutputImageType, StructuringElementType> BinaryDilateImageFilterType;
  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(airwayImage);
  radius.Fill( 5 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
  dilateFilter->SetForegroundValue( lungsColor );
  dilateFilter->SetKernel(structElement);
  dilateFilter->Update();

  OutputImageType::Pointer pasteImage = OutputImageType::New(); 

  pasteImage->SetRegions( inputImage->GetRequestedRegion() );                                   
  pasteImage->SetBufferedRegion( inputImage->GetBufferedRegion() );
  pasteImage->SetLargestPossibleRegion( inputImage->GetLargestPossibleRegion() );
  pasteImage->CopyInformation( inputImage );
  pasteImage->Allocate();
  pasteImage->FillBuffer(itk::NumericTraits< OutputPixelType >::Zero); 

  typedef itk::PasteImageFilter< OutputImageType, OutputImageType > PasteImageFilterType;
  PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New ();  
      
  pasteFilter->SetSourceImage(dilateFilter->GetOutput());                                                 
  pasteFilter->SetDestinationImage(pasteImage);                                               
  pasteFilter->SetSourceRegion(dilateFilter->GetOutput()->GetLargestPossibleRegion());                    
  pasteFilter->SetDestinationIndex(cropIndex);

  OutputImageType::Pointer lungImage = lungDuplicatorFilter->GetOutput();

  typedef itk::NotImageFilter<OutputImageType,OutputImageType> NotImageFilterType;
  NotImageFilterType::Pointer notFilter = NotImageFilterType::New();
  notFilter->SetInput(pasteFilter->GetOutput());

  typedef itk::AndImageFilter<OutputImageType> AndImageFilterType;
  AndImageFilterType::Pointer andFilter  = AndImageFilterType::New();
  andFilter->SetInput(0, lungImage);
  andFilter->SetInput(1, notFilter->GetOutput());

  typedef itk::ConnectedComponentImageFilter <OutputImageType, OutputImageType>  ConnectedComponentImageFilterType;

  ConnectedComponentImageFilterType::Pointer connectedComponentFilter = ConnectedComponentImageFilterType::New ();
  connectedComponentFilter->SetInput(andFilter->GetOutput());

  typedef itk::RelabelComponentImageFilter<OutputImageType, OutputImageType> RelabelImageFilterType;
  RelabelImageFilterType::Pointer relabelFilter = RelabelImageFilterType::New();

  relabelFilter->SetInput(connectedComponentFilter->GetOutput());
  relabelFilter->SetMinimumObjectSize(5000);  

  typedef itk::ThresholdImageFilter<OutputImageType> ThresholdImageFilterType;
  ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();

  thresholdFilter->SetInput(relabelFilter->GetOutput());
  thresholdFilter->ThresholdOutside(0,2);
  thresholdFilter->SetOutsideValue(0);

  ConnectedComponentImageFilterType::Pointer lungsCCFilter = ConnectedComponentImageFilterType::New ();
  lungsCCFilter->SetInput(thresholdFilter->GetOutput());
  lungsCCFilter->Update();
  OutputImageType::Pointer connCompImage = lungsCCFilter->GetOutput();

  //SbSCCFilter->SetInput(connCompImage);
  typedef itk::SliceBySliceImageFilter<OutputImageType, OutputImageType> SliceBySliceFilterType;
  typedef itk::ConnectedComponentImageFilter <SliceBySliceFilterType::InternalInputImageType, SliceBySliceFilterType::InternalOutputImageType>  					SbSCCImageFilterType;
  SbSCCImageFilterType::Pointer SbSCCFilter = SbSCCImageFilterType::New();

  SliceBySliceFilterType::Pointer sliceBySliceFilter = SliceBySliceFilterType::New();
  sliceBySliceFilter->SetInput(connCompImage);
  sliceBySliceFilter->SetFilter(SbSCCFilter);
  sliceBySliceFilter->SetDimension(2);
  sliceBySliceFilter->Update();
  OutputImageType::Pointer SbSImage = sliceBySliceFilter->GetOutput();

  typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > SbSROIFilterType;  
  SbSROIFilterType::Pointer SbSExtractor = SbSROIFilterType::New();
  SbSExtractor->SetInput(SbSImage);
  InputImageType::SizeType SbSSize = SbSImage->GetLargestPossibleRegion().GetSize();
  SbSSize[2] = 1;
  InputImageType::IndexType SbSIndex;
  SbSIndex.Fill(0);
  InputImageType::RegionType SbSRegion;
  SbSRegion.SetSize(  SbSSize  );                                                                

  typedef itk::MinimumMaximumImageCalculator<OutputImageType> MinMaxCalculatorType;
  MinMaxCalculatorType::Pointer minMaxCalculator;

  std::cout<<SbSImage->GetLargestPossibleRegion().GetSize(2)<<std::endl;
  for(unsigned int i = 0; i < SbSImage->GetLargestPossibleRegion().GetSize(2); ++i)
  {
    /*SbSIndex[2] = i;
    SbSRegion.SetIndex( SbSIndex );
    SbSExtractor->SetRegionOfInterest( SbSRegion );		 
    SbSExtractor->Update();
    minMaxCalculator->SetImage(SbSExtractor->GetOutput());
    minMaxCalculator->ComputeMaximum();
    if(minMaxCalculator->GetMaximum()==1)
    {*/
      std::cout<<"slice: "<<i<<std::endl;
    //}

  }
  /*typedef itk::ImageMomentsCalculator<OutputImageType> ImageMomentsCalculatorType;
  ImageMomentsCalculatorType::Pointer imageMomentsCalculator = ImageMomentsCalculatorType::New();
  imageMomentsCalculator->SetImage(connCompImage);
  imageMomentsCalculator->Compute();
  ImageMomentsCalculatorType::VectorType CoG = imageMomentsCalculator->GetCenterOfGravity();

  std::cout<<"Axial slices range: "<<CoG[0]-10<<"->"<<CoG[0]+10<<std::endl;*/


  /*radius.Fill( 15 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
  ClosingFilterType::Pointer finalClosing = ClosingFilterType::New();
  finalClosing->SetKernel( structElement );
  //finalClosing->SetSafeBorder( 1 );

  for(unsigned int i = 1; i <= 2; ++i)
  {
    finalClosing->SetInput( connCompImage );
    finalClosing->SetForegroundValue( i );
    finalClosing->Update();

    connCompImage = finalClosing->GetOutput();
  }  

  typedef itk::SubtractImageFilter< OutputImageType,OutputImageType,OutputImageType > SubtractLabelImageType; 
  SubtractLabelImageType::Pointer subtractFilter = SubtractLabelImageType::New();
  subtractFilter->SetInput1(connCompImage);
  subtractFilter->SetInput2(lungsCCFilter->GetOutput());
  subtractFilter->Update();
  OutputImageType::Pointer diffImage = subtractFilter->GetOutput();

  typedef itk::ImageRegionConstIteratorWithIndex<OutputImageType> labelDifferenceIterator;
  labelDifferenceIterator ldIt(diffImage, diffImage->GetLargestPossibleRegion());
  OutputImageType::Pointer ccImage = lungsCCFilter->GetOutput();

  for( ldIt.GoToBegin(); !ldIt.IsAtEnd(); ++ldIt)
  {
    if(ldIt.Get() != 0 && ldIt.Get() != 1 && ldIt.Get() != 2)
    {
      InputImageType::IndexType idx = ldIt.GetIndex();
      connCompImage->SetPixel(idx, ccImage->GetPixel(idx));
    }
  }*/

  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( OutputVolume.c_str() );
  writer->SetInput( sliceBySliceFilter->GetOutput() );
  writer->SetUseCompression(1);
  try
  {   
        writer->Update();
  }
  catch ( itk::ExceptionObject & e )
  {
    	std::cerr << "exception in file writer " << std::endl;
	std::cerr << e.GetDescription() << std::endl;
	std::cerr << e.GetLocation() << std::endl;
	return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
