/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $HeadURL: $
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
#include "itkBinaryErodeImageFilter.h"

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

#include "itkSliceBySliceImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "itkLabelContourImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkAddImageFilter.h"

#include "itkTimeProbe.h"

#include "cipChestConventions.h"
#include "cipChestConventions.cxx"

#include "LungSegmentationCLICLP.h"

const unsigned int Dim = 3;

typedef signed short 	InputPixelType;
typedef unsigned int 	LabelPixelType;
typedef float			MSOutputPixelType;

typedef itk::Image< InputPixelType,  Dim > InputImageType;
typedef itk::Image< LabelPixelType, Dim >  LabelImageType;

LabelImageType::Pointer AirwaySegmentation( InputImageType::Pointer image, InputImageType::IndexType indexFiducialSlice, 
                                             InputImageType::IndexType index, int labelColor = 1 )
{
    LabelImageType::Pointer airways 		= LabelImageType::New(); 
    LabelImageType::Pointer airwaysPrev 	= LabelImageType::New();
	 		
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

    typedef itk::CastImageFilter<InputImageType, LabelImageType> CastingFilterType;  
    CastingFilterType::Pointer  caster = CastingFilterType::New();		
	
    caster->SetInput( thresholdConnected->GetOutput() );  
    caster->Update();
    airways = caster->GetOutput();
    
    /** COMPUTING THE LABEL SIZES */ 	                  
    LabelImageType::Pointer airwaysAxialCopy = LabelImageType::New(); 
    LabelImageType::Pointer airwaysCoronalCopy = LabelImageType::New(); 

    typedef itk::ImageDuplicator<LabelImageType> DuplicatorFilterType;

    DuplicatorFilterType::Pointer duplicatorFilter = DuplicatorFilterType::New();
    duplicatorFilter->SetInputImage(airways);
    duplicatorFilter->Update();

    // Extracting the axial slice containing the airways fiducial point
    LabelImageType::SizeType  oneAxialSliceSize;
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
	
    LabelImageType::RegionType axialSlice;
 
    typedef itk::RegionOfInterestImageFilter< LabelImageType, LabelImageType > ROIFilterType;
    ROIFilterType::Pointer axialAirwaysFilter = ROIFilterType::New();

    typedef itk::BinaryImageToShapeLabelMapFilter< LabelImageType > ShapeLabelType;	
    ShapeLabelType::Pointer axialLabelSizeFilter = ShapeLabelType::New();

    axialLabelSizeFilter->SetInputForegroundValue( labelColor );
    axialLabelSizeFilter->SetFullyConnected(1);
	
    // Extracting the coronal slice containing the airways fiducial point
    LabelImageType::SizeType oneCoronalSliceSize;
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
    LabelImageType::RegionType coronalSlice;

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

    typedef itk::SubtractImageFilter< LabelImageType,LabelImageType,LabelImageType > SubtractLabelImageType; 
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
  itk::TimeProbe clockWriting;
  clockWriting.Start();

  typedef itk::Image< MSOutputPixelType, Dim > 						MSOutputImageType;

  typedef itk::ImageFileReader<InputImageType>						ReaderType;

  typedef itk::BinaryImageToShapeLabelMapFilter< LabelImageType >	I2LType;
  typedef I2LType::OutputImageType 									LabelMapType;
  typedef LabelMapType::LabelObjectType								ShapeLabelObjectType;

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
  typedef itk::BinaryThresholdImageFilter< InputImageType, LabelImageType > BinaryThresholdFilterType;
  BinaryThresholdFilterType::Pointer binaryThresholdFilter = BinaryThresholdFilterType::New();
  binaryThresholdFilter->SetInput( 0, medianFilter->GetOutput() );
  binaryThresholdFilter->SetInsideValue( lungsColor );
  binaryThresholdFilter->SetOutsideValue( 0 );
  binaryThresholdFilter->SetLowerThreshold(Lower);
  binaryThresholdFilter->SetUpperThreshold(Upper);
  
  // Remove holes not connected to the boundary of the image.
  typedef itk::BinaryFillholeImageFilter< LabelImageType > FillHolesFilterType;
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

  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, LabelImageType > LabelMapToBinaryImageType; 
  LabelMapToBinaryImageType::Pointer labelToBinaryFilter = LabelMapToBinaryImageType::New();
  labelToBinaryFilter->SetInput( labelMap );
  labelToBinaryFilter->SetBackgroundValue( 0 );
  labelToBinaryFilter->SetForegroundValue( lungsColor );  
  labelToBinaryFilter->Update();

  typedef itk::ImageDuplicator<LabelImageType> LungDuplicatorFilterType;

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
  LabelImageType::SizeType sliceSize = clonedInput->GetLargestPossibleRegion().GetSize();
  sliceSize[2] = 1;

  LabelImageType::IndexType sliceIndex;
  sliceIndex.Fill(0);
  sliceIndex[2] = Zposition;

  LabelImageType::RegionType sliceRegion( sliceIndex, sliceSize );

  typedef itk::RegionOfInterestImageFilter< LabelImageType, LabelImageType > ROIFilterType;

  ROIFilterType::Pointer sliceExtractor = ROIFilterType::New();
  sliceExtractor->SetRegionOfInterest( sliceRegion );
  sliceExtractor->SetInput( labelToBinaryFilter->GetOutput() );

  I2LType::Pointer slice2labels = I2LType::New();
  slice2labels->SetInput( sliceExtractor->GetOutput() );
  slice2labels->SetInputForegroundValue( lungsColor );
  slice2labels->SetFullyConnected( 1 );
  slice2labels->SetComputeFeretDiameter( 1 );
  slice2labels->Update();

  //std::cout << "Number of labels on slice: " << slice2labels->GetOutput()->GetNumberOfLabelObjects() << std::endl;

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
  //std::cout<<"centroid: "<<centroid<<std::endl;

  clonedInput = inputDuplicator->GetOutput();
  
  /** Use airways segmentation method to remove main bronchi from lungs segmentation */
  InputImageType::IndexType airwayIndex;
  clonedInput->TransformPhysicalPointToIndex(centroid, airwayIndex);

  //std::cout<<"index: "<<airwayIndex<<std::endl;

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

  //std::cout<<"index: "<<airwayIndex<<std::endl;

  LabelImageType::Pointer airwayImage = AirwaySegmentation(airwayROIFilter->GetOutput(), fiducialSlice, croppedAirwayIndex, lungsColor);

  // Remove small (i.e., smaller than the structuring element) holes and tube like structures in the interior or at the boundaries of the image.
  typedef itk::BinaryBallStructuringElement< LabelPixelType, Dim > StructuringElementType;  	
  StructuringElementType structElement;
  StructuringElementType::SizeType radius;
  radius.Fill( 10 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
	
  typedef itk::BinaryMorphologicalClosingImageFilter < LabelImageType, LabelImageType, StructuringElementType > ClosingFilterType;
  ClosingFilterType::Pointer airwayClosing = ClosingFilterType::New();

  airwayClosing->SetInput( airwayImage );
  airwayClosing->SetKernel( structElement );
  airwayClosing->SetForegroundValue( lungsColor );
  airwayClosing->SetSafeBorder( 1 );
  airwayClosing->Update();
  airwayImage = airwayClosing->GetOutput();
 
  typedef itk::BinaryDilateImageFilter <LabelImageType, LabelImageType, StructuringElementType> BinaryDilateImageFilterType;
  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(airwayImage);
  radius.Fill( 5 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
  dilateFilter->SetForegroundValue( lungsColor );
  dilateFilter->SetKernel(structElement);
  dilateFilter->Update();

  LabelImageType::Pointer pasteImage = LabelImageType::New(); 

  pasteImage->SetRegions( inputImage->GetRequestedRegion() );                                   
  pasteImage->SetBufferedRegion( inputImage->GetBufferedRegion() );
  pasteImage->SetLargestPossibleRegion( inputImage->GetLargestPossibleRegion() );
  pasteImage->CopyInformation( inputImage );
  pasteImage->Allocate();
  pasteImage->FillBuffer(itk::NumericTraits< LabelPixelType >::Zero); 

  typedef itk::PasteImageFilter< LabelImageType, LabelImageType > PasteImageFilterType;
  PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New();  
      
  pasteFilter->SetSourceImage(dilateFilter->GetOutput());                                                 
  pasteFilter->SetDestinationImage(pasteImage);                                               
  pasteFilter->SetSourceRegion(dilateFilter->GetOutput()->GetLargestPossibleRegion());                    
  pasteFilter->SetDestinationIndex(cropIndex);

  LabelImageType::Pointer lungImage = lungDuplicatorFilter->GetOutput();

  typedef itk::NotImageFilter<LabelImageType,LabelImageType> NotImageFilterType;
  NotImageFilterType::Pointer notFilter = NotImageFilterType::New();
  notFilter->SetInput(pasteFilter->GetOutput());

  typedef itk::AndImageFilter<LabelImageType> AndImageFilterType;
  AndImageFilterType::Pointer andFilter  = AndImageFilterType::New();
  andFilter->SetInput(0, lungImage);
  andFilter->SetInput(1, notFilter->GetOutput());

  typedef itk::ConnectedComponentImageFilter <LabelImageType, LabelImageType> ConnectedComponentImageFilterType;

  ConnectedComponentImageFilterType::Pointer connectedComponentFilter = ConnectedComponentImageFilterType::New ();
  connectedComponentFilter->SetInput(andFilter->GetOutput());

  typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType> RelabelImageFilterType;
  RelabelImageFilterType::Pointer relabelFilter = RelabelImageFilterType::New();

  relabelFilter->SetInput(connectedComponentFilter->GetOutput());
  relabelFilter->SetMinimumObjectSize(5000);  

  typedef itk::ThresholdImageFilter<LabelImageType> ThresholdImageFilterType;
  ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();

  thresholdFilter->SetInput(relabelFilter->GetOutput());
  thresholdFilter->ThresholdOutside(0,2);
  thresholdFilter->SetOutsideValue(0);

  ConnectedComponentImageFilterType::Pointer lungsCCFilter = ConnectedComponentImageFilterType::New();
  lungsCCFilter->SetInput(thresholdFilter->GetOutput());
  lungsCCFilter->Update();
  LabelImageType::Pointer connCompImage = lungsCCFilter->GetOutput();

  clockWriting.Stop();
  std::cout << "Lung Segm.: " << clockWriting.GetTotal() << std::endl;

  itk::TimeProbe clockSeparation;
  clockSeparation.Start();

  typedef itk::Image< LabelPixelType, Dim-1 > SliceOutputImageType;
  typedef itk::ConnectedComponentImageFilter <SliceOutputImageType, SliceOutputImageType>  SbSCCImageFilterType;
  SbSCCImageFilterType::Pointer SbSCCFilter = SbSCCImageFilterType::New();

  typedef itk::RelabelComponentImageFilter <SliceOutputImageType, SliceOutputImageType> SbSRelabelImageFilterType;
  SbSRelabelImageFilterType::Pointer SbSRelabelFilter = SbSRelabelImageFilterType::New();

  SbSRelabelFilter->SetInput(SbSCCFilter->GetOutput());
  SbSRelabelFilter->SetMinimumObjectSize(2000); 

  typedef itk::SliceBySliceImageFilter<LabelImageType, LabelImageType, SbSCCImageFilterType, SbSRelabelImageFilterType > CCSliceBySliceFilterType;

  CCSliceBySliceFilterType::Pointer ccSliceBySliceFilter = CCSliceBySliceFilterType::New();
  ccSliceBySliceFilter->SetInput(connCompImage);
  ccSliceBySliceFilter->SetInputFilter(SbSCCFilter);
  ccSliceBySliceFilter->SetOutputFilter(SbSRelabelFilter);
  ccSliceBySliceFilter->SetDimension(2);
  ccSliceBySliceFilter->Update();
  LabelImageType::Pointer ccSbSImage = ccSliceBySliceFilter->GetOutput();

  typedef itk::ImageDuplicator<LabelImageType> DuplicatorFilterType;
  DuplicatorFilterType::Pointer ccSbSDuplicatorFilter = DuplicatorFilterType::New();
  ccSbSDuplicatorFilter->SetInputImage(ccSbSImage);
  ccSbSDuplicatorFilter->Update();
  LabelImageType::Pointer ccSbSImageCopy = ccSbSDuplicatorFilter->GetOutput();

  typedef itk::RegionOfInterestImageFilter< LabelImageType, LabelImageType > SbSROIFilterType;  
  SbSROIFilterType::Pointer SbSExtractor = SbSROIFilterType::New();
  SbSExtractor->SetInput(ccSbSImage);
  InputImageType::SizeType SbSSize = ccSbSImage->GetLargestPossibleRegion().GetSize();
  SbSSize[2] = 1;
  InputImageType::IndexType SbSIndex;
  SbSIndex.Fill(0);
  InputImageType::RegionType SbSRegion;
  SbSRegion.SetSize( SbSSize );                                                                

  typedef itk::MinimumMaximumImageCalculator<LabelImageType> MinMaxCalculatorType;

  unsigned int imageSize = ccSbSImage->GetLargestPossibleRegion().GetSize(2);

  slice2labels->SetInputForegroundValue( 1 );
  slice2labels->SetFullyConnected( 1 );
  slice2labels->SetComputeFeretDiameter( 0 );

  typedef itk::LabelContourImageFilter<LabelImageType,LabelImageType> LabelContourFilterType;

  SbSROIFilterType::Pointer SbSPrevSliceExtractor = SbSROIFilterType::New();
  SbSPrevSliceExtractor->SetInput(ccSbSImageCopy);
  InputImageType::IndexType SbSPrevIndex;
  SbSPrevIndex.Fill(0);
  InputImageType::RegionType SbSPrevRegion;
  SbSPrevRegion.SetSize( SbSSize );
 
  typedef itk::ImageLinearConstIteratorWithIndex<LabelImageType> LinearIterator;

  typedef itk::BinaryErodeImageFilter <LabelImageType, LabelImageType, StructuringElementType> BinaryErodeImageFilterType;

  typedef itk::AddImageFilter<LabelImageType, LabelImageType> AddImageFilterType;

  ConnectedComponentImageFilterType::Pointer finalCCFilter = ConnectedComponentImageFilterType::New();
  RelabelImageFilterType::Pointer finalRelabelFilter = RelabelImageFilterType::New();

  LabelImageType::PixelType pixelValue;
  LabelImageType::PixelType minDist;
  LabelImageType::PixelType prevIdx;  
  LabelImageType::IndexType minDistIdx;
  bool firstIteration;

  MinMaxCalculatorType::Pointer minMaxCalculator = MinMaxCalculatorType::New();
  DuplicatorFilterType::Pointer extractorPrevDuplicatorFilter = DuplicatorFilterType::New();
  LabelContourFilterType::Pointer labelContourFilter = LabelContourFilterType::New();
  SbSROIFilterType::Pointer connectionExtractor = SbSROIFilterType::New();
  BinaryErodeImageFilterType::Pointer connectionErodeFilter = BinaryErodeImageFilterType::New();
  ConnectedComponentImageFilterType::Pointer connectionCCFilter = ConnectedComponentImageFilterType::New();
  BinaryDilateImageFilterType::Pointer connectionDilateFilterOne = BinaryDilateImageFilterType::New();
  BinaryDilateImageFilterType::Pointer connectionDilateFilterTwo = BinaryDilateImageFilterType::New();
  AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  BinaryDilateImageFilterType::Pointer connectionDilateFilter = BinaryDilateImageFilterType::New();
  PasteImageFilterType::Pointer maskPasteFilter = PasteImageFilterType::New();
  DuplicatorFilterType::Pointer extractedImageDuplicator = DuplicatorFilterType::New();
  PasteImageFilterType::Pointer correctedPasteFilter = PasteImageFilterType::New();
  
  PasteImageFilterType::Pointer connCompImagePasteFilter = PasteImageFilterType::New();

  for(unsigned int i = 40; i < imageSize; ++i)
  {
    ccSbSDuplicatorFilter->SetInputImage(ccSbSImage);
    ccSbSDuplicatorFilter->Update();
    ccSbSImageCopy = ccSbSDuplicatorFilter->GetOutput();

    SbSIndex[2] = i;
    SbSRegion.SetIndex( SbSIndex );
    SbSExtractor->SetRegionOfInterest( SbSRegion );
    SbSExtractor->SetInput(ccSbSImageCopy);		 
    SbSExtractor->Update();

    LabelImageType::Pointer extractedImage = SbSExtractor->GetOutput();

    minMaxCalculator->SetImage(extractedImage);
    minMaxCalculator->ComputeMaximum();
    LabelImageType::PixelType maxValue = minMaxCalculator->GetMaximum();

    firstIteration = 1;

    LabelImageType::Pointer extractedSbSImageCopy;
    LabelImageType::Pointer prevSlice;

    if( maxValue == 1 )
    {
      slice2labels->SetInput( extractedImage );
      slice2labels->Update();
      ShapeLabelObjectType::Pointer labelObject = slice2labels->GetOutput()->GetNthLabelObject(0);
      if( labelObject->GetNumberOfPixels() > 10000 || labelObject->GetPerimeter() > 15000 )
      {
        while( maxValue == 1 )
        {
          if( firstIteration )
          {
            SbSPrevIndex[2] = i-1;
            SbSPrevRegion.SetIndex( SbSPrevIndex );
            SbSPrevSliceExtractor->SetRegionOfInterest( SbSPrevRegion );

            ccSbSDuplicatorFilter->SetInputImage( ccSbSImage );
            ccSbSDuplicatorFilter->Update();
            ccSbSImageCopy = ccSbSDuplicatorFilter->GetOutput();

            SbSPrevSliceExtractor->SetInput( ccSbSImageCopy );
            SbSPrevSliceExtractor->Update();
            prevSlice = SbSPrevSliceExtractor->GetOutput();
          }

          extractorPrevDuplicatorFilter->SetInputImage(prevSlice);
          extractorPrevDuplicatorFilter->Update();

          labelContourFilter->SetInput(extractorPrevDuplicatorFilter->GetOutput());
          labelContourFilter->SetFullyConnected(0);
          labelContourFilter->Update();
          LabelImageType::Pointer contoursImage = labelContourFilter->GetOutput();      

          pixelValue = 1;
          minDist = itk::NumericTraits< float >::max();;
          prevIdx = 0;  
          minDistIdx.Fill(0); 

          LinearIterator linIt(contoursImage, contoursImage->GetRequestedRegion());
          linIt.SetDirection(0);
          linIt.GoToBegin();

          while( !linIt.IsAtEnd() )
          {
            while( !linIt.IsAtEndOfLine() )
            {
              if(linIt.Get() > 0 && linIt.Get() != pixelValue)
              {
                // compute distance               
                if( prevIdx != 0 )
                {
                  if( linIt.GetIndex()[0] - prevIdx != 0 && linIt.GetIndex()[0] - prevIdx < minDist )
                  {
                    minDist = linIt.GetIndex()[0] - prevIdx;
                    minDistIdx[0] = prevIdx + int(minDist/2);
                    minDistIdx[1] = linIt.GetIndex()[1];
                   } 
                 }
                 prevIdx = linIt.GetIndex()[0];
                 pixelValue = linIt.Get();
               }
               else if( linIt.Get() > 0 )
               {
                 pixelValue = linIt.Get();
                 //update prev index
                 prevIdx = linIt.GetIndex()[0];
               }
               ++linIt;
             }
             pixelValue = 1;
             prevIdx = 0;
             linIt.NextLine();
          }

          InputImageType::SizeType connectionSize;
          connectionSize[0] = 50;
          connectionSize[1] = 80;
          connectionSize[2] = 1;

          InputImageType::IndexType connectionIndex;
          connectionIndex[0] = minDistIdx[0]-25;
          connectionIndex[1] = minDistIdx[1]-40;
          connectionIndex[2] = 0;

          InputImageType::RegionType connectionRegion;
          connectionRegion.SetSize( connectionSize );
          connectionRegion.SetIndex( connectionIndex );  
 
          connectionExtractor->SetInput(extractedImage);
          connectionExtractor->SetRegionOfInterest( connectionRegion );		 
          connectionExtractor->Update(); 
          LabelImageType::Pointer connectionImage = connectionExtractor->GetOutput();

          // Connection erosion
          StructuringElementType connectionErodeStructElement;
          StructuringElementType::SizeType connectionErodeRadius;
          unsigned int erodeRadiusSize = 3;
          connectionErodeRadius.Fill( erodeRadiusSize );
          connectionErodeRadius[2] = 1;
          connectionErodeStructElement.SetRadius( connectionErodeRadius );
          connectionErodeStructElement.CreateStructuringElement();

          connectionErodeFilter->SetKernel(connectionErodeStructElement);
          connectionErodeFilter->SetForegroundValue( 1 );
          connectionErodeFilter->SetInput(connectionExtractor->GetOutput());
          connectionErodeFilter->Update();

          LabelImageType::Pointer erodeImage = connectionErodeFilter->GetOutput();

          connectionCCFilter->SetInput(erodeImage);
          connectionCCFilter->Update();
            
          while( connectionCCFilter->GetObjectCount() == 1 )
          {
            connectionErodeRadius[0] += 1;
            connectionErodeRadius[1] += 1;        
            connectionErodeRadius[2] = 1;
            connectionErodeStructElement.SetRadius( connectionErodeRadius );
            connectionErodeStructElement.CreateStructuringElement();
            connectionErodeFilter->Update();
            erodeImage = connectionErodeFilter->GetOutput();
            connectionCCFilter->SetInput(erodeImage);
            connectionCCFilter->Update();
          }

          // Connection dilation
          StructuringElementType connectionDilateStructElement;
          StructuringElementType::SizeType connectionDilateRadius;
          unsigned int dilateRadiusSize = erodeRadiusSize;
          connectionDilateRadius.Fill( dilateRadiusSize );
          connectionDilateRadius[2] = 1;
          connectionDilateStructElement.SetRadius( connectionDilateRadius );
          connectionDilateStructElement.CreateStructuringElement();

          connectionDilateFilterOne->SetForegroundValue( 1 );
          connectionDilateFilterOne->SetKernel(connectionDilateStructElement);
          connectionDilateFilterOne->SetInput(connectionCCFilter->GetOutput());
          connectionDilateFilterOne->Update();
          LabelImageType::Pointer dilateImageOne = connectionDilateFilterOne->GetOutput();
        
          connectionDilateFilterTwo->SetForegroundValue( 2 );
          connectionDilateFilterTwo->SetKernel(connectionDilateStructElement);
          connectionDilateFilterTwo->SetInput(connectionCCFilter->GetOutput());
          connectionDilateFilterTwo->Update();      

          LabelImageType::Pointer dilateImageTwo = connectionDilateFilterTwo->GetOutput();
   
          addImageFilter->SetInput1(dilateImageOne);
          addImageFilter->SetInput2(dilateImageTwo);
          addImageFilter->Update();
 
          LabelImageType::Pointer sumImage = addImageFilter->GetOutput();

          connectionDilateFilter->SetForegroundValue( 2 );
          connectionDilateFilter->SetKernel(connectionDilateStructElement);
          connectionDilateFilter->SetInput(connectionDilateFilterOne->GetOutput());
          connectionDilateFilter->Update();
 
          LabelImageType::Pointer dilateImage = connectionDilateFilter->GetOutput();

          LinearIterator sumIt(sumImage, sumImage->GetRequestedRegion());
          sumIt.SetDirection(0);
          sumIt.GoToBegin();
          while( !sumIt.IsAtEnd() )
          {
            while( !sumIt.IsAtEndOfLine() )
            {
              if(sumIt.Get() == 3)
              {
                dilateImage->SetPixel(sumIt.GetIndex(), 0);
              }
              ++sumIt;
            }
            sumIt.NextLine();
          }
           
          pixelValue = 0;
          prevIdx = 0;
          LinearIterator connIt(dilateImage, dilateImage->GetRequestedRegion());
          connIt.SetDirection(0);
          connIt.GoToBegin();
          while( !connIt.IsAtEnd() )
          {
            while( !connIt.IsAtEndOfLine() )
            {
              if(pixelValue > 0 && connIt.Get() > 0 && connIt.Get() != pixelValue)
              {
                LabelImageType::IndexType prIdx = connIt.GetIndex();
                prIdx[0] = prevIdx;
                dilateImage->SetPixel(prIdx, 0);               
                dilateImage->SetPixel(connIt.GetIndex(), 0);
               }
               pixelValue = connIt.Get();
               prevIdx = connIt.GetIndex()[0];
               ++connIt;
             }
             pixelValue = 0;
             prevIdx = 0;
             connIt.NextLine();
          }

          ThresholdImageFilterType::Pointer threshFilter = ThresholdImageFilterType::New();
          threshFilter->SetInput(dilateImage); 
          threshFilter->ThresholdOutside(0,2);
          threshFilter->SetOutsideValue(0);
          threshFilter->Update();

          dilateImage = threshFilter->GetOutput();

          PasteImageFilterType::Pointer dilatePasteFilter = PasteImageFilterType::New();
          dilatePasteFilter->SetSourceImage(dilateImage);                                                 
          dilatePasteFilter->SetDestinationImage(extractedImage);                                               
          dilatePasteFilter->SetSourceRegion(dilateImage->GetLargestPossibleRegion());   
          dilatePasteFilter->SetDestinationIndex(connectionIndex);

          finalCCFilter->SetInput(dilatePasteFilter->GetOutput());
          finalCCFilter->Update();

          extractedImage = finalCCFilter->GetOutput();

          minMaxCalculator->SetImage(extractedImage);
          minMaxCalculator->ComputeMaximum();
          maxValue = minMaxCalculator->GetMaximum();

          if( maxValue == 1 )
          {
            connectionImage->FillBuffer(itk::NumericTraits< LabelPixelType >::Zero);
            maskPasteFilter->SetSourceImage(connectionImage);                                                 
            maskPasteFilter->SetDestinationImage(prevSlice);                                               
            maskPasteFilter->SetSourceRegion(connectionImage->GetLargestPossibleRegion());   
            maskPasteFilter->SetDestinationIndex(connectionIndex);
            maskPasteFilter->Update();
            prevSlice = maskPasteFilter->GetOutput();
          }
 
          connectionIndex.Fill(0);
          connectionIndex[2] = i;

          extractedImageDuplicator->SetInputImage(extractedImage);
          extractedImageDuplicator->Update();
          LabelImageType::Pointer extractedImageCopy = extractedImageDuplicator->GetOutput();
          
          ccSbSDuplicatorFilter->Update();
          ccSbSImageCopy = ccSbSDuplicatorFilter->GetOutput();

          correctedPasteFilter->SetSourceImage(extractedImageCopy);                                                 
          correctedPasteFilter->SetDestinationImage(ccSbSImageCopy);                                               
          correctedPasteFilter->SetSourceRegion(extractedImageCopy->GetLargestPossibleRegion());   
          correctedPasteFilter->SetDestinationIndex(connectionIndex);
          correctedPasteFilter->Update();

          ccSbSImage = correctedPasteFilter->GetOutput();

          connCompImagePasteFilter->SetSourceImage(extractedImageCopy);                                                 
          connCompImagePasteFilter->SetDestinationImage(connCompImage);                                               
          connCompImagePasteFilter->SetSourceRegion(extractedImageCopy->GetLargestPossibleRegion());   
          connCompImagePasteFilter->SetDestinationIndex(connectionIndex);
          connCompImagePasteFilter->Update();
          connCompImage = connCompImagePasteFilter->GetOutput();

          firstIteration = 0;
        }
      }
    }
  }

  ConnectedComponentImageFilterType::Pointer ConnCompFilter = ConnectedComponentImageFilterType::New();
  ConnCompFilter->SetInput(connCompImage);
  RelabelImageFilterType::Pointer relFilter = RelabelImageFilterType::New();
  relFilter->SetInput(ConnCompFilter->GetOutput());
  relFilter->SetMinimumObjectSize(5000);
  relFilter->Update();

  connCompImage = relFilter->GetOutput();

  typedef itk::SliceBySliceImageFilter<LabelImageType, LabelImageType> ClosingSliceBySliceFilterType;
  typedef itk::BinaryBallStructuringElement< LabelPixelType, Dim-1 > lungsStructuringElementType;  	
  lungsStructuringElementType lungsStructElement;
  lungsStructuringElementType::SizeType lungsRadius;
  lungsRadius.Fill( 15 );
  lungsStructElement.SetRadius( lungsRadius );
  lungsStructElement.CreateStructuringElement();
	
  typedef itk::BinaryMorphologicalClosingImageFilter < ClosingSliceBySliceFilterType::InternalInputImageType, ClosingSliceBySliceFilterType::InternalOutputImageType, lungsStructuringElementType > lungsClosingFilterType;
  lungsClosingFilterType::Pointer finalClosingFilter = lungsClosingFilterType::New();

  finalClosingFilter->SetKernel( lungsStructElement );

  ClosingSliceBySliceFilterType::Pointer finalClosingSbSFilter = ClosingSliceBySliceFilterType::New();
  finalClosingSbSFilter->SetDimension(2);

  for(unsigned int i = 1; i <= 2; ++i)
  {
    finalClosingSbSFilter->SetInput(connCompImage);
    finalClosingFilter->SetForegroundValue( i );
    finalClosingSbSFilter->SetFilter(finalClosingFilter);
    finalClosingSbSFilter->Update();

    connCompImage = finalClosingSbSFilter->GetOutput();
  }
  
  clockSeparation.Stop();
  std::cout << "Lung Separation.: " << clockSeparation.GetTotal() << std::endl;

  //
  // First set all types to 'UNDEFINEDTYPE'. This is necessary in the
  // case that, e.g., airways are present, connecting the left and
  // right lungs
  //

  cip::ChestConventions conventions;
  
  typedef itk::ImageRegionIteratorWithIndex< LabelImageType > LabelMapIteratorType;
  LabelMapIteratorType mIt( connCompImage, connCompImage->GetBufferedRegion() );

  unsigned char cipRegion;

  mIt.GoToBegin();
  while ( !mIt.IsAtEnd() )
  {
	  if ( mIt.Get() != 0 )
	  {
		  cipRegion = conventions.GetChestRegionFromValue( mIt.Get() );
		  mIt.Set( conventions.GetValueFromChestRegionAndType( cipRegion, static_cast< unsigned char >( cip::UNDEFINEDTYPE ) ) );
	  }
	  ++mIt;
  }

  connectedComponentFilter->SetInput( connCompImage );
  connectedComponentFilter->Update();

  relabelFilter->SetInput( connectedComponentFilter->GetOutput() );
  try
  {
	  relabelFilter->Update();
  }
  catch ( itk::ExceptionObject &excp )
  {
	  std::cerr << "Exception caught relabeling:";
	  std::cerr << excp << std::endl;
  }

  unsigned int total = 0;
  for ( unsigned int i = 0; i < relabelFilter->GetNumberOfObjects(); i++ )
  {
	  total += relabelFilter->GetSizeOfObjectsInPixels()[i];
  }

  LabelImageType::Pointer relabeledLungLabelMap = relabelFilter->GetOutput();
  //
  // If we're here, we assume that the left and right have been
  // separated, so label them. First, we need to get the relabel
  // component corresponding to the left and the right. We assume that
  // the relabel component value = 1 corresponds to one of the two
  // lungs and a value of 2 corresponds to the other. Find the
  // left-most and right-most component value. Assuming the scan is
  // supine, head-first, the component value corresponding to the
  // smallest x-index will be the left lung and the other major
  // component will be the right lung.
  //
  unsigned int minX = relabeledLungLabelMap->GetBufferedRegion().GetSize()[0];
  unsigned int maxX = 0;

  unsigned int smallIndexComponentLabel, largeIndexComponentLabel;

  LabelMapIteratorType rIt( relabeledLungLabelMap, relabeledLungLabelMap->GetBufferedRegion() );

  rIt.GoToBegin();
  while ( !rIt.IsAtEnd() )
  {
	  if ( rIt.Get() == 1 || rIt.Get() == 2 )
	  {
		  if ( rIt.GetIndex()[0] < minX )
		  {
			  smallIndexComponentLabel = rIt.Get();
			  minX = rIt.GetIndex()[0];
		  }
		  if ( rIt.GetIndex()[0] > maxX )
		  {
			  largeIndexComponentLabel = rIt.Get();
			  maxX = rIt.GetIndex()[0];
		  }
	  }	  
	  ++rIt;
  }

  unsigned int leftLungComponentLabel, rightLungComponentLabel;
  //if ( (this->HeadFirst && this->Supine) || (this->FeetFirst && this->Prone) )
    //{
    leftLungComponentLabel  = largeIndexComponentLabel;
    rightLungComponentLabel = smallIndexComponentLabel;
    //}
  //else
    //{
    //leftLungComponentLabel  = smallIndexComponentLabel;
    //rightLungComponentLabel = largeIndexComponentLabel;
    //}

  mIt.GoToBegin();
  rIt.GoToBegin();
  while ( !mIt.IsAtEnd() )
  {
	  if ( rIt.Get() == leftLungComponentLabel )
	  {
		  mIt.Set( static_cast< unsigned short >( cip::LEFTLUNG ) );
	  }
	  if ( rIt.Get() == rightLungComponentLabel )
      {
		  mIt.Set( static_cast< unsigned short >( cip::RIGHTLUNG ) );
      }
	  ++rIt;
	  ++mIt;
  }

  //
  // Get the number of voxels in the label map
  //
  unsigned int totalVoxelCount = 0;
 
  mIt.GoToBegin();
  while ( !mIt.IsAtEnd() )
  {
	  if ( mIt.Get() != 0 )
	  {
		  totalVoxelCount++;
      }
	  ++mIt;
  }

  //
  // Label by thirds the lung label map
  //
  bool foundLeftLung  = false;
  bool foundRightLung = false;
  unsigned int voxelCount = 0;

  mIt.GoToBegin();
  while ( !mIt.IsAtEnd() )
  {
	  if ( mIt.Get() != 0 )
	  {
		  voxelCount++;
		  
		  if ( static_cast< double >( voxelCount ) < static_cast< double >( totalVoxelCount )/3.0 )
		  {
			  if ( mIt.Get() == static_cast< unsigned short >( cip::LEFTLUNG ) )
			  {
				  foundLeftLung = true;
				  mIt.Set( static_cast< unsigned short >( cip::LEFTLOWERTHIRD ) );
			  }
			  else
			  {
				  foundRightLung = true;
				  mIt.Set( static_cast< unsigned short >( cip::RIGHTLOWERTHIRD ) );
			  }
		  }
		  else if ( static_cast< double >( voxelCount ) < 2.0*static_cast< double >( totalVoxelCount )/3.0 )
		  {
			  if ( mIt.Get() == static_cast< unsigned short >( cip::LEFTLUNG ) )
			  {
				  foundLeftLung = true;
				  mIt.Set( static_cast< unsigned short >( cip::LEFTMIDDLETHIRD ) );
			  }
			  else
			  {
				  foundRightLung = true;
				  mIt.Set( static_cast< unsigned short >( cip::RIGHTMIDDLETHIRD ) );
			  }
		  }
		  else
		  {
			  if ( mIt.Get() == static_cast< unsigned short >( cip::LEFTLUNG ) )
			  {
				  foundLeftLung = true;
				  mIt.Set( static_cast< unsigned short >( cip::LEFTUPPERTHIRD ) );
			  }
			  else
			  {
				  foundRightLung = true;
				  mIt.Set( static_cast< unsigned short >( cip::RIGHTUPPERTHIRD ) );
			  }
		  }
	  }
	  ++mIt;
  }   

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( OutputVolume.c_str() );
  writer->SetInput( connCompImage );
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
