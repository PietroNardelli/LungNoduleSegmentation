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

#include "itkImageDuplicator.h"

#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"

#include "itkTimeProbe.h"

#include "LungSegmentationCLICLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const unsigned int Dim = 3;

  typedef signed short							InputPixelType;
  typedef unsigned short						OutputPixelType;
  typedef unsigned int 							LabelType;
  typedef float								MSOutputPixelType;

  typedef itk::Image< InputPixelType,  Dim > 				InputImageType;
  typedef itk::Image< OutputPixelType, Dim > 				OutputImageType;
  typedef itk::Image< LabelType,       Dim > 				LabelImageType;
  typedef itk::Image< MSOutputPixelType, Dim > 				MSOutputImageType;

  typedef itk::ImageFileReader<InputImageType>  			ReaderType;

  typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > 	I2LType;
  typedef I2LType::OutputImageType 					LabelMapType;
  typedef LabelMapType::LabelObjectType					ShapeLabelObjectType;

  itk::TimeProbe lungSegmentationClock;
  lungSegmentationClock.Start();

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( InputVolume.c_str() );
  reader->Update();
  InputImageType::Pointer inputImage = reader->GetOutput();  

  typedef itk::ImageDuplicator< InputImageType >InputDuplicatorType;
  InputDuplicatorType::Pointer inputDuplicator = InputDuplicatorType::New();
  inputDuplicator->SetInputImage(inputImage);
  inputDuplicator->Update();
  InputImageType::Pointer clonedInput = inputDuplicator->GetOutput();

  // Binarize the image based on the specified threshold
  typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > thresholdFilterType;
  thresholdFilterType::Pointer thresholdFilter = thresholdFilterType::New();
  thresholdFilter->SetInput( 0, inputImage );
  thresholdFilter->SetInsideValue( lungsColor );
  thresholdFilter->SetOutsideValue( 0 );
  thresholdFilter->SetLowerThreshold(Lower);
  thresholdFilter->SetUpperThreshold(Upper);
  
  // Remove holes not connected to the boundary of the image.
  typedef itk::BinaryFillholeImageFilter< OutputImageType > FillHolesFilterType;
  FillHolesFilterType::Pointer holeFillingFilter = FillHolesFilterType::New();
  holeFillingFilter->SetInput( thresholdFilter->GetOutput() );
  holeFillingFilter->SetForegroundValue( lungsColor );
  holeFillingFilter->SetFullyConnected( 1 );
 
  // Remove small (i.e., smaller than the structuring element) holes and tube like structures in the interior or at the boundaries of the image.
  typedef itk::BinaryBallStructuringElement< OutputPixelType, Dim > StructuringElementType;  	
  StructuringElementType structElement;
  StructuringElementType::SizeType radius;
  radius.Fill( 1 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
	
  typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > ClosingFilterType;
  ClosingFilterType::Pointer closing = ClosingFilterType::New();
  closing->SetInput( holeFillingFilter->GetOutput() );
  closing->SetKernel( structElement );
  closing->SetForegroundValue( lungsColor );

  // Convert the binary image to a label map to valuate the shape attributes 
  I2LType::Pointer i2l = I2LType::New();
  i2l->SetInput( closing->GetOutput() );
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
  labelToBinaryFilter->GetOutput();

  radius.Fill( 20 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
  closing->SetInput( labelToBinaryFilter->GetOutput() );
  closing->SetKernel( structElement );
  closing->SetForegroundValue( lungsColor );
  closing->SetSafeBorder( 1 );
  
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( OutputVolume.c_str() );
  writer->SetInput( closing->GetOutput() );
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

  lungSegmentationClock.Stop();
  std::cout << "Total time for lung segmentation: " << lungSegmentationClock.GetTotal() << std::endl;

  return EXIT_SUCCESS;
}
