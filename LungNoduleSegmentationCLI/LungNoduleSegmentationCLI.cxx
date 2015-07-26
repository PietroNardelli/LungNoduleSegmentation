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

#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"

#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"

#include "itkRegionOfInterestImageFilter.h"

#include "itkNumericTraits.h"
#include "itkMergeLabelMapFilter.h"

#include "itkMaskImageFilter.h"

#include "itkMultiScaleGaussianEnhancementImageFilter.h"
#include "itkGeneralEigenValueFunctor.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "itkConnectedThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkTimeProbe.h"

#include "LungNoduleSegmentationCLICLP.h"

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

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( InputVolume.c_str() );

  itk::TimeProbe clock;
  // Binarize the image based on the specified threshold
  clock.Start();
  typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > thresholdFilterType;
  thresholdFilterType::Pointer thresholdFilter = thresholdFilterType::New();
  thresholdFilter->SetInput( 0, reader->GetOutput() );
  thresholdFilter->SetInsideValue( 1 );
  thresholdFilter->SetOutsideValue( 0 );
  thresholdFilter->SetLowerThreshold(Lower);
  thresholdFilter->SetUpperThreshold(Upper);

  // Remove holes not connected to the boundary of the image.
  typedef itk::BinaryFillholeImageFilter< OutputImageType > FillHolesFilterType;
  FillHolesFilterType::Pointer holeFillingFilter = FillHolesFilterType::New();
  holeFillingFilter->SetInput( thresholdFilter->GetOutput() );
  holeFillingFilter->SetForegroundValue( 1 );
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
  closing->SetForegroundValue( 1 );
  closing->Update();

  // Convert the binary image to a label map to valuate the shape attributes 
  I2LType::Pointer i2l = I2LType::New();
  i2l->SetInput( closing->GetOutput() );
  i2l->SetInputForegroundValue( 1 );
  i2l->SetFullyConnected( 1 );
  i2l->Update();

  LabelMapType::Pointer labelMap = i2l->GetOutput();
  std::cout << "Number of objects beforehand: " << labelMap->GetNumberOfLabelObjects() << std::endl;
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

  std::cout << "Number of objects afterwards: " << labelMap->GetNumberOfLabelObjects() << std::endl;

  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, OutputImageType > LabelMapToBinaryImageType; 
  LabelMapToBinaryImageType::Pointer labelToBinaryFilter = LabelMapToBinaryImageType::New();
  labelToBinaryFilter->SetInput( labelMap );
  labelToBinaryFilter->SetBackgroundValue( 0 );
  labelToBinaryFilter->SetForegroundValue( 1 );

  radius.Fill( 7 );
  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
  closing->SetInput( labelToBinaryFilter->GetOutput() );
  closing->SetKernel( structElement );
  closing->SetForegroundValue( 1 );
  closing->SetSafeBorder( 1 );
  closing->Update();

  // Mask the input image with the lung label
  typedef itk::MaskImageFilter< InputImageType, OutputImageType, InputImageType > MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  maskFilter->SetInput(reader->GetOutput());
  maskFilter->SetMaskImage(closing->GetOutput());
  maskFilter->SetOutsideValue(-1024);

  clock.Stop();
  std::cout << "Total time for lung segmentation: " << clock.GetTotal() << std::endl;

  // Extract a ROI containing the nodule
  
  double roi[6];
  clock.Start();
  for (unsigned int i = 0; i < Dim; i++)
  {
     roi[2*i]= seed[0][i] - 20;
     roi[2*i+1]= seed[0][i] + 20;
  }

  roi[0] = roi[0] * (-reader->GetOutput()->GetDirection()[0][0]);
  roi[1] = roi[1] * (-reader->GetOutput()->GetDirection()[0][0]);
  roi[2] = roi[2] * (-reader->GetOutput()->GetDirection()[1][1]);
  roi[3] = roi[3] * (-reader->GetOutput()->GetDirection()[1][1]);
  roi[4] = roi[4] *   reader->GetOutput()->GetDirection()[2][2];
  roi[5] = roi[5] *   reader->GetOutput()->GetDirection()[2][2];

  // Convert bounds into region indices
  OutputImageType::PointType p1, p2;
  OutputImageType::IndexType pi1, pi2;
  OutputImageType::IndexType startIndex;
  for (unsigned int i = 0; i < Dim; i++)
  {
    p1[i] = roi[2*i];
    p2[i] = roi[2*i+1];
  }

  reader->GetOutput()->TransformPhysicalPointToIndex(p1, pi1);
  reader->GetOutput()->TransformPhysicalPointToIndex(p2, pi2);

  OutputImageType::SizeType roiSize;
  for (unsigned int i = 0; i < Dim; i++)
    {
    roiSize[i] = fabs(pi2[i] - pi1[i]);
    startIndex[i] = (pi1[i]<pi2[i])?pi1[i]:pi2[i];
    }

  OutputImageType::RegionType roiRegion( startIndex, roiSize );

  typedef itk::RegionOfInterestImageFilter< InputImageType, InputImageType > ROIFilterType;

  ROIFilterType::Pointer inputROIFilter = ROIFilterType::New();
  inputROIFilter->SetRegionOfInterest( roiRegion );
  inputROIFilter->SetInput( reader->GetOutput() );

  // Enhance blob-like structures in the ROI
  typedef itk::MultiScaleGaussianEnhancementImageFilter<InputImageType,MSOutputImageType> MultiScaleFilterType;
  MultiScaleFilterType::Pointer multiScaleFilter = MultiScaleFilterType::New();

  multiScaleFilter->SetInput( inputROIFilter->GetOutput() );
  multiScaleFilter->SetSigmaMinimum( 1 );
  multiScaleFilter->SetSigmaMaximum( 6 );
  multiScaleFilter->SetNumberOfSigmaSteps( 5 );
  multiScaleFilter->SetNonNegativeHessianBasedMeasure( true );
  multiScaleFilter->SetSigmaStepMethod( 1 );
  multiScaleFilter->SetRescale( 0 );

  typedef itk::Functor::GeneralEigenValueFunctor<MultiScaleFilterType::EigenValueArrayType, MSOutputPixelType > FunctorType;
  FunctorType::Pointer functor = FunctorType::New();

  multiScaleFilter->SetUnaryFunctor( functor );

  // Find the max value in the enhanced image
  typedef itk::MinimumMaximumImageFilter<MSOutputImageType> MinMaxFilterType;
  MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
  minMaxFilter->SetInput(multiScaleFilter->GetOutput());
  minMaxFilter->Update();

  // Separate nodule candidate region from other structures
  typedef itk::BinaryThresholdImageFilter< MSOutputImageType, OutputImageType > noduleThresholdFilterType;
  noduleThresholdFilterType::Pointer noduleThresholdFilter = noduleThresholdFilterType::New();
  noduleThresholdFilter->SetInput( 0, multiScaleFilter->GetOutput() );
  noduleThresholdFilter->SetInsideValue( 1 );
  noduleThresholdFilter->SetOutsideValue( 0 );
  noduleThresholdFilter->SetLowerThreshold(40);
  noduleThresholdFilter->SetUpperThreshold(minMaxFilter->GetMaximum()+1);

  // Exclude objects with a Feret diameter smaller than 2.5 mm 
  I2LType::Pointer nodule2label = I2LType::New();
  nodule2label->SetInput( noduleThresholdFilter->GetOutput() );
  nodule2label->SetInputForegroundValue( 1 );
  nodule2label->SetFullyConnected( 1 );
  nodule2label->SetComputeFeretDiameter( 1 );
  nodule2label->Update();

  LabelMapType::Pointer noduleLabelMap = nodule2label->GetOutput();
  std::cout << "Number of objects in the nodule label map: " << noduleLabelMap->GetNumberOfLabelObjects() << std::endl;
  std::vector<unsigned long> noduleLabelsToRemove;

  for(unsigned int n = 0; n < noduleLabelMap->GetNumberOfLabelObjects(); ++n)
  {
     ShapeLabelObjectType::Pointer noduleLabelObject = noduleLabelMap->GetNthLabelObject(n);
     if( noduleLabelObject->GetFeretDiameter() <= 2.5 )
     {
       noduleLabelsToRemove.push_back(noduleLabelObject->GetLabel());
     }
  }

  for(unsigned int i = 0; i < noduleLabelsToRemove.size(); ++i)
  {
      noduleLabelMap->RemoveLabel(noduleLabelsToRemove[i]);
  }

  // Convert nodule label map back into a binary image
  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, OutputImageType > noduleLabelMapToBinaryImageType; 
  noduleLabelMapToBinaryImageType::Pointer noduleLabelToBinaryFilter = noduleLabelMapToBinaryImageType::New();
  noduleLabelToBinaryFilter->SetInput( noduleLabelMap );
  noduleLabelToBinaryFilter->SetBackgroundValue( 0 );
  noduleLabelToBinaryFilter->SetForegroundValue( 1 );

  clock.Stop();
  std::cout << "Total time to enhance the ROI: " << clock.GetTotal() << std::endl;

  // Use the ROI without the pleura
  clock.Start();
  ROIFilterType::Pointer maskedImageROIFilter = ROIFilterType::New();
  maskedImageROIFilter->SetRegionOfInterest( roiRegion );
  maskedImageROIFilter->SetInput( maskFilter->GetOutput() );
  maskedImageROIFilter->Update();
 
  InputImageType::Pointer roiImage = InputImageType::New();
  roiImage = maskedImageROIFilter->GetOutput();

  // Exclude all the voxels with an intensity smaller than -824 HU from the following analysis and segmentation 
  thresholdFilterType::Pointer roiThresholdFilter = thresholdFilterType::New();
  roiThresholdFilter->SetInput( 0, roiImage );
  roiThresholdFilter->SetInsideValue( 1 );
  roiThresholdFilter->SetOutsideValue( 0 );
  roiThresholdFilter->SetLowerThreshold(-824);
  roiThresholdFilter->SetUpperThreshold(5000);

  // And filtering to exclude voxels with an intensity smaller than -824 HU (in case of airways within the nodule)
  typedef itk::AndImageFilter< OutputImageType > andFilterType;
  andFilterType::Pointer andFilter = andFilterType::New();
  andFilter->SetInput(0, noduleLabelToBinaryFilter->GetOutput());
  andFilter->SetInput(1, roiThresholdFilter->GetOutput());

  // Compute mean and std. dev. of the region underneath the image label
  typedef itk::LabelStatisticsImageFilter< InputImageType, OutputImageType > LabelStatisticsImageFilterType;
  LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetLabelInput( andFilter->GetOutput() );
  labelStatisticsImageFilter->SetInput(roiImage);
  labelStatisticsImageFilter->Update();

  LabelStatisticsImageFilterType::RealType mean = labelStatisticsImageFilter->GetMean( 1 );
  LabelStatisticsImageFilterType::RealType sigma = labelStatisticsImageFilter->GetSigma( 1 );

  LabelStatisticsImageFilterType::RealType upperThreshold = mean + 2*sigma;
  LabelStatisticsImageFilterType::RealType lowerThreshold = mean - 2*sigma;

  /*std::cout << "mean: " << mean << std::endl;
  std::cout << "sigma: " << sigma << std::endl;
  std::cout << "range: " << lowerThreshold << " -> " << upperThreshold << std::endl;*/

  // Segment the nodule
  typedef itk::ConnectedThresholdImageFilter< InputImageType, OutputImageType > ConnectedThresholdFilterType;
  ConnectedThresholdFilterType::Pointer connectedThresholdFilter = ConnectedThresholdFilterType::New();
  connectedThresholdFilter->SetInput( roiImage );			  
  connectedThresholdFilter->SetReplaceValue( noduleColor ); 	
  connectedThresholdFilter->SetUpper( upperThreshold );
  connectedThresholdFilter->SetLower( lowerThreshold ); 

  // Seed comes in ras, convert to lps
  InputImageType::PointType lpsPoint;
  InputImageType::IndexType index;

  lpsPoint[0] = seed[0][0] * (-roiImage->GetDirection()[0][0]);
  lpsPoint[1] = seed[0][1] * (-roiImage->GetDirection()[1][1]);
  lpsPoint[2] = seed[0][2] *   roiImage->GetDirection()[2][2];

  roiImage->TransformPhysicalPointToIndex(lpsPoint, index);
  connectedThresholdFilter->AddSeed( index );
  connectedThresholdFilter->Update();

  OutputImageType::Pointer noduleLabel = OutputImageType::New();
  noduleLabel = connectedThresholdFilter->GetOutput();

  // If the starting threshold gives an empty segmentation, neighbours are checked
  typedef itk::StatisticsImageFilter< OutputImageType > StatisticsImageFilterType;
  StatisticsImageFilterType::Pointer statisticsFilter = StatisticsImageFilterType::New();

  statisticsFilter->SetInput(noduleLabel);
  statisticsFilter->Update();
  double n_voxels = ( statisticsFilter->GetSum() );

  OutputImageType::SizeType   neighborhoodRadius, neighborhoodRegionSize;
  OutputImageType::IndexType  neighborhoodRegionIndex;
  OutputImageType::RegionType neighborhoodRegion;	                                                        
    
  if( n_voxels == 0 )
  {
      bool isMinor = 0;
      neighborhoodRegionSize.Fill(3);
      neighborhoodRegionIndex = index;
      neighborhoodRadius.Fill(3);

      neighborhoodRegion.SetSize(neighborhoodRegionSize);
      neighborhoodRegion.SetIndex(neighborhoodRegionIndex);
	
      typedef itk::ConstNeighborhoodIterator< InputImageType > neighborhoodIterator;
      neighborhoodIterator iterator(neighborhoodRadius, roiImage, neighborhoodRegion);
		
      unsigned int counter = 0;

      while( counter < iterator.Size() && !isMinor )
      {
          if( iterator.GetPixel(counter) > lowerThreshold )
          {
                index = iterator.GetIndex( counter );				
			
                connectedThresholdFilter->ClearSeeds();
                connectedThresholdFilter->AddSeed( index );
                connectedThresholdFilter->Update();
				
                noduleLabel = connectedThresholdFilter->GetOutput();	
			
                isMinor = 1;
            }
            counter++;  
        }
        if ( !isMinor )
        {
            std::cout<<"Please move the seed point in a different position."<<std::endl;
            return EXIT_FAILURE;
        }
  }
 
  for( unsigned int i = 0; i < 5; i++)
  {
    // Compute mean and std. dev. of the region underneath the image label
    labelStatisticsImageFilter->SetLabelInput( noduleLabel );
    labelStatisticsImageFilter->SetInput( roiImage );
    labelStatisticsImageFilter->Update();

    mean = labelStatisticsImageFilter->GetMean( noduleColor );
    sigma = labelStatisticsImageFilter->GetSigma( noduleColor );

    upperThreshold = mean + 2*sigma;
    lowerThreshold = mean - 2*sigma;

    /*std::cout << "mean: " << mean << std::endl;
    std::cout << "sigma: " << sigma << std::endl;
    std::cout << "range: " << lowerThreshold << " -> " << upperThreshold << std::endl;*/

    // Segment the nodule
    connectedThresholdFilter->SetUpper( upperThreshold );
    connectedThresholdFilter->SetLower( lowerThreshold ); 
    connectedThresholdFilter->AddSeed( index );
    connectedThresholdFilter->Update();
    noduleLabel = connectedThresholdFilter->GetOutput();
  }

  typedef itk::BinaryBallStructuringElement< OutputPixelType, Dim > finalStructuringElementType;  	
  finalStructuringElementType finalStructElement;
  finalStructuringElementType::SizeType finalRadius;
  finalRadius.Fill( 1 );

  finalStructElement.SetRadius( finalRadius );
  finalStructElement.CreateStructuringElement();
	
  typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > finalClosingFilterType;
  finalClosingFilterType::Pointer finalClosing = finalClosingFilterType::New();
  finalClosing->SetInput( noduleLabel );
  finalClosing->SetKernel( finalStructElement );
  finalClosing->SetForegroundValue( noduleColor );

  try
  {   
        finalClosing->Update();
  }
  catch ( itk::ExceptionObject & e )
  {
    	std::cerr << "exception in closing " << std::endl;
	std::cerr << e.GetDescription() << std::endl;
	std::cerr << e.GetLocation() << std::endl;
	return EXIT_FAILURE;
  }
  clock.Stop();
  std::cout << "Total time for nodule segmentation: " << clock.GetTotal() << std::endl;
  
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( OutputVolume.c_str() );
  writer->SetInput( finalClosing->GetOutput() );
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
