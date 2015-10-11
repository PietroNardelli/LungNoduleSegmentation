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

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"

#include "itkRegionOfInterestImageFilter.h"

#include "itkNumericTraits.h"

#include "itkMaskImageFilter.h"

#include "itkMultiScaleGaussianEnhancementImageFilter.h"
#include "itkGeneralEigenValueFunctor.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include "itkConnectedThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkImageDuplicator.h"

#include "itkLabelStatisticsImageFilter.h"

#include "itkPoint.h"

#include "cipChestConventions.h"

#include "LungNoduleSegmentationCLICLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const unsigned int Dim = 3;

  typedef signed short							InputPixelType;
  typedef unsigned short						OutputPixelType;
  typedef unsigned int 							LabelType;
  typedef float									MSOutputPixelType;

  typedef itk::Image< InputPixelType,  Dim > 				InputImageType;
  typedef itk::Image< OutputPixelType, Dim > 				OutputImageType;
  typedef itk::Image< LabelType,       Dim > 				LabelImageType;
  typedef itk::Image< MSOutputPixelType, Dim > 				MSOutputImageType;

  typedef itk::ImageFileReader<InputImageType>  			InputReaderType;
  typedef itk::ImageFileReader<OutputImageType>  			LungLabelReaderType;

  typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType > 	I2LType;
  typedef I2LType::OutputImageType 					LabelMapType;
  typedef LabelMapType::LabelObjectType					ShapeLabelObjectType;

  InputReaderType::Pointer inputReader = InputReaderType::New();
  inputReader->SetFileName( InputVolume.c_str() );
  inputReader->Update();
  InputImageType::Pointer inputImage = inputReader->GetOutput();  

  LungLabelReaderType::Pointer lungLabelReader = LungLabelReaderType::New();
  lungLabelReader->SetFileName( InputLungLabel.c_str() );
  lungLabelReader->Update();
  OutputImageType::Pointer lungLabelImage = lungLabelReader->GetOutput();

  typedef itk::ImageDuplicator< InputImageType > DuplicatorType;
  DuplicatorType::Pointer inputDuplicator = DuplicatorType::New();
  inputDuplicator->SetInputImage(inputImage);
  inputDuplicator->Update();
  InputImageType::Pointer clonedInput = inputDuplicator->GetOutput();

  // Extract a ROI containing the nodule
  double roi[6];
  for (unsigned int i = 0; i < Dim; i++)
  {
     roi[2*i]= seed[0][i] - 20;
     roi[2*i+1]= seed[0][i] + 20;
  }

  roi[0] = roi[0] * (-inputReader->GetOutput()->GetDirection()[0][0]);
  roi[1] = roi[1] * (-inputReader->GetOutput()->GetDirection()[0][0]);
  roi[2] = roi[2] * (-inputReader->GetOutput()->GetDirection()[1][1]);
  roi[3] = roi[3] * (-inputReader->GetOutput()->GetDirection()[1][1]);
  roi[4] = roi[4] *   inputReader->GetOutput()->GetDirection()[2][2];
  roi[5] = roi[5] *   inputReader->GetOutput()->GetDirection()[2][2];

  // Convert bounds into region indices
  InputImageType::PointType p1, p2;
  InputImageType::IndexType pi1, pi2;
  InputImageType::IndexType startIndex;
  for (unsigned int i = 0; i < Dim; i++)
  {
    p1[i] = roi[2*i];
    p2[i] = roi[2*i+1];
  }

  inputReader->GetOutput()->TransformPhysicalPointToIndex(p1, pi1);
  inputReader->GetOutput()->TransformPhysicalPointToIndex(p2, pi2);

  InputImageType::SizeType roiSize;
  for(unsigned int i = 0; i < Dim; i++)
    {
    roiSize[i] = fabs(static_cast<float>(pi2[i] - pi1[i]));
    startIndex[i] = (pi1[i]<pi2[i])?pi1[i]:pi2[i];
    }

  InputImageType::RegionType roiRegion( startIndex, roiSize );

  typedef itk::RegionOfInterestImageFilter< InputImageType, InputImageType > ROIFilterType;

  ROIFilterType::Pointer inputROIFilter = ROIFilterType::New();
  inputROIFilter->SetRegionOfInterest( roiRegion );
  inputROIFilter->SetInput( inputImage );

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

  //std::cout << "Number of objects in the nodule label map before: " << noduleLabelMap->GetNumberOfLabelObjects() << std::endl;
  std::vector<unsigned long> noduleLabelsToRemove;

  itk::Point< double, Dim > lps;
  lps[0] = seed[0][0] * (-inputReader->GetOutput()->GetDirection()[0][0]);
  lps[1] = seed[0][1] * (-inputReader->GetOutput()->GetDirection()[1][1]);
  lps[2] = seed[0][2] *   inputReader->GetOutput()->GetDirection()[2][2];

  double minDist = itk::NumericTraits< double >::max();
  unsigned long minLabel = noduleLabelMap->GetNthLabelObject(0)->GetLabel();

  for(unsigned int n = 0; n < noduleLabelMap->GetNumberOfLabelObjects(); ++n)
  {
     ShapeLabelObjectType::Pointer noduleLabelObject = noduleLabelMap->GetNthLabelObject(n);
     double currentDist = noduleLabelObject->GetCentroid().EuclideanDistanceTo(lps);
     if( noduleLabelObject->GetFeretDiameter() > 1 && currentDist < minDist )
     {
		 minDist = currentDist;
        //std::cout<<minDist<<std::endl;
        if( n != 0 )
        {
          noduleLabelsToRemove.push_back(minLabel);
        }
        minLabel = noduleLabelObject->GetLabel();
     }
     else
     {
       noduleLabelsToRemove.push_back(noduleLabelObject->GetLabel());
     }
  
     //if( noduleLabelObject->GetFeretDiameter() <= 2.5 )
     //{     
       //noduleLabelsToRemove.push_back(noduleLabelObject->GetLabel());
     //}
  }

  for(unsigned int i = 0; i < noduleLabelsToRemove.size(); ++i)
  {
      noduleLabelMap->RemoveLabel(noduleLabelsToRemove[i]);
  }

  //std::cout << "Number of objects in the nodule label map after: " << noduleLabelMap->GetNumberOfLabelObjects() << std::endl;

  // Convert nodule label map back into a binary image
  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, OutputImageType > noduleLabelMapToBinaryImageType; 
  noduleLabelMapToBinaryImageType::Pointer noduleLabelToBinaryFilter = noduleLabelMapToBinaryImageType::New();
  noduleLabelToBinaryFilter->SetInput( noduleLabelMap );
  noduleLabelToBinaryFilter->SetBackgroundValue( 0 );
  noduleLabelToBinaryFilter->SetForegroundValue( 1 );

  // Use the ROI without the pleura

  // Mask the input image with the lung label
  clonedInput = inputDuplicator->GetOutput();
  typedef itk::MaskImageFilter< InputImageType, OutputImageType, InputImageType > MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  maskFilter->SetInput(clonedInput);
  maskFilter->SetMaskImage(lungLabelImage);
  maskFilter->SetOutsideValue(-1024);

  ROIFilterType::Pointer maskedImageROIFilter = ROIFilterType::New();
  maskedImageROIFilter->SetRegionOfInterest( roiRegion );
  maskedImageROIFilter->SetInput( maskFilter->GetOutput() );
  maskedImageROIFilter->Update();

  InputImageType::Pointer roiImage = InputImageType::New();
  roiImage = maskedImageROIFilter->GetOutput();

  // Duplicate roiImage to avoid delays!
  DuplicatorType::Pointer roiDuplicator = DuplicatorType::New();
  roiDuplicator->SetInputImage(roiImage);
  roiDuplicator->Update();
  InputImageType::Pointer clonedROI = roiDuplicator->GetOutput();   

  // Exclude all the voxels with an intensity smaller than -824 HU from the following analysis and segmentation 
  typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > thresholdFilterType;
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
  labelStatisticsImageFilter->SetInput(clonedROI);

  labelStatisticsImageFilter->Update();

  LabelStatisticsImageFilterType::RealType mean = labelStatisticsImageFilter->GetMean( 1 );
  LabelStatisticsImageFilterType::RealType sigma = labelStatisticsImageFilter->GetSigma( 1 );

  LabelStatisticsImageFilterType::RealType upperThreshold = mean + 2*sigma;
  LabelStatisticsImageFilterType::RealType lowerThreshold = mean - 2*sigma;

  std::cout<<"Range: "<<lowerThreshold<<"->"<<upperThreshold<<std::endl;
  // Segment the nodule
  clonedROI = roiDuplicator->GetOutput();

  int noduleColor = cip::NODULAR;
  
  typedef itk::ConnectedThresholdImageFilter< InputImageType, OutputImageType > ConnectedThresholdFilterType;
  ConnectedThresholdFilterType::Pointer connectedThresholdFilter = ConnectedThresholdFilterType::New();
  connectedThresholdFilter->SetInput( clonedROI );			  
  connectedThresholdFilter->SetReplaceValue( noduleColor ); 	
  connectedThresholdFilter->SetUpper( upperThreshold );
  connectedThresholdFilter->SetLower( lowerThreshold );
  //connectedThresholdFilter->SetConnectivity(ConnectedThresholdFilterType::FullConnectivity);

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
  double n_voxels = statisticsFilter->GetSum();

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
			
                connectedThresholdFilter->SetInput( clonedROI );
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

  for( unsigned int i = 0; i < 10; i++)
  {
    // Compute mean and std. dev. of the region underneath the image label
    clonedROI = roiDuplicator->GetOutput();

    labelStatisticsImageFilter->SetLabelInput( noduleLabel );
    labelStatisticsImageFilter->SetInput( clonedROI );
    labelStatisticsImageFilter->Update();
    
    mean = labelStatisticsImageFilter->GetMean( noduleColor );
    sigma = labelStatisticsImageFilter->GetSigma( noduleColor );

    upperThreshold = mean + 2*sigma;
    lowerThreshold = mean - 2*sigma;

    // Segment the nodule
    clonedROI = roiDuplicator->GetOutput();   

    connectedThresholdFilter->SetInput( clonedROI );			  
    connectedThresholdFilter->SetReplaceValue( noduleColor ); 	
    connectedThresholdFilter->SetUpper( upperThreshold );
    connectedThresholdFilter->SetLower( lowerThreshold );    
    connectedThresholdFilter->AddSeed( index );
    connectedThresholdFilter->Update();

    noduleLabel = connectedThresholdFilter->GetOutput();
  }

  typedef itk::BinaryBallStructuringElement< OutputPixelType, Dim > StructuringElementType;  	
  StructuringElementType structElement;
  StructuringElementType::SizeType radius;
  radius.Fill( 1 );

  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();
	
  typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > ClosingFilterType;
  ClosingFilterType::Pointer closing = ClosingFilterType::New();
  closing->SetInput( noduleLabel );
  closing->SetKernel( structElement );
  closing->SetForegroundValue( noduleColor );
  closing->Update();

  I2LType::Pointer n2iFilter = I2LType::New();
  n2iFilter->SetInput( closing->GetOutput() );
  n2iFilter->SetInputForegroundValue( noduleColor );
  n2iFilter->Update();

  ShapeLabelObjectType::Pointer noduleLabelObject = n2iFilter->GetOutput()->GetNthLabelObject(0);
  std::cout<<"Roundness: "<<noduleLabelObject->GetRoundness()<<std::endl;
  std::cout<<"Area mean: "<<mean<<std::endl;
  std::cout<<"sigma: "<<sigma<<std::endl;

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

  return EXIT_SUCCESS;
}
