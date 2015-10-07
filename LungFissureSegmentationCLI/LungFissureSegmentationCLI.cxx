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

#include "itkPasteImageFilter.h"

#include "itkNumericTraits.h"

#include "itkMaskImageFilter.h"
#include "itkMaskNegatedImageFilter.h"

#include "itkGaussianEnhancementImageFilter.h"
#include "itkFissuresEigenValueFunctor.h"

#include "itkSymmetricEigenAnalysis.h"

#include "itkThresholdImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"

#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageSliceIteratorWithIndex.h"

#include "itkMergeLabelMapFilter.h"

#include "itkStatisticsImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "LungFissureSegmentationCLICLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const unsigned int Dim = 3;

  typedef signed int												InputPixelType;
  typedef unsigned int												OutputPixelType;
  typedef unsigned int 												LabelType;
  typedef float														GEOutputPixelType;
  typedef unsigned char												VoronoiPixelType;

  typedef itk::Image< InputPixelType,  Dim > 						InputImageType;
  typedef itk::Image< OutputPixelType, Dim > 						OutputImageType;
  typedef itk::Image< LabelType,       Dim > 						LabelImageType;
  typedef itk::Image< GEOutputPixelType, Dim > 						GEOutputImageType;
  typedef itk::Image< VoronoiPixelType, Dim > 						VoronoiImageType;

  typedef itk::ImageFileReader<InputImageType>						InputReaderType;
  typedef itk::ImageFileReader<OutputImageType>						LungLabelReaderType;

  typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType >	I2LType;
  typedef I2LType::OutputImageType									LabelMapType;
  typedef LabelMapType::LabelObjectType								ShapeLabelObjectType;

  InputReaderType::Pointer inputReader = InputReaderType::New();
  inputReader->SetFileName( InputVolume.c_str() );
  inputReader->Update();
  InputImageType::Pointer inputImage = inputReader->GetOutput();  

  LungLabelReaderType::Pointer lungLabelReader = LungLabelReaderType::New();
  lungLabelReader->SetFileName( InputLungLabel.c_str() );
  lungLabelReader->Update();
  OutputImageType::Pointer lungLabelImage = lungLabelReader->GetOutput();

  // Mask the input image with the lung label
  typedef itk::MaskImageFilter< InputImageType, OutputImageType, InputImageType > MaskFilterType;
  MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  maskFilter->SetInput(inputImage);
  maskFilter->SetMaskImage(lungLabelImage);
  maskFilter->SetOutsideValue(300);
  //maskFilter->Update();

  //InputImageType::Pointer hilumImage = maskFilter->GetOutput();

  typedef itk::ThresholdImageFilter<InputImageType> ThresholdImageFilterType;
  ThresholdImageFilterType::Pointer thresholdImageFilter = ThresholdImageFilterType::New();

  thresholdImageFilter->SetInput(maskFilter->GetOutput());
  thresholdImageFilter->ThresholdAbove(0);
  thresholdImageFilter->SetOutsideValue(0);
  thresholdImageFilter->Update();

  InputImageType::Pointer parenchymaImage = thresholdImageFilter->GetOutput();

  /** Segment Vessels */
  /*typedef itk::ImageRegionIteratorWithIndex<InputImageType> HilumRegionIteratorType;
  HilumRegionIteratorType hilumIt(hilumImage, hilumImage->GetRequestedRegion());

  typedef itk::ImageRegionIterator<OutputImageType> LungLabelIteratorType;
  LungLabelIteratorType lungLabelIt(lungLabelImage, lungLabelImage->GetRequestedRegion());

  hilumIt.GoToReverseBegin();
  lungLabelIt.GoToEnd();

  InputImageType::IndexType firstFiducial;
  InputImageType::IndexType secondFiducial;

  bool firstFiducialFound = 0;
  bool secondFiducialFound = 0;

  while( !hilumIt.IsAtReverseEnd() )
  {
	  if( hilumIt.Get() > 20 && hilumIt.Get() < 200 )
	  {
		  if( lungLabelIt.Get() == 1 && !firstFiducialFound )
		  {
			  firstFiducial = hilumIt.GetIndex();
			  firstFiducialFound = 1;
		  }
		  if( lungLabelIt.Get() == 2 && !secondFiducialFound )
		  {
			  secondFiducial = hilumIt.GetIndex();
			  secondFiducialFound = 1;
		  }
	  }
	  
	  if( firstFiducialFound && secondFiducialFound )
	  {
		  hilumIt.GoToBegin();
	  }
	  --hilumIt;
	  --lungLabelIt;
  }

  typedef itk::ConnectedThresholdImageFilter<InputImageType,OutputImageType> ConnectedThresholdImageFilterType;
  ConnectedThresholdImageFilterType::Pointer connectedThresholdFilter = ConnectedThresholdImageFilterType::New();

  connectedThresholdFilter->SetInput(hilumImage);
  connectedThresholdFilter->AddSeed(firstFiducial);
  connectedThresholdFilter->AddSeed(secondFiducial);
  connectedThresholdFilter->SetReplaceValue(1);
  connectedThresholdFilter->SetUpper(200);
  connectedThresholdFilter->SetLower(-500);
  connectedThresholdFilter->Update();

  OutputImageType::Pointer vesselsLabel = connectedThresholdFilter->GetOutput();

  typedef itk::DanielssonDistanceMapImageFilter< OutputImageType, InputImageType, VoronoiImageType > DanielssonDMFilterType;
  DanielssonDMFilterType::Pointer distanceMapFilter = DanielssonDMFilterType::New();

  distanceMapFilter->SetInput( vesselsLabel );
  distanceMapFilter->SetInputIsBinary(1);
  distanceMapFilter->SetSquaredDistance(1);
  //distanceMapFilter->Update();

  typedef itk::MaskImageFilter<InputImageType, OutputImageType> DistanceMapMaskFilterType;
  DistanceMapMaskFilterType::Pointer DMMaskFilter = DistanceMapMaskFilterType::New();
  DMMaskFilter->SetInput( distanceMapFilter->GetOutput() );
  DMMaskFilter->SetMaskImage( lungLabelImage );
  DMMaskFilter->SetOutsideValue( 0.0 );

  ThresholdImageFilterType::Pointer vesselsDMThresholdFilter = ThresholdImageFilterType::New();
  vesselsDMThresholdFilter->SetInput( DMMaskFilter->GetOutput() );
  vesselsDMThresholdFilter->ThresholdBelow( 50 );
  vesselsDMThresholdFilter->SetOutsideValue( 0 );
  vesselsDMThresholdFilter->Update();

  InputImageType::Pointer vesselsDMImage = vesselsDMThresholdFilter->GetOutput();*/

  /** Enhance Plate-Like Structures And Segment Fissures */
  typedef itk::GaussianEnhancementImageFilter<InputImageType,GEOutputImageType> MSGaussianEnhancementFilterType;
  MSGaussianEnhancementFilterType::Pointer GEFilter = MSGaussianEnhancementFilterType::New();

  GEFilter->SetInput( parenchymaImage );
  GEFilter->SetNormalizeAcrossScale(0);
  GEFilter->SetSigma(1.0);
  GEFilter->SetRescale( 0 );

  typedef itk::Functor::FissuresEigenValueFunctor< MSGaussianEnhancementFilterType::EigenValueArrayType, GEOutputPixelType > FunctorType;
  FunctorType::Pointer functor = FunctorType::New();

  GEFilter->SetUnaryFunctor( functor );
  GEFilter->Update();
  GEOutputImageType::Pointer GEImage = GEFilter->GetOutput();

  typedef itk::ImageRegionIterator<GEOutputImageType> GERegionIteratorType;
  GERegionIteratorType GEIt(GEImage, GEImage->GetRequestedRegion());

  typedef itk::ImageRegionIterator<InputImageType> inputRegionIteratorType;
  inputRegionIteratorType parenchymaIt(parenchymaImage, parenchymaImage->GetRequestedRegion());

  GEIt.GoToBegin();
  parenchymaIt.GoToBegin();

  while( !GEIt.IsAtEnd() )
  {
    if( GEIt.Get() <= 0.1 || parenchymaIt.Get() >= -300 || parenchymaIt.Get() <= -900)
    {
      GEIt.Set(itk::NumericTraits<GEOutputPixelType>::Zero);
    }
    ++parenchymaIt;
    ++GEIt;
  }

  typedef itk::Image< itk::SymmetricSecondRankTensor< float, Dim >, Dim > HessianTensorImageType;  	
  typedef HessianTensorImageType::ConstPointer HessianTensorConstPointer;

  HessianTensorConstPointer hessianImage = GEFilter->GetHessianImage();

  typedef itk::Matrix< float, Dim, Dim > EigenVectorMatrixType;

  typedef itk::SymmetricEigenAnalysis<MSGaussianEnhancementFilterType::HessianTensorImageType::PixelType,
     MSGaussianEnhancementFilterType::EigenValueArrayType, EigenVectorMatrixType> SymmetricEigenAnalysisType;

  SymmetricEigenAnalysisType symmetricEigenSystem(Dim);
  symmetricEigenSystem.SetOrderEigenMagnitudes( true );

  MSGaussianEnhancementFilterType::EigenValueArrayType eigenValues;
  EigenVectorMatrixType eigenVectors;

  EigenVectorMatrixType inMatrix;
  inMatrix.Fill(0.0);

  typedef MSGaussianEnhancementFilterType::EigenValueArrayType PrincipalEigenVectorType; 

  typedef itk::Image<PrincipalEigenVectorType, Dim> PrincipalEigenVectorImageType;
  PrincipalEigenVectorImageType::Pointer principalEVImage = PrincipalEigenVectorImageType::New();

  principalEVImage->SetRegions( hessianImage->GetRequestedRegion() );                                   
  principalEVImage->SetBufferedRegion( hessianImage->GetBufferedRegion() );
  principalEVImage->SetLargestPossibleRegion( hessianImage->GetLargestPossibleRegion() );
  principalEVImage->CopyInformation( hessianImage );
  principalEVImage->Allocate();
  principalEVImage->FillBuffer(itk::NumericTraits<PrincipalEigenVectorType>::Zero);

  typedef itk::ImageRegionConstIterator<MSGaussianEnhancementFilterType::HessianTensorImageType> RegionConstIteratorIndexType;
  RegionConstIteratorIndexType hessianIt(hessianImage, hessianImage->GetRequestedRegion());

  typedef itk::ImageRegionIterator<PrincipalEigenVectorImageType> RegionIteratorType;
  RegionIteratorType EVIt(principalEVImage, principalEVImage->GetRequestedRegion());

  hessianIt.GoToBegin();
  EVIt.GoToBegin();
  GEIt.GoToBegin();

  while( !GEIt.IsAtEnd() )
  {
    if( GEIt.Get() >= 0.1 )
    {
      symmetricEigenSystem.ComputeEigenValuesAndVectors(hessianIt.Get(), eigenValues, eigenVectors);
      EVIt.Set(eigenVectors[2]);
    }
    ++GEIt;
    ++hessianIt;
    ++EVIt;
  }

  OutputImageType::SizeType   neighborhoodRadius;

  neighborhoodRadius.Fill(3); 
  //neighborhoodRadius[2]=0;

  typedef itk::ConstNeighborhoodIterator<GEOutputImageType> neighborhoodIterator;
  neighborhoodIterator neighborIterator(neighborhoodRadius, GEImage, GEImage->GetRequestedRegion());

  OutputImageType::Pointer fissuresImage = OutputImageType::New();
  fissuresImage->SetRegions( GEImage->GetRequestedRegion() );                                   
  fissuresImage->SetBufferedRegion( GEImage->GetBufferedRegion() );
  fissuresImage->SetLargestPossibleRegion( GEImage->GetLargestPossibleRegion() );
  fissuresImage->CopyInformation( GEImage );
  fissuresImage->Allocate();
  fissuresImage->FillBuffer(itk::NumericTraits<OutputPixelType>::Zero);

  typedef itk::ImageRegionIterator<OutputImageType> OutputRegionIteratorType;
  OutputRegionIteratorType outputIt(fissuresImage, fissuresImage->GetRequestedRegion());

  typedef itk::Vector<float, Dim> VectorType;
  VectorType centralEV;
  VectorType adjacentEV;
  VectorType offset;
  VectorType outerProduct;

  neighborIterator.GoToBegin();
  outputIt.GoToBegin();
  while(!neighborIterator.IsAtEnd())
  {
	  if( neighborIterator.GetCenterPixel() != 0.0 )
	  {
		  for(unsigned int i = 0; i < neighborIterator.Size(); ++i)
		  {
			  if(neighborIterator.GetPixel(i) > 0.1)
			  {
				  centralEV[0] = principalEVImage->GetPixel(neighborIterator.GetIndex())[0];
				  centralEV[1] = principalEVImage->GetPixel(neighborIterator.GetIndex())[1];
				  centralEV[2] = principalEVImage->GetPixel(neighborIterator.GetIndex())[2];
				  centralEV.Normalize();
				  adjacentEV[0] = principalEVImage->GetPixel(neighborIterator.GetIndex()-neighborIterator.GetOffset(i))[0];
				  adjacentEV[1] = principalEVImage->GetPixel(neighborIterator.GetIndex()-neighborIterator.GetOffset(i))[1];
				  adjacentEV[2] = principalEVImage->GetPixel(neighborIterator.GetIndex()-neighborIterator.GetOffset(i))[2];
				  adjacentEV.Normalize();
				  float innerProduct = centralEV * adjacentEV;
				  
				  offset[0] = neighborIterator.GetOffset(i)[0];
				  offset[1] = neighborIterator.GetOffset(i)[1];
				  offset[2] = neighborIterator.GetOffset(i)[2];   
				  offset.Normalize();
				  
				  outerProduct[0] = (offset[1]*centralEV[2]) - (offset[2]*centralEV[1]);
				  outerProduct[1] = (offset[2]*centralEV[0]) - (offset[0]*centralEV[2]);
				  outerProduct[2] = (offset[0]*centralEV[1]) - (offset[1]*centralEV[0]);
				  
				  float finalOuterProduct = sqrt(outerProduct[0]*outerProduct[0] + outerProduct[1]*outerProduct[1] + outerProduct[2]*outerProduct[2]);
				  if(innerProduct >= 0.985 && finalOuterProduct >= 0.985)
				  {
					  fissuresImage->SetPixel(neighborIterator.GetIndex()-neighborIterator.GetOffset(i), itk::NumericTraits<OutputPixelType>::One);
				  }
			  }
		  }
	  }
	  ++neighborIterator;
	  ++outputIt;
  }
  
  typedef itk::ConnectedComponentImageFilter <OutputImageType, OutputImageType> ConnectedComponentImageFilterType;

  ConnectedComponentImageFilterType::Pointer connectedComponentFilter = ConnectedComponentImageFilterType::New ();
  connectedComponentFilter->SetInput(fissuresImage);

  typedef itk::RelabelComponentImageFilter<OutputImageType, OutputImageType> RelabelImageFilterType;
  RelabelImageFilterType::Pointer relabelFilter = RelabelImageFilterType::New();

  relabelFilter->SetInput(connectedComponentFilter->GetOutput());
  relabelFilter->SetMinimumObjectSize(8000); 
  relabelFilter->Update();
  fissuresImage = relabelFilter->GetOutput();

  typedef itk::RegionOfInterestImageFilter<OutputImageType, OutputImageType> ROIFilterType;
  ROIFilterType::Pointer roiFilter = ROIFilterType::New();

  OutputImageType::SizeType regionSize;
  OutputImageType::IndexType regionIndex;

  OutputImageType::RegionType ROI;

  regionSize[0] = 1;
  regionSize[1] = fissuresImage->GetLargestPossibleRegion().GetSize()[1];
  regionSize[2] = fissuresImage->GetLargestPossibleRegion().GetSize()[2];

  regionIndex.Fill(0);

  ROI.SetSize(regionSize);
  //roiFilter->Update();

  OutputImageType::Pointer cleanedFissuresImage = OutputImageType::New();
  cleanedFissuresImage->SetRegions( fissuresImage->GetRequestedRegion() );                                   
  cleanedFissuresImage->SetBufferedRegion( fissuresImage->GetBufferedRegion() );
  cleanedFissuresImage->SetLargestPossibleRegion( fissuresImage->GetLargestPossibleRegion() );
  cleanedFissuresImage->CopyInformation( fissuresImage );
  cleanedFissuresImage->Allocate();
  cleanedFissuresImage->FillBuffer(itk::NumericTraits<OutputPixelType>::Zero);

  typedef itk::BinaryImageToShapeLabelMapFilter<OutputImageType> ImageToShapeLabelMapFilterType;
  typedef I2LType::OutputImageType LabelMapType;
  typedef LabelMapType::LabelObjectType ShapeLabelObjectType;

  typedef itk::MergeLabelMapFilter<LabelMapType> MergeLabelMapFilterType;
  MergeLabelMapFilterType::Pointer mergeLabelFilter = MergeLabelMapFilterType::New();

  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, OutputImageType > LabelMapToBinaryImageType;
  LabelMapToBinaryImageType::Pointer L2BFilter = LabelMapToBinaryImageType::New();

  typedef itk::PasteImageFilter<OutputImageType, OutputImageType> PasteImageFilterType;
  PasteImageFilterType::Pointer pasteImageFilter = PasteImageFilterType::New();

  for( unsigned int sliceN = 0; sliceN < fissuresImage->GetLargestPossibleRegion().GetSize()[0]; sliceN++ ) //
  {
	  regionIndex[0] = sliceN;
	  ROI.SetIndex(regionIndex);

	  roiFilter->SetRegionOfInterest(ROI);
	  roiFilter->SetInput(fissuresImage);
	  roiFilter->Update();
	  OutputImageType::Pointer sagittalSlice = roiFilter->GetOutput();

	  LabelMapType::Pointer totalLabelMap = LabelMapType::New();

	  for( unsigned int label = 1; label < 3; label++ )
	  {
		  ImageToShapeLabelMapFilterType::Pointer I2SLFilter = ImageToShapeLabelMapFilterType::New();
		  I2SLFilter->SetInputForegroundValue(label);
		  I2SLFilter->SetInput(roiFilter->GetOutput());
		  I2SLFilter->SetFullyConnected(1);
		  I2SLFilter->Update();
		  
		  LabelMapType::Pointer labelMap = I2SLFilter->GetOutput();
		  
		  std::vector<unsigned long> labelsToRemove;
		  
		  for(unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); i++)
		  {
			  ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);
			  if( labelObject->GetRoundness() > 0.5 || labelObject->GetNumberOfPixels() < 50 )
			  {
				  labelsToRemove.push_back(labelObject->GetLabel());
			  }
			  else
			  {
				  std::cout<<"centroid: "<<labelObject->GetCentroid()<<" round: "<<labelObject->GetRoundness()<<" pixels: "<<labelObject->GetNumberOfPixels()<<std::endl;
			  }
		  }
		  for(unsigned int n = 0; n < labelsToRemove.size(); n++)
		  {		  
			  labelMap->RemoveLabel(labelsToRemove[n]);
		  }

		  if( label == 2 )
		  {
			  mergeLabelFilter->SetMethod(MergeLabelMapFilterType::PACK);
			  mergeLabelFilter->SetInput(totalLabelMap);
			  mergeLabelFilter->SetInput(1,labelMap);
			  mergeLabelFilter->Update();
			  totalLabelMap = mergeLabelFilter->GetOutput();
		  }
		  else
		  {
			  totalLabelMap = labelMap;
		  }
	  }
	  
	  L2BFilter->SetInput(totalLabelMap);
	  L2BFilter->SetBackgroundValue( 0 );
	  L2BFilter->SetForegroundValue( 1 );
	  L2BFilter->Update();

	  pasteImageFilter->SetSourceImage(L2BFilter->GetOutput());
	  pasteImageFilter->SetDestinationImage(cleanedFissuresImage);
	  pasteImageFilter->SetSourceRegion(L2BFilter->GetOutput()->GetLargestPossibleRegion());
	  pasteImageFilter->SetDestinationIndex(regionIndex);
	  pasteImageFilter->Update();

	  cleanedFissuresImage = pasteImageFilter->GetOutput();
  }

  /*typedef itk::BinaryBallStructuringElement< OutputPixelType, Dim > StructuringElementType;  	
  StructuringElementType structElement;
  StructuringElementType::SizeType radius;
  radius.Fill( 5 );

  structElement.SetRadius( radius );
  structElement.CreateStructuringElement();

  typedef itk::BinaryMorphologicalClosingImageFilter < OutputImageType, OutputImageType, StructuringElementType > ClosingFilterType;
  ClosingFilterType::Pointer closing = ClosingFilterType::New();
  closing->SetInput( cleanedFissuresImage );
  closing->SetKernel( structElement );
  closing->SetForegroundValue( 1 );*/

  /*typedef itk::DanielssonDistanceMapImageFilter< OutputImageType, InputImageType, VoronoiImageType > DanielssonDMFilterType;
  DanielssonDMFilterType::Pointer distanceMapFilter = DanielssonDMFilterType::New();*/

  /*distanceMapFilter->SetInput( cleanedFissuresImage );
  distanceMapFilter->SetInputIsBinary(1);
  distanceMapFilter->Update();
  
  InputImageType::Pointer fissuresDMImage = distanceMapFilter->GetOutput();*/

  /*typedef itk::SubtractImageFilter< InputImageType > SubtractImageType;
  SubtractImageType::Pointer subtractFilter = SubtractImageType::New();
  subtractFilter->SetConstant1(10);
  subtractFilter->SetInput2( distanceMapFilter->GetOutput() );
  subtractFilter->Update();*/
  /*typedef itk::ImageDuplicator<OutputImageType> ImageDuplicatorType;
  ImageDuplicatorType::Pointer imageDuplicator = ImageDuplicatorType::New();
  imageDuplicator->SetInputImage(lungLabelImage);
  imageDuplicator->Update();

  OutputImageType::Pointer fissureLungLabelImage = imageDuplicator->GetOutput();*/
  
  /*inputRegionIteratorType fissureDMIt( fissuresDMImage, fissuresDMImage->GetRequestedRegion() );
  OutputRegionIteratorType lungLabelIt( fissureLungLabelImage, fissureLungLabelImage->GetRequestedRegion() );

  fissureDMIt.GoToBegin();
  lungLabelIt.GoToBegin();

  while( !fissureDMIt.IsAtEnd() )
  {
	  if( fissureDMIt.Get() <= 3 )
	  {
		  OutputPixelType value = lungLabelIt.Get();
		  value += 2;
		  lungLabelIt.Set(value);

    }
    ++fissureDMIt;
	++lungLabelIt;
  }

  typedef itk::ImageSliceIteratorWithIndex<OutputImageType> SliceIterator;
  SliceIterator sIt(fissureLungLabelImage, fissureLungLabelImage->GetLargestPossibleRegion());

  sIt.SetFirstDirection(0);
  sIt.SetSecondDirection(1);

  sIt.GoToBegin();

  while( !sIt.IsAtEnd() )
  {
	  while( !sIt.IsAtEndOfSlice() )
	  {
		  int count = 0;
		  OutputPixelType previous;
		  while( !sIt.IsAtEndOfLine() )
		  {
			  if( sIt.Get() == 3 )
			  {
				  count++;
			  }
			  if( sIt.Get() == 1 && previous == 3 && count >= 50 )
			  {
				  sIt.Set( 3 );
			  }
			  previous = sIt.Get();			  
			  ++sIt;
		  }
		  sIt.NextLine();
	  }
	  sIt.NextSlice();
  }*/

  //fissuresPasteImageFilter->Update();
 

  
  /*typedef itk::AddImageFilter< InputImageType, InputImageType > AddImageFilterType;
  AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

  addFilter->SetInput1( fissureDMImage );
  addFilter->SetInput2( parenchymaImage );
  //addFilter->Update();*/

  /*typedef itk::MaskNegatedImageFilter< InputImageType, OutputImageType > MaskNegatedFilterType;
  MaskNegatedFilterType::Pointer maskNegatedFilter = MaskNegatedFilterType::New();

  maskNegatedFilter->SetInput(parenchymaImage);
  maskNegatedFilter->SetMaskImage(vesselsLabel);
  maskNegatedFilter->SetOutsideValue(0);
  maskNegatedFilter->Update();
  
  InputImageType::Pointer maskNegatedImage = maskNegatedFilter->GetOutput();

  //InputImageType::Pointer fissuresOriginalImage = addFilter->GetOutput();*/
  
  /*OutputRegionIteratorType fissuresImageIt( cleanedFissuresImage, cleanedFissuresImage->GetRequestedRegion());
  inputRegionIteratorType vesselsDMIt( vesselsDMImage, vesselsDMImage->GetRequestedRegion() );
  inputRegionIteratorType fissureDMIt( fissuresDMImage, fissuresDMImage->GetRequestedRegion() );

  fissuresImageIt.GoToBegin();
  vesselsDMIt.GoToBegin();
  fissureDMIt.GoToBegin();

  while(!fissuresImageIt.IsAtEnd())
  {
	  if( fissuresImageIt.Get() == 0 && vesselsDMIt.Get() > 0 && fissureDMIt.Get() <= 5 )
	  {
		  fissuresImageIt.Set(1);
	  }
  	  ++fissuresImageIt;
	  ++vesselsDMIt;
	  ++fissureDMIt;
  }*/

  /*InputImageType::Pointer vesselsDMImage = distanceMapFilter->GetOutput();

  typedef itk::ImageRegionIterator<InputImageType> inputRegionIteratorType;  
  inputRegionIteratorType vesselsDMIt(vesselsDMImage, vesselsDMImage->GetRequestedRegion());

  vesselsDMIt.GoToBegin();

  while( !vesselsDMIt.IsAtEnd() )
  {
    if( vesselsDMIt.Get() >= 6 )
    {
		InputImageType::PixelType value = vesselsDMIt.Get();
		value *= value;
		vesselsDMIt.Set(value);
    }
	else
	{
		vesselsDMIt.Set(itk::NumericTraits<InputPixelType>::Zero);
	}
    ++vesselsDMIt;
  }*/

  /*typedef itk::MultiplyImageFilter<InputImageType,InputImageType> MultiplyFilterType;
  MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
  multiplyFilter->SetInput1(distanceMapFilter->GetOutput());
  multiplyFilter->SetConstant2(5);
  //distanceMapFilter->Update();*/

  /*typedef itk::AddImageFilter< InputImageType, InputImageType > AddImageFilterType;
  AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

  addFilter->SetInput1( maskNegatedImage );
  addFilter->SetInput2( distanceMapFilter->GetOutput() );*/
  //addFilter->Update();

  /*typedef itk::MaskImageFilter<InputImageType, OutputImageType> DistanceMapMaskFilterType;
  DistanceMapMaskFilterType::Pointer DMMaskFilter = DistanceMapMaskFilterType::New();
  DMMaskFilter->SetInput( distanceMapFilter->GetOutput() );
  DMMaskFilter->SetMaskImage( lungLabelImage );
  DMMaskFilter->SetOutsideValue( 0.0 );*/
  
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( OutputVolume.c_str() );
  writer->SetInput( cleanedFissuresImage );
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
