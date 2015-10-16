#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "cipLabelMapToLungLobeLabelMapImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "itkImageSliceIteratorWithIndex.h"

#include "LungLobeSegmentationCLICLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  cip::ChestConventions conventions;
  
  const unsigned int Dim = 3;

  typedef signed int						InputPixelType;
  typedef unsigned short					LabelPixelType;

  typedef itk::Image< InputPixelType, Dim > InputImageType;
  typedef itk::Image< LabelPixelType, Dim > LabelImageType;
  
  typedef itk::ImageFileReader< InputImageType > GreyScaleReaderType;
  GreyScaleReaderType::Pointer greyScaleReader = GreyScaleReaderType::New();     
  greyScaleReader->SetFileName( InputVolume.c_str() );
  greyScaleReader->Update();

  InputImageType::Pointer inputImage = greyScaleReader->GetOutput();

  typedef itk::ImageFileReader< LabelImageType > LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( LabelMapVolume.c_str() );
  labelReader->Update();

  LabelImageType::Pointer lungLabelMap = labelReader->GetOutput();

  std::vector< LabelImageType::IndexType > leftOIdxVec;
  std::vector< LabelImageType::IndexType > rightOIdxVec;
  std::vector< LabelImageType::IndexType > rightHIdxVec;

  LabelReaderType::Pointer fissuresReader = LabelReaderType::New();
  fissuresReader->SetFileName(FissuresVolume.c_str() );
  fissuresReader->Update();

  LabelImageType::Pointer fissuresLabelMap = fissuresReader->GetOutput();

  typedef itk::ImageSliceIteratorWithIndex<LabelImageType> LabelMapSliceIteratorType;
  LabelMapSliceIteratorType fIt( fissuresLabelMap, fissuresLabelMap->GetBufferedRegion() );
  LabelMapSliceIteratorType mIt( lungLabelMap, lungLabelMap->GetBufferedRegion() );

  fIt.SetFirstDirection(0);
  fIt.SetSecondDirection(2);

  mIt.SetFirstDirection(0);
  mIt.SetSecondDirection(2);

  fIt.GoToBegin();
  mIt.GoToBegin();

  LabelImageType::IndexType rhIdx;
  LabelImageType::IndexType roIdx;
  LabelImageType::IndexType idx;

  while( !fIt.IsAtEnd() )
  {
	  while( !fIt.IsAtEndOfSlice() )
	  {
		  while( !fIt.IsAtEndOfLine() )
		  {
			  idx = fIt.GetIndex();
			  if( fIt.Get() != 0 && mIt.Get() >= cip::LEFTUPPERTHIRD && mIt.Get() <= cip::LEFTLOWERTHIRD )
			  {
				  leftOIdxVec.push_back( idx );
			  }
			  else if( fIt.Get() != 0 && mIt.Get() >= cip::RIGHTUPPERTHIRD && mIt.Get() <= cip::RIGHTLOWERTHIRD )
			  {
				  if( mIt.Get() == cip::RIGHTLOWERTHIRD )
				  {
					  rightOIdxVec.push_back( idx );
				  }
				  else
				  {
					  if( rightHIdxVec.size() > 0 && rightOIdxVec.size() > 0 )
					  {
						  rhIdx = rightHIdxVec.back();
						  roIdx = rightOIdxVec.back();
						  
						  if( abs(static_cast<int>(rhIdx[2]-idx[2])) <= abs(static_cast<int>(roIdx[2]-idx[2])) )
						  {
							  rightHIdxVec.push_back(idx);
						  }
						  else
						  {
								  rightOIdxVec.push_back(idx);
								  fIt.Set( 5 );
						  }
					  }
					  else if( rightHIdxVec.size() == 0 && rightOIdxVec.size() > 0 )
					  {
						  roIdx = rightOIdxVec.back();
						  if( abs(static_cast<int>(roIdx[2]-idx[2])) <= 5 )
						  {
							  rightOIdxVec.push_back(idx);
						  }
						  else
						  {
							  rightHIdxVec.push_back(idx);
						  }
					  }
					  else if( rightOIdxVec.size() == 0 && rightHIdxVec.size() > 0 )
					  {
						  rhIdx = rightOIdxVec.back();
						  if( abs(static_cast<int>(rhIdx[2]-idx[2])) <= 5 )
						  {
							  rightHIdxVec.push_back(idx);
						  }
						  else
						  {
							  rightOIdxVec.push_back(idx);
						  }
					  }
					  else //TODO if both empty choose based on position within the lung third
					  {
					  }
				  }
			  }
			  ++fIt;
			  ++mIt;
		  }
		  fIt.NextLine();
		  mIt.NextLine();		  
	  }
	  fIt.NextSlice();
	  mIt.NextSlice();
  }

  std::vector< LabelImageType::IndexType > rightObliqueIndicesVec;
  std::vector< LabelImageType::IndexType > rightHorizontalIndicesVec;
  std::vector< LabelImageType::IndexType > leftObliqueIndicesVec;

  unsigned int leftStep, rightOStep, rightHStep;
  if( leftOIdxVec.size() > 100 )
  {
	  leftStep = int(leftOIdxVec.size()/100);
  }
  else
  {
	  leftStep = 1;
  }
  if( rightOIdxVec.size() > 200 )
  {
	  rightOStep = int(rightOIdxVec.size()/200);
  }
  else
  {
	  rightOStep = 1;
  }
  if( rightHIdxVec.size() > 200 )
  {
	  rightHStep = int(rightHIdxVec.size()/200);
  }
  else
  {
	  rightHStep = 1;
  }

  std::cout<<rightOIdxVec.size()<<" "<<rightHIdxVec.size()<<std::endl;

  for( unsigned int i = 0; i < leftOIdxVec.size(); i += leftStep )
  {
	  leftObliqueIndicesVec.push_back( leftOIdxVec.at(i) );
  }

  for( unsigned int i = 0; i < rightOIdxVec.size(); i += rightOStep )
  {
	  rightObliqueIndicesVec.push_back( rightOIdxVec.at(i) );
  }

  for( unsigned int j = 0; j < rightHIdxVec.size(); j += rightHStep )
  {
	  rightHorizontalIndicesVec.push_back( rightHIdxVec.at(j) );
  }
  
  /*LabelImageType::PointType point;
  LabelImageType::IndexType idx;

  for( unsigned int i = 0; i < leftObliqueSeeds.size(); i++ )
  { 		
        // Convert to lps the seed point
        point[0] = leftObliqueSeeds[i][0] * (-inputImage->GetDirection()[0][0]);	
        point[1] = leftObliqueSeeds[i][1] * (-inputImage->GetDirection()[1][1]);
        point[2] = leftObliqueSeeds[i][2] *   inputImage->GetDirection()[2][2];

        // Convert the lps physical point to index
        inputImage->TransformPhysicalPointToIndex( point, idx );
		leftOIdxVec.push_back(idx);
  }
  
  for( unsigned int i = 0; i < rightObliqueSeeds.size(); i++ )
  { 	
	  // Convert to lps the seed point
	  point[0] = rightObliqueSeeds[i][0] * (-inputImage->GetDirection()[0][0]);	
	  point[1] = rightObliqueSeeds[i][1] * (-inputImage->GetDirection()[1][1]);
	  point[2] = rightObliqueSeeds[i][2] *   inputImage->GetDirection()[2][2];

	  // Convert the lps physical point to index
	  inputImage->TransformPhysicalPointToIndex( point, idx );
	  rightObliqueIndicesVec.push_back(idx);
  }

  for( unsigned int i = 0; i < rightHorizontalSeeds.size(); i++ )
  { 
	  // Convert to lps the seed point
	  point[0] = rightHorizontalSeeds[i][0] * (-inputImage->GetDirection()[0][0]);	
	  point[1] = rightHorizontalSeeds[i][1] * (-inputImage->GetDirection()[1][1]);
	  point[2] = rightHorizontalSeeds[i][2] *   inputImage->GetDirection()[2][2];

	  // Convert the lps physical point to index
	  inputImage->TransformPhysicalPointToIndex( point, idx );
	  rightHorizontalIndicesVec.push_back(idx);
  }*/
  
  typedef cipLabelMapToLungLobeLabelMapImageFilter LobeSegmentationType;
  LobeSegmentationType::Pointer lobeSegmenter = LobeSegmentationType::New();
  lobeSegmenter->SetInput(lungLabelMap);

  if( rightHorizontalIndicesVec.size() > 0 && rightObliqueIndicesVec.size() > 0 )
  {
	  lobeSegmenter->SetRightHorizontalFissureIndices( rightHorizontalIndicesVec );
	  lobeSegmenter->SetRightObliqueFissureIndices( rightObliqueIndicesVec );
  }
  if( leftOIdxVec.size() > 0 )
  {
	  lobeSegmenter->SetLeftObliqueFissureIndices( leftObliqueIndicesVec );
  }

  typedef itk::ImageFileWriter<LobeSegmentationType::OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( OutputVolume.c_str() );
  writer->SetInput( lobeSegmenter->GetOutput() );
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
