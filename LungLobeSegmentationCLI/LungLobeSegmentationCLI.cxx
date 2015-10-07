#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "cipLabelMapToLungLobeLabelMapImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"

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
  greyScaleReader->SetFileName( inputVolume.c_str() );
  greyScaleReader->Update();

  InputImageType::Pointer inputImage = greyScaleReader->GetOutput();

  typedef itk::ImageFileReader< LabelImageType > LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( labelMapVolume.c_str() );
  labelReader->Update();

  LabelImageType::Pointer lungLabelMap = labelReader->GetOutput();

  std::vector< LabelImageType::IndexType > leftObliqueIndicesVec;
  std::vector< LabelImageType::IndexType > rightObliqueIndicesVec;
  std::vector< LabelImageType::IndexType > rightHorizontalIndicesVec;

  LabelReaderType::Pointer fissuresReader = LabelReaderType::New();
  fissuresReader->SetFileName(fissuresVolume.c_str() );
  fissuresReader->Update();

  LabelImageType::Pointer fissuresLabelMap = fissuresReader->GetOutput();

  typedef itk::ImageRegionIteratorWithIndex<LabelImageType> LabelMapIteratorType;
  LabelMapIteratorType fIt( fissuresLabelMap, fissuresLabelMap->GetBufferedRegion() );
  LabelMapIteratorType mIt( lungLabelMap, lungLabelMap->GetBufferedRegion() );

  fIt.GoToBegin();
  mIt.GoToBegin();

  unsigned int leftCount = 0;

  LabelImageType::IndexType rhIdx;
  LabelImageType::IndexType roIdx;
  LabelImageType::IndexType idx;

  while( !fIt.IsAtEnd() )
  {
	  if( leftCount == 500 && fIt.Get() != 0 && mIt.Get() >= 9 && mIt.Get() <= 11 )
	  {
		  leftObliqueIndicesVec.push_back( fIt.GetIndex() );
		  leftCount = 0;
	  }
	  else if(fIt.Get() != 0 && mIt.Get() >= 9 && mIt.Get() <= 11)
	  {
		  leftCount++;
	  }
	  else if( fIt.Get() != 0 && mIt.Get() >= 12 && mIt.Get() <= 14 )
	  {
		  if( mIt.Get() == 12 )
		  {
			  rightHorizontalIndicesVec.push_back( fIt.GetIndex() );
		  }
		  else if( mIt.Get() == 14 )
		  {
			  rightObliqueIndicesVec.push_back( fIt.GetIndex() );
		  }
		  /*else
		  {
			  idx = fIt.GetIndex();
			  if( rightHorizontalIndicesVec.size() > 0 && rightObliqueIndicesVec.size() > 0 )
			  {
				  rhIdx = rightHorizontalIndicesVec.back();
				  roIdx = rightObliqueIndicesVec.back();

				  if( abs(static_cast<int>(rhIdx[0]-idx[0])) <= abs(static_cast<int>(roIdx[0]-idx[0])) )
				  {
					  rightHorizontalIndicesVec.push_back(idx);
				  }
				  else
				  {
					  rightObliqueIndicesVec.push_back(idx);
				  }
			  }
			  else if( rightHorizontalIndicesVec.size() == 0 && rightObliqueIndicesVec.size() > 0 )
			  {
				  roIdx = rightObliqueIndicesVec.back();
				  if( abs(static_cast<int>(roIdx[0]-idx[0])) <= 50 )
				  {
					  rightObliqueIndicesVec.push_back(idx);
				  }
				  else
				  {
					  rightHorizontalIndicesVec.push_back(idx);
				  }
			  }
			  else if( rightObliqueIndicesVec.size() == 0 && rightHorizontalIndicesVec.size() > 0 )
			  {
				  rhIdx = rightObliqueIndicesVec.back();
				  if( abs(static_cast<int>(rhIdx[0]-idx[0])) <= 50 )
				  {
					  rightHorizontalIndicesVec.push_back(idx);
				  }
				  else
				  {
					  rightObliqueIndicesVec.push_back(idx);
				  }
			  }
			  else //TODO if both empty choose based on position within the lung third
			  {
			  }
		  }*/
	  }
	  ++fIt;
	  ++mIt;
  }

  std::cout<<rightHorizontalIndicesVec.front()<<std::endl;
  std::cout<<rightObliqueIndicesVec.front()<<std::endl;

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
		leftObliqueIndicesVec.push_back(idx);
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

  /*if( rightHorizontalIndicesVec.size() > 0 && rightObliqueIndicesVec.size() > 0 )
  {
	  lobeSegmenter->SetRightHorizontalFissureIndices( rightHorizontalIndicesVec );
	  lobeSegmenter->SetRightObliqueFissureIndices( rightObliqueIndicesVec );
  }*/
  if( leftObliqueIndicesVec.size() > 0 )
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
