/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEnhanceImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/27 17:01:10 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEnhanceImageFilter_h
#define __itkEnhanceImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

#include "itkMultiScaleGaussianEnhancementImageFilter.h"
#include "itkGeneralEigenValueFunctor.h"
#include "itkSymmetricSecondRankTensor.h"

namespace itk
{

/** \class itkEnhanceImageFilter
 *
 * DO NOT assume a particular image or pixel type, which is, the input image
 * may be a VectorImage as well as an Image obeject with vectorial pixel type.
 *
 * \sa Image
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT EnhanceImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef EnhanceImageFilter                    			Self;

	/** Convenient typedefs for simplifying declarations. */
	typedef TInputImage                           			InputImageType;
	typedef typename InputImageType::Pointer      			InputImagePointer;
	typedef typename InputImageType::ConstPointer 			InputImageConstPointer;
	typedef TOutputImage                          			OutputImageType;
	typedef typename OutputImageType::Pointer     			OutputImagePointer;

	/** Standard class typedefs. */
	typedef ImageToImageFilter< InputImageType, OutputImageType>	Superclass;
	typedef SmartPointer<Self>                                  	Pointer;
	typedef SmartPointer<const Self>                             	ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods). */
	itkTypeMacro( EnhanceImageFilter, ImageToImageFilter );
  
	/** Image typedef support. */
	typedef typename InputImageType::PixelType           		InputPixelType;
	typedef typename OutputImageType::PixelType          		OutputPixelType;
	typedef typename InputImageType::RegionType          		InputImageRegionType;
	typedef typename InputImageType::SizeType            		InputImageSizeType;
	typedef typename InputImageType::IndexType          		InputImageIndexType;
	typedef typename InputImageType::SpacingType	   		InputSpacingType;
	typedef typename OutputImageType::RegionType        		OutputImageRegionType;

	/** EigenValue image support */
	typedef itk::Image< double,TInputImage::ImageDimension > 	EigenValueImageType;
	typedef typename EigenValueImageType::ConstPointer       	EigenValueImageConstPointer;
	typedef typename EigenValueImageType::PixelType          	EigenValuePixelType;

	/** Scales image support 
	typedef itk::Image< double,TInputImage::ImageDimension > ScalesImageType;
	typedef typename ScalesImageType::ConstPointer      	 ScalesImageConstPointer;
	typedef typename ScalesImageType::PixelType              ScalesPixelType;*/

	/** Mean typedef support 
	typedef itk::Image< double,TInputImage::ImageDimension >  MeanImageType;
	typedef typename MeanImageType::ConstPointer	     	  MeanImageConstPointer;*/

	/** Gradient typedef support 
 	typedef itk::Image< FixedArray< double,
    	  TInputImage::ImageDimension >,
          TInputImage::ImageDimension >                           GradientImageType;       
	typedef typename GradientImageType::ConstPointer          GradientImageConstPointer;*/

	/** Hessian typedef support 
	typedef itk::Image< SymmetricSecondRankTensor< double, 
	   TInputImage::ImageDimension >, 
	   TInputImage::ImageDimension >    		          HessianTensorImageType;  	
	typedef typename HessianTensorImageType::ConstPointer     HessianTensorConstPointer;*/

	/** Features map
	typedef itk::VectorImage< double, TInputImage::ImageDimension >  FeaturesMapType;
	typedef typename FeaturesMapType::Pointer                        FeaturesMapPointer;
	typedef itk::VariableLengthVector<double> VariableVectorType;*/

	/** Multiscale filter */
 	typedef itk::MultiScaleGaussianEnhancementImageFilter<
          InputImageType, 
          EigenValueImageType >						  MultiScaleFilterType;
	typedef typename MultiScaleFilterType::Pointer			  MultiScaleFilterPointer;
  	typedef typename MultiScaleFilterType::EigenValueArrayType 	  EigenValueArrayType;

  	/** Supported functors. */
  	typedef itk::Functor::GeneralEigenValueFunctor<EigenValueArrayType, 
          ScalesPixelType >						  FunctorType;
	typedef typename FunctorType::Pointer				  FunctorPointer;

	/*typedef std::map< double, vnl_matrix<double> >  		  MatrixMapType;
	typedef std::pair< double, vnl_matrix<double> > 		  MatrixPairType;
	typedef std::map< std::vector<double>, vnl_matrix<double> >  	  BijMapType;
	typedef std::pair< std::vector<double>, vnl_matrix<double> > 	  BijPairType;

	typedef std::map< double, double >			  	  BiOrd0MapType;
	typedef std::pair< double, double >			  	  BiOrd0PairType;

	typedef std::map< std::vector<double>, double>		  	  BijOrd0MapType;
	typedef std::pair< std::vector<double>, double > 		  BijOrd0PairType;

	typedef std::map< double, std::map< unsigned int, 
	  std::vector<unsigned int> > >  		                  IndexVectorMapType;
	typedef std::map< std::vector<double>, std::map< unsigned int, 
	  std::vector<unsigned int> > >  		                  IndexBijMapType;

	typedef std::pair< double, std::map< unsigned int, 
	  std::vector<unsigned int> > > 				  IndexVectorPairType;
	typedef std::pair< std::vector<double>, std::map< unsigned int, 
	  std::vector<unsigned int> > >  		                  IndexBijPairType;

	typedef std::map< double, double >  		    		  TraceMapType;
	typedef std::pair< double, double > 				  TracePairType;*/
	
        //here you have to add the sigma values for the multiscale computation
	itkSetMacro( Sigma,      float              );//noise sigma!! 
	itkGetMacro( Sigma,      float              );
	/*itkSetMacro( H,          float              );
	itkGetMacro( H,          float              );
	itkSetMacro( PSTh,       float              );
    	itkGetMacro( PSTh,       float              );
	itkSetMacro( RSearch,    InputImageSizeType );
	itkGetMacro( RSearch,    InputImageSizeType );
	itkSetMacro( RComp,      InputImageSizeType );
	itkGetMacro( RComp,      InputImageSizeType );*/


        itkSetMacro( MinimumLevel, 	 double       );
	itkGetMacro( MinimumLevel,       double       );
	itkSetMacro( MaximumLevel,       double       );
	itkGetMacro( MaximumLevel,       double       );
	itkSetMacro( NumberOfLevelSteps, unsigned int );
	itkGetMacro( NumberOfLevelSteps, unsigned int );
	itkSetMacro( Order, unsigned int );
	itkGetMacro( Order, unsigned int );
	//itkSetMacro( Alpha, unsigned int );
	//itkGetMacro( Alpha, unsigned int );	

	/** Methods to turn on/off flag to generate an image with scale values at
  	*  each pixel for the best vesselness response */
  	itkSetMacro( GenerateStrengthOutput, bool );
  	itkGetConstMacro( GenerateStrengthOutput, bool );
  	itkBooleanMacro( GenerateStrengthOutput );

  	/*itkSetMacro( GenerateScalesOutput, bool );
  	itkGetConstMacro( GenerateScalesOutput, bool );
  	itkBooleanMacro( GenerateScalesOutput );

  	itkSetMacro( UseEstimatedDistanceMean, bool );
  	itkGetConstMacro( UseEstimatedDistanceMean, bool );
  	itkBooleanMacro( UseEstimatedDistanceMean );

  	itkSetMacro( NormalizedByFeatureStrength, bool );
  	itkGetConstMacro( NormalizedByFeatureStrength, bool );
  	itkBooleanMacro( NormalizedByFeatureStrength );

  	itkSetMacro( UseDeltaFeatureStrength, bool );
  	itkGetConstMacro( UseDeltaFeatureStrength, bool );
  	itkBooleanMacro( UseDeltaFeatureStrength );*/

        //void ComputeFeatures(void);
	//void ComputeLSWeightsAndTrace( const InputImageSizeType&, unsigned int );

  	/** Get the image containing the scales at which each pixel gave the best response */
  	//const ScalesImageType * GetScalesOutput( void ) const;

  	/** Get the image containing the scales at which each pixel gave the best response */
  	//const FeaturesMapType * GetFeaturesOutput( void ) const;

  	/** Get the image containing the strength image */
  	const EigenValueImageType * GetFunctorOutput( void ) const;
	


protected:
	EnhanceImageFilter();
	virtual ~EnhanceImageFilter();
	// Threaded filter:
#if ITK_VERSION_MAJOR < 4
    void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, int threadId );
    
#else
    void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId );
    
#endif
	void GenerateInputRequestedRegion(void);
	void BeforeThreadedGenerateData( void );
	void PrintSelf( std::ostream &os, Indent indent ) const;
private:
	EnhanceImageFilter(const Self&);         // purposely not implemented
	void operator=(const Self&);    // purposely not implemented
	//float ComputeTraceMO0( const InputImageSizeType& rcomp );
	//float ComputeTraceMO1( const InputImageSizeType& rcomp );


	// The standard deviation of noise (in the complex domain)
	float                m_Sigma;
	// The true parameteres of EnhanceImage:
	/*float                m_H;
	float                m_PSTh;
	InputImageSizeType   m_RSearch;
	InputImageSizeType   m_RComp;
	
	// the new Features image
	FeaturesMapPointer   m_Features;*/

	double 		     m_MinimumLevel;
  	double 		     m_MaximumLevel;
  	unsigned int	     m_NumberOfLevelSteps;
	//unsigned int         m_Order;

  	//ScalesImageConstPointer m_scalesImage;
	//bool                    m_GenerateScalesOutput;

	EigenValueImageConstPointer m_strengthImage;
	bool			    m_GenerateStrengthOutput;
	//double		     	    m_Alpha;

	/*bool 		     m_UseEstimatedDistanceMean;
	bool 		     m_NormalizedByFeatureStrength;
	bool 		     m_UseDeltaFeatureStrength;

	MatrixMapType	     m_BMatrixMap;
	MatrixMapType	     m_RXMatrixMap;
	MatrixMapType	     m_rvectorMap;

	IndexVectorMapType   m_IndexVectorMap;
	TraceMapType  	     m_TraceOrder0Map;
	TraceMapType         m_TraceMaxOrderMap;

	BijMapType	     m_BijMatrixMap;	
	IndexBijMapType	     m_BijIndexMap;
	
	BiOrd0MapType        m_BiOrd0Map;
	BijOrd0MapType       m_BijOrd0Map;*/
};


 
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEnhanceImageFilter.txx"
#endif

#endif
