/*=========================================================================
*
* Copyright Marius Staring, Stefan Klein, David Doria. 2011.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0.txt
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
*=========================================================================*/
#ifndef __itkComputeFeaturesFilter_h
#define __itkComputeFeaturesFilter_h

// ITK include files
#include "itkImageToImageFilter.h"

#include "itkUnaryFunctorBase.h"
#include "itkUnaryFunctorImageFilter2.h"

#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"

//#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

namespace itk
{

/** \class ComputeFeaturesFilter
 * \brief A filter to compute Hessian
 *      and Gradient measures of an image in a single scale framework.
 */

template < typename TInputImage, typename TOutputImage >
class ComputeFeaturesFilter
  : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef ComputeFeaturesFilter                    	    Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage >   Superclass;
  typedef SmartPointer<Self>                                Pointer;
  typedef SmartPointer<const Self>                          ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro( ComputeFeaturesFilter, ImageToImageFilter );

  /** Method for creation through the object factory.*/
  itkNewMacro( Self );

  /** Typedef's. */
  typedef TInputImage                               	InputImageType;
  typedef TOutputImage                              	OutputImageType;
  typedef typename InputImageType::PixelType        	InputPixelType;
  typedef typename OutputImageType::PixelType       	OutputPixelType;

  typedef typename NumericTraits<OutputPixelType>::RealType RealType;

  /** Image dimension = 3. */
  itkStaticConstMacro( ImageDimension, unsigned int, InputImageType::ImageDimension );

  /** Gradient features filter type, since e.g. the strain energy vesselness
   * is not only a function of the Hessian, but also of the first order derivatives.
   
  typedef OutputPixelType                         	GradientPixelType;
  typedef Image< 
    FixedArray< OutputPixelType,
    itkGetStaticConstMacro( ImageDimension ) >,
    itkGetStaticConstMacro( ImageDimension ) > 		GradientImageType;
  typedef GradientRecursiveGaussianImageFilter<
    InputImageType, GradientImageType >  		GradientFilterType;*/
  
  /** Hessian filter type */
  typedef Image<SymmetricSecondRankTensor<
    OutputPixelType,  3 >,  3 >    			HessianTensorImageType;
  typedef HessianRecursiveGaussianImageFilter<
    InputImageType, HessianTensorImageType >      	HessianFilterType;

  /** EigenValue analysis filter */
  typedef FixedArray< OutputPixelType,  3 >    		EigenValueArrayType;
  typedef Image< EigenValueArrayType,   3 >    		EigenValueImageType;
  typedef SymmetricEigenAnalysisImageFilter<
    HessianTensorImageType, EigenValueImageType > 	EigenAnalysisFilterType;

  /** Rescale filter type */
  typedef RescaleIntensityImageFilter<
    OutputImageType, OutputImageType >            	RescaleFilterType;

  /** Unary functor filter type */
  typedef UnaryFunctorImageFilter2<
    EigenValueImageType,
    OutputImageType >                             		UnaryFunctorImageFilterType;
  typedef typename UnaryFunctorImageFilterType::FunctorType 	UnaryFunctorBaseType;

  /** Set/Get unary functor filter */
  virtual void SetUnaryFunctor( UnaryFunctorBaseType * _arg );
  itkGetObjectMacro( UnaryFunctor, UnaryFunctorBaseType );

  /** Set/Get macros for Sigma. The current scale used.
   * Sigma should be positive. */
  itkSetClampMacro( Sigma, double, 0.0, NumericTraits<double>::max() );
  itkGetConstReferenceMacro( Sigma, double );

  /** Methods to turn on/off flag to rescale function output */
  itkSetMacro( Rescale, bool );
  itkGetConstMacro( Rescale, bool );
  itkBooleanMacro( Rescale );

  /** Define whether or not normalization factor will be used for the Gaussian. default true */
  void SetNormalizeAcrossScale( bool normalize );
  itkGetConstMacro( NormalizeAcrossScale, bool );

  /** Set the number of threads to create when executing. */
  void SetNumberOfThreads( ThreadIdType nt );

protected:
  ComputeFeaturesFilter();
  virtual ~ComputeFeaturesFilter() {};

  virtual void PrintSelf(std::ostream& os, Indent indent) const;
  virtual void GenerateData( void );

private:
  ComputeFeaturesFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Member variables. */
  //typename GradientFilterType::Pointer		m_GradientFilter;

  typename HessianFilterType::Pointer		m_HessianFilter;
  typename EigenAnalysisFilterType::Pointer	m_SymmetricEigenValueFilter;
  typename RescaleFilterType::Pointer		m_RescaleFilter;

  typename UnaryFunctorBaseType::Pointer 	m_UnaryFunctor;
  typename UnaryFunctorImageFilterType::Pointer m_UnaryFunctorFilter;

  double  m_Sigma;
  bool    m_Rescale;
  bool    m_NormalizeAcrossScale; // Normalize the image across scale space
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComputeFeaturesFilter.hxx"
#endif

#endif
