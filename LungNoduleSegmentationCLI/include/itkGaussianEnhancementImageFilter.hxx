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
#ifndef __itkGaussianEnhancementImageFilter_hxx
#define __itkGaussianEnhancementImageFilter_hxx

#include "itkGaussianEnhancementImageFilter.h"
#include "itkImageRegionIterator.h"

namespace itk
{

/**
 * ********************* Constructor ****************************
 */

template < typename TInPixel, typename TOutPixel >
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::GaussianEnhancementImageFilter()
{
  this->m_UnaryFunctor = NULL;
  this->m_UnaryFunctorFilter = UnaryFunctorImageFilterType::New();//needed to be global?
  this->m_Sigma = 1.0;
  this->m_Rescale = true;
  this->m_NormalizeAcrossScale = true;

  // Construct the gradient magnitude filter
  this->m_GradientMagnitudeFilter = GradientMagnitudeFilterType::New();
  this->m_GradientMagnitudeFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

  // Construct the gaussian filter
  this->m_GaussianFilter = GaussianFilterType::New();
  this->m_GaussianFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );
//  this->m_GaussianFilterX = GaussianFilterType::New();
//  this->m_GaussianFilterX->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );
//  this->m_GaussianFilterY = GaussianFilterType::New();
//  this->m_GaussianFilterY->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );
//  this->m_GaussianFilterZ = GaussianLastFilterType::New();
//  this->m_GaussianFilterZ->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

  // Construct the gradient filter
  this->m_GradientFilter = GradientFilterType::New();
  this->m_GradientFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );


  // Construct the Hessian filter
  this->m_HessianFilter = HessianFilterType::New();
  this->m_HessianFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );


  // Construct the eigenvalue filter
  this->m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
  this->m_SymmetricEigenValueFilter->SetDimension( ImageDimension );
  this->m_SymmetricEigenValueFilter->OrderEigenValuesBy(
    EigenAnalysisFilterType::FunctorType::OrderByValue);//OrderByMagnitude?

  // Construct the rescale filter
  this->m_RescaleFilter = RescaleFilterType::New();
  this->m_RescaleFilter->SetOutputMinimum( 0.0 );
  this->m_RescaleFilter->SetOutputMaximum( 1.0 );

  // Allow progressive memory release
  this->m_GradientMagnitudeFilter->ReleaseDataFlagOn();

  this->m_GaussianFilter->ReleaseDataFlagOn();

//  this->m_GaussianFilterX->ReleaseDataFlagOn();
//  this->m_GaussianFilterY->ReleaseDataFlagOn();
//  this->m_GaussianFilterZ->ReleaseDataFlagOn();

  this->m_GradientFilter->ReleaseDataFlagOn();

  this->m_HessianFilter->ReleaseDataFlagOn();

  this->m_SymmetricEigenValueFilter->ReleaseDataFlagOn();
  this->m_RescaleFilter->ReleaseDataFlagOn();

  this->ProcessObject::SetNumberOfRequiredOutputs( 4 );
 
  typename GaussianImageType::Pointer gaussianImage = GaussianImageType::New();
  //  typename GaussianLastImageType::Pointer gaussianImage = GaussianLastImageType::New();
  this->ProcessObject::SetNthOutput( 1, gaussianImage.GetPointer() );

  typename GradientImageType::Pointer gradientImage = GradientImageType::New();
  this->ProcessObject::SetNthOutput( 2, gradientImage.GetPointer() );

  typename HessianTensorImageType::Pointer hessianImage = HessianTensorImageType::New();
  this->ProcessObject::SetNthOutput( 3, hessianImage.GetPointer() );


} // end Constructor


/**
 * ********************* SetUnaryFunctor ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::SetUnaryFunctor( UnaryFunctorBaseType * _arg )
{
  if ( this->m_UnaryFunctor != _arg )
  {
    // Only one of them should be initialized
    this->m_UnaryFunctor = _arg;
    this->m_UnaryFunctorFilter->SetFunctor( _arg );
    this->Modified();
  }
} // end SetUnaryFunctor()

/**
 * ********************* SetNumberOfThreads ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::SetNumberOfThreads( ThreadIdType nt )
{
  Superclass::SetNumberOfThreads( nt );

  this->m_GradientMagnitudeFilter->SetNumberOfThreads( nt );

  this->m_GaussianFilter->SetNumberOfThreads( nt );

//  this->m_GaussianFilterX->SetNumberOfThreads( nt );
//  this->m_GaussianFilterY->SetNumberOfThreads( nt );
//  this->m_GaussianFilterZ->SetNumberOfThreads( nt );

  this->m_GradientFilter->SetNumberOfThreads( nt );

  this->m_HessianFilter->SetNumberOfThreads( nt );

  this->m_SymmetricEigenValueFilter->SetNumberOfThreads( nt );
  this->m_RescaleFilter->SetNumberOfThreads( nt );

  if ( this->m_UnaryFunctorFilter.IsNotNull() )
  {
    this->m_UnaryFunctorFilter->SetNumberOfThreads( nt );
  }

  if ( this->GetNumberOfThreads() != ( nt < 1 ? 1 : ( nt > ITK_MAX_THREADS ? ITK_MAX_THREADS : nt ) ) )
  {
    this->Modified();
  }
} // end SetNumberOfThreads()


/**
 * ********************* SetNormalizeAcrossScale ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::SetNormalizeAcrossScale( bool normalize )
{
  itkDebugMacro( "Setting NormalizeAcrossScale to " << normalize );
  if( this->m_NormalizeAcrossScale != normalize )
  {
    this->m_NormalizeAcrossScale = normalize;

    this->m_GradientMagnitudeFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

    this->m_GaussianFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

//    this->m_GaussianFilterX->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );
 //   this->m_GaussianFilterY->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );
//    this->m_GaussianFilterZ->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

    this->m_GradienteFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

    this->m_HessianFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

    this->Modified();
  }
} // end SetNormalizeAcrossScale()

/**
 * ********************* GenerateData ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::GenerateData( void )
{
  if ( this->m_UnaryFunctor.IsNull() )
  {
    itkExceptionMacro( << "ERROR: Missing Functor. "
      << "Please provide functor for multi scale framework." );
  }
  
// Calculate the gradient image.
//  this->m_GaussianFilterX->SetInput( this->GetInput() );
//  this->m_GaussianFilterX->SetSigma( this->m_Sigma );
//  this->m_GaussianFilterX->SetDirection( 0 ); // X direction
//
//  this->m_GaussianFilterY->SetInput( this->m_GaussianFilterX->GetOutput() );
//  this->m_GaussianFilterY->SetSigma( this->m_Sigma );
//  this->m_GaussianFilterY->SetDirection( 1 ); // Y direction
//
//  this->m_GaussianFilterZ->SetInput( this->m_GaussianFilterY->GetOutput() );
//  this->m_GaussianFilterZ->SetSigma( this->m_Sigma );
//  this->m_GaussianFilterZ->SetDirection( 2 ); // Z direction
//
//  this->m_GaussianFilterZ->Update();
  this->m_GaussianFilter->SetInput( this->GetInput() );
  this->m_GaussianFilter->SetSigma( this->m_Sigma );
  this->m_GaussianFilter->Update();

  this->GraftNthOutput(1, this->m_GaussianFilter->GetOutput() );

  // Calculate the gradient image.
  this->m_GradientFilter->SetInput( this->GetInput() );
  this->m_GradientFilter->SetSigma( this->m_Sigma );
  this->m_GradientFilter->Update();
  this->GraftNthOutput(2, this->m_GradientFilter->GetOutput() );

  // Calculate the hessian image.
  this->m_HessianFilter->SetInput( this->GetInput() );
  this->m_HessianFilter->SetSigma( this->m_Sigma );
  this->m_HessianFilter->Update();
  this->GraftNthOutput(3, this->m_HessianFilter->GetOutput() );

  this->m_SymmetricEigenValueFilter->SetInput( this->m_HessianFilter->GetOutput() );  
  this->m_SymmetricEigenValueFilter->Update();
  
  // Calculate unary functor filter.
  this->m_UnaryFunctorFilter->SetInput( this->m_SymmetricEigenValueFilter->GetOutput() );
  this->m_UnaryFunctorFilter->Update();

  // Apply rescale
  if( this->m_Rescale )
  {
    // Rescale the output to [0,1].
    this->m_RescaleFilter->SetInput( this->m_UnaryFunctorFilter->GetOutput() );
    this->m_RescaleFilter->Update();

    // Put the output of the rescale filter to this filter's output.
    this->GraftOutput( this->m_RescaleFilter->GetOutput() );
  }
  else
  {
    this->GraftOutput( this->m_UnaryFunctorFilter->GetOutput() );
  }
} // end GenerateData()


/**
 * ********************* GetGaussianImage ****************************
 */
template < typename TInPixel, typename TOutPixel >
const typename GaussianEnhancementImageFilter< TInPixel, TOutPixel >::GaussianImageType *
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::GetGaussianImage( void ) const
{ 
	return static_cast<const GaussianImageType*>(this->ProcessObject::GetOutput(1));
} 



/**
 * ********************* GetGradientImage ****************************
 */
template < typename TInPixel, typename TOutPixel >
const typename GaussianEnhancementImageFilter< TInPixel, TOutPixel >::GradientImageType *
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::GetGradientImage( void ) const
{ 
	return static_cast<const GradientImageType*>(this->ProcessObject::GetOutput(2));
} 



/**
 * ********************* GetHessianImage ****************************
 */
template < typename TInPixel, typename TOutPixel >
const typename GaussianEnhancementImageFilter< TInPixel, TOutPixel >::HessianTensorImageType *
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::GetHessianImage( void ) const
{ 
	return static_cast<const HessianTensorImageType*>(this->ProcessObject::GetOutput(3));
} 



/**
 * ********************* PrintSelf ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
GaussianEnhancementImageFilter< TInPixel, TOutPixel >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Sigma: " << this->m_Sigma << std::endl;
  os << indent << "Rescale: " << this->m_Rescale << std::endl;
  os << indent << "NormalizeAcrossScale: " << this->m_NormalizeAcrossScale << std::endl;

  Indent nextIndent = indent.GetNextIndent();

  this->m_UnaryFunctorFilter->Print( os, nextIndent );

} // end PrintSelf()


} // end namespace itk

#endif
