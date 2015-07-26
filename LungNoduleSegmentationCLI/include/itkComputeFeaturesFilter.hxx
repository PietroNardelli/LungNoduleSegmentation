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
#ifndef __itkComputeFeaturesFilter_hxx
#define __itkComputeFeaturesFilter_hxx

#include "itkComputeFeaturesFilter.h"
#include "itkImageRegionIterator.h"

namespace itk
{

/**
 * ********************* Constructor ****************************
 */

template < typename TInPixel, typename TOutPixel >
ComputeFeaturesFilter< TInPixel, TOutPixel >
::ComputeFeaturesFilter()
{
  this->m_UnaryFunctor = NULL;
  this->m_UnaryFunctorFilter = UnaryFunctorImageFilterType::New();//needed to be global?
  this->m_Sigma = 1.0;
  this->m_Rescale = true;
  this->m_NormalizeAcrossScale = true;

  // Construct the gradient filter
  //this->m_GradientFilter = GradientFilterType::New();  //need of this??
  //this->m_GradientFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

  // Construct the Hessian filter
  this->m_HessianFilter = HessianFilterType::New();
  this->m_HessianFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

  // Construct the eigenvalue filter
  this->m_SymmetricEigenValueFilter = EigenAnalysisFilterType::New();
  this->m_SymmetricEigenValueFilter->SetDimension( 3 );
  this->m_SymmetricEigenValueFilter->OrderEigenValuesBy(
    EigenAnalysisFilterType::FunctorType::OrderByValue );//OrderByMagnitude?

  // Construct the rescale filter
  this->m_RescaleFilter = RescaleFilterType::New();
  this->m_RescaleFilter->SetOutputMinimum( 0.0 );
  this->m_RescaleFilter->SetOutputMaximum( 1.0 );

  // Allow progressive memory release
  this->m_HessianFilter->ReleaseDataFlagOn();
  //this->m_GradientFilter->ReleaseDataFlagOn();
  this->m_SymmetricEigenValueFilter->ReleaseDataFlagOn();
  this->m_RescaleFilter->ReleaseDataFlagOn();

} // end Constructor


/**
 * ********************* SetUnaryFunctor ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
ComputeFeaturesFilter< TInPixel, TOutPixel >
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
ComputeFeaturesFilter< TInPixel, TOutPixel >
::SetNumberOfThreads( ThreadIdType nt )
{
  Superclass::SetNumberOfThreads( nt );

  //this->m_GradientFilter->SetNumberOfThreads( nt );
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
ComputeFeaturesFilter< TInPixel, TOutPixel >
::SetNormalizeAcrossScale( bool normalize )
{
  itkDebugMacro( "Setting NormalizeAcrossScale to " << normalize );
  if( this->m_NormalizeAcrossScale != normalize )
  {
    this->m_NormalizeAcrossScale = normalize;

    //this->m_GradienteFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

    this->m_HessianFilter->SetNormalizeAcrossScale( this->m_NormalizeAcrossScale );

    this->Modified();
  }
} // end SetNormalizeAcrossScale()


/**
 * ********************* GenerateData ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
ComputeFeaturesFilter< TInPixel, TOutPixel >
::GenerateData( void )
{
  if ( this->m_UnaryFunctor.IsNull() )
  {
    itkExceptionMacro( << "ERROR: Missing Functor. "
      << "Please provide functor for multi scale framework." );
  }

  /* Gradient features */
  /*this->m_GradientFilter->SetInput( this->GetInput() );
  this->m_GradientFilter->SetSigma( this->m_Sigma );
  this->m_GradientFilter->Update();

  typedef typename GradientImageType::RegionType::IndexType inputStartType1;
  typedef typename GradientImageType::RegionType::SizeType inputSizeType1;
  inputStartType1 inputStart1;
  inputSizeType1 size1;
  inputStart1[0] = 198;
  inputStart1[1] = 231;
  inputStart1[2] = 8;
  size1[0]=5;
  size1[1]=5;
  size1[2]=5;
  typedef typename GradientImageType::RegionType inputRegionType1;
  inputRegionType1 inputRegion1;
  inputRegion1.SetSize(size1);
  inputRegion1.SetIndex(inputStart1);
  itk::ImageRegionConstIterator<GradientImageType>      it1;
  it1  = itk::ImageRegionConstIterator<GradientImageType>( this->m_GradientFilter->GetOutput() ,    inputRegion1 );
  it1.GoToBegin();
  std::cout<<"value: "<<it1.Get()<<std::endl;*/

  // Calculate the eigenvalue vector image.
  this->m_HessianFilter->SetInput( this->GetInput() );
  this->m_HessianFilter->SetSigma( this->m_Sigma );
  this->m_HessianFilter->Update();

  this->m_SymmetricEigenValueFilter->SetInput( this->m_HessianFilter->GetOutput() );  
  this->m_SymmetricEigenValueFilter->Update();

  //trying
  /*typedef typename HessianTensorImageType::RegionType::IndexType inputStartType;
  typedef typename HessianTensorImageType::RegionType::SizeType inputSizeType;
  inputStartType inputStart;
  inputSizeType size;
  inputStart[0] = 198;
  inputStart[1] = 231;
  inputStart[2] = 8;
  size[0]=5;
  size[1]=5;
  size[2]=5;
  typedef typename HessianTensorImageType::RegionType inputRegionType;
  inputRegionType inputRegion;
  inputRegion.SetSize(size);
  inputRegion.SetIndex(inputStart);
  itk::ImageRegionConstIterator<HessianTensorImageType>      it;
  this->m_HessianFilter->Update(); 
  it  = itk::ImageRegionConstIterator<HessianTensorImageType>( this->m_HessianFilter->GetOutput() ,    inputRegion );
  it.GoToBegin();
  std::cout<<"value: "<<it.Get()<<std::endl;*/

  if ( this->m_UnaryFunctor.IsNotNull() )
  {
    // Calculate unary functor filter.
    this->m_UnaryFunctorFilter->SetInput(
      this->m_SymmetricEigenValueFilter->GetOutput() );
    this->m_UnaryFunctorFilter->Update();
  }
  // Apply rescale
  if( this->m_Rescale )
  {
    // Rescale the output to [0,1].
    if ( this->m_UnaryFunctor.IsNotNull() )
    {
      this->m_RescaleFilter->SetInput( this->m_UnaryFunctorFilter->GetOutput() );
    }
    this->m_RescaleFilter->Update();

    // Put the output of the rescale filter to this filter's output.
    this->GraftOutput( this->m_RescaleFilter->GetOutput() );
  }
  else
  {
    if ( this->m_UnaryFunctor.IsNotNull() )
    {
      this->GraftOutput( this->m_UnaryFunctorFilter->GetOutput() );
    }
  }
} // end GenerateData()


/**
 * ********************* PrintSelf ****************************
 */

template < typename TInPixel, typename TOutPixel >
void
ComputeFeaturesFilter< TInPixel, TOutPixel >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Sigma: " << this->m_Sigma << std::endl;
  os << indent << "Rescale: " << this->m_Rescale << std::endl;
  os << indent << "NormalizeAcrossScale: " << this->m_NormalizeAcrossScale << std::endl;

  Indent nextIndent = indent.GetNextIndent();
  if ( this->m_UnaryFunctorFilter.IsNotNull() )
  {
    this->m_UnaryFunctorFilter->Print( os, nextIndent );
  }
} // end PrintSelf()


} // end namespace itk

#endif
