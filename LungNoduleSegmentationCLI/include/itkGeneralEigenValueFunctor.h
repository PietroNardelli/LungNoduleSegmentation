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
#ifndef __itkGeneralEigenValueFunctor_h
#define __itkGeneralEigenValueFunctor_h

#include "itkUnaryFunctorBase.h"
#include "itkComparisonOperators.h"
#include "vnl/vnl_math.h"

namespace itk
{
namespace Functor
{

/** \class GeneralEigenValueFunctor
 * \brief Computes a measure the Hessian eigenvalues.
 *
 *
 * \sa FrangiVesselnessImageFilter
 * \ingroup IntensityImageFilters Multithreaded
 */

template< class TInput, class TOutput >
class GeneralEigenValueFunctor
  : public UnaryFunctorBase< TInput, TOutput >
{
public:
  /** Standard class typedefs. */
  typedef GeneralEigenValueFunctor            Self;
  typedef UnaryFunctorBase< TInput, TOutput > Superclass;
  typedef SmartPointer< Self >                Pointer;
  typedef SmartPointer< const Self >          ConstPointer;

  /** New macro for creation of through a smart pointer. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GeneralEigenValueFunctor, UnaryFunctorBase );

  /** Typedef's. */
  typedef typename NumericTraits<TOutput>::RealType RealType;
  typedef TInput                                    EigenValueArrayType;
  typedef typename EigenValueArrayType::ValueType   EigenValueType;

  /** This does the real computation */
  virtual TOutput Evaluate( const TInput & eigenValues ) const
  {
    /** Sort the eigenvalues by their absolute value, such that |l1| < |l2| < |l3|. */
    EigenValueArrayType sortedEigenValues = eigenValues;
    std::sort( sortedEigenValues.Begin(), sortedEigenValues.End(),
    Functor::AbsLessCompare<EigenValueType>() );

    const RealType l1 = sortedEigenValues[ 0 ];
    const RealType l2 = sortedEigenValues[ 1 ];
    const RealType l3 = sortedEigenValues[ 2 ];
    
    /** Take the absolute values and abbreviate. */
    const RealType abs_l1 = vnl_math_abs( sortedEigenValues[ 0 ] );
    const RealType abs_l2 = vnl_math_abs( sortedEigenValues[ 1 ] );
    const RealType abs_l3 = vnl_math_abs( sortedEigenValues[ 2 ] );

    /** Check Frobenius Norm. This relates to the background noise*/

    //const RealType S  = vcl_sqrt( l1 * l1 + l2 * l2 + l3 * l3 );


    /** Avoid divisions by zero (or close to zero). */
    RealType S;
    if( l1 < 0.0 && l2 < 0.0 && l3 < 0.0 )    
    {
        S = (abs_l1*abs_l1)/abs_l3; // Li's method
        //S = abs_l1/vcl_sqrt( abs_l2 * abs_l3 ); //Frangi's method
        //S = abs_l3*vcl_sqrt(l1/l3); //Sato's method
    }  
    else
    {
      	S = NumericTraits<TOutput>::Zero;
    }

    return S;

  } // end operator ()

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( DimensionIs3Check,
    ( Concept::SameDimension< EigenValueArrayType::Dimension, 3 > ) );
  /** End concept checking */
#endif

protected:
  /** Constructor */
  GeneralEigenValueFunctor(){};
  virtual ~GeneralEigenValueFunctor(){};

private:
  GeneralEigenValueFunctor(const Self &); // purposely not implemented
  void operator=(const Self &);    // purposely not implemented
}; // end class GeneralEigenValueFunctor

} // end namespace itk::Functor
} // end namespace itk

#endif
