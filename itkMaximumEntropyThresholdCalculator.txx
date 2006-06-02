/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMaximumEntropyThresholdCalculator.txx,v $
  Language:  C++
  Date:      $Date: 2005/01/13 15:36:46 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMaximumEntropyThresholdCalculator_txx
#define _itkMaximumEntropyThresholdCalculator_txx

#include "itkMaximumEntropyThresholdCalculator.h"

namespace itk
{


/*
 * Compute Otsu's thresholds
 */                    
template<class TInputHistogram>
void
MaximumEntropyThresholdCalculator<TInputHistogram>
::GenerateData()
{
  typename TInputHistogram::ConstPointer histogram = this->GetInputHistogram();

  // TODO: as an improvement, the class could accept multi-dimensional histograms
  // and the user could specify the dimension to apply the algorithm to.
  if (histogram->GetSize().GetSizeDimension() != 1)
    {
    itkExceptionMacro(<<"Histogram must be 1-dimensional.");
    }

  // create a vector to store the entropy, and fill it with zeros
  typedef typename std::vector< double > EntropyListType;
  EntropyListType entropyList;
  entropyList.resize( histogram->Size() );

  // entropy for white and black part of the histogram
  double pt = 0;
  for( unsigned long i=0; i<histogram->Size(); i++)
    {
    double prob = histogram->GetFrequency( i ) / histogram->GetTotalFrequency();
    pt += prob;

    if( pt > 0 )
      {
      double hhB = 0;
      for( int j=0; j<=i; j++ )
        {
        double probj = histogram->GetFrequency( j ) / histogram->GetTotalFrequency();
        if( probj > 0 )
          {
          hhB -= probj / pt * vcl_log( probj / pt );
          }
        }
        entropyList[i] = hhB;
      }
    else
      { entropyList[i] = 0; }

    double ptw = 1 - pt;
    if( ptw > 0 )
      {
      double hhW = 0;
      for( int j=i+1; j<histogram->Size(); j++ )
        {
        double probj = histogram->GetFrequency( j ) / histogram->GetTotalFrequency();
        if( probj > 0 )
          {
          hhW -= probj / ptw * vcl_log( probj / ptw );
          }
        }
        entropyList[i] += hhW;
      }
    
    }


// std::cout << std::endl<< histogram->GetTotalFrequency()<< std::endl<< std::endl;
// std::cout << std::endl<< histogram->GetSize()<< std::endl<< std::endl;
// std::cout << std::endl<< histogram->GetSize(0)<< std::endl<< std::endl;

  // find index with maximum entropy
  double max = entropyList[0];
  unsigned long idx = 0;
  for(unsigned long i = 0; i<entropyList.size(); i++)
    {
/*std::cout << histogram->GetMeasurementVector(i)[0]+0.0 << ":\t" << entropyList[i] << "\t" << histogram->GetFrequency(i) << "\t" << histogram->GetBinMin(0,i)+0.0 << "\t" << histogram->GetBinMax(0,i)+0.0 << std::endl;*/
    if( entropyList[i] > max )
      {
      max = entropyList[i];
      idx = i;
      }
    }

  m_Output = histogram->GetMeasurementVector( idx )[0];
}

template<class TInputHistogram>
void
MaximumEntropyThresholdCalculator<TInputHistogram>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Output: " << static_cast<typename NumericTraits<MeasurementType>::PrintType>(m_Output) << std::endl;

}

} // end namespace itk

#endif
