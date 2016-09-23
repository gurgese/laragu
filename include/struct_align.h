// ==========================================================================
//               LaRAgu - Lagrangian Relaxation Aligner GU
// ==========================================================================
// Copyright (c) 2015-2016, Gianvito Urgese
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Gianvito Urgese nor the names of its contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL GIANVITO URGESE OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================
// This file contains
// ==========================================================================
#ifndef _INCLUDE_STRUCT_ALIGN_H_
#define _INCLUDE_STRUCT_ALIGN_H_

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

//#include "vienna_rna.h"

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function alignVectorBuild()
// ----------------------------------------------------------------------------
// Used to generate the alignments from a single input file
template <typename TRnaAligns, typename TRnaSeqs, typename TOption>
void alignVectorBuild(TRnaAligns & rnaAligns, TRnaSeqs const & rnaSeqs,
                      TOption const & options)
{
    for (unsigned i = 0; i < length(rnaSeqs) - 1; ++i) {
        TRnaAlign rnaAlign;
        for (unsigned j = i + 1; j < length(rnaSeqs); ++j) {
            if (length(rnaSeqs[i].seq) <
                length(rnaSeqs[j].seq)) // in this way the alignment map structure will be always created with the maximum size
            {
                rnaAlign.rna1 = rnaSeqs[j];
                rnaAlign.rna2 = rnaSeqs[i];
            } else {
                rnaAlign.rna1 = rnaSeqs[i];
                rnaAlign.rna2 = rnaSeqs[j];
            }
            if(options.verbose > 2)
            {
                std::cout << rnaAlign.rna1.seq << std::endl;
                std::cout << rnaAlign.rna2.seq << std::endl;
            }
            rnaAligns.push_back(rnaAlign);
        }
    }
}


// ----------------------------------------------------------------------------
// Function alignVectorBuild()
// ----------------------------------------------------------------------------
// Used to generate the alignments from two different input files
template <typename TRnaAligns, typename TRnaSeqs, typename TOption>
void alignVectorBuild(TRnaAligns & rnaAligns, TRnaSeqs const & rnaSeqs,
                      TRnaSeqs const & rnaSeqsRef, TOption const & options)
{
    for(unsigned i=0;i<length(rnaSeqs); ++i)
    {
        TRnaAlign rnaAlign;
        for(unsigned j=0;j<length(rnaSeqsRef); ++j)
        {
            if(length(rnaSeqs[i].seq) < length(rnaSeqsRef[j].seq)) // in this way the alignment map structure will be always created with the maximum size
            {
                rnaAlign.rna1 = rnaSeqsRef[j];
                rnaAlign.rna2 = rnaSeqs[i];
            }else {
                rnaAlign.rna1 = rnaSeqs[i];
                rnaAlign.rna2 = rnaSeqsRef[j];
            }
            if(options.verbose > 2)
            {
                std::cout << rnaAlign.rna1.seq << std::endl;
                std::cout << rnaAlign.rna2.seq << std::endl;
            }
        }
        rnaAligns.push_back(rnaAlign);
    }
}


#endif //_INCLUDE_STRUCT_ALIGN_H_