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
// This file contains the
// ==========================================================================

#ifndef _INCLUDE_VIENNA_RNA_H_
#define _INCLUDE_VIENNA_RNA_H_

// ----------------------------------------------------------------------------
// Vienna headers
// ----------------------------------------------------------------------------

extern "C" {
    #include  <ViennaRNA/data_structures.h>
    #include  <ViennaRNA/params.h>
    #include  <ViennaRNA/utils.h>
    #include  <ViennaRNA/eval.h>
    #include  <ViennaRNA/fold.h>
    #include  <ViennaRNA/part_func.h>
    #include <ViennaRNA/PS_dot.h>
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function computeBppMatrix()
// ----------------------------------------------------------------------------

template <typename TOption, typename TRnaStruct>
void computeBppMatrix(TOption const & options, TRnaStruct & rnaSeq)
{
    std::cout << length(rnaSeq.seq) << " " << rnaSeq.seq << std::endl;
    char *structure = new char[length(rnaSeq.seq) + 1];
    vrna_md_t md_p;

//  apply default model details
    set_model_details(&md_p);

    // Create a c-style string object for str:
    String<char, CStyle> seq;

    seq = rnaSeq.seq;
//  get a vrna_fold_compound with MFE and PF DP matrices and default model details
    vrna_fold_compound_t *vc = vrna_fold_compound(toCString(seq), &md_p, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    double gibbs = (double)vrna_pf(vc, structure); //FIXME the structure is not well saved

    vrna_plist_t *pl1, *ptr;
    pl1 = vrna_plist_from_probs(vc, options.thrBppm);
/*
    if(options.verbose>1)
    {
        vrna_plist_t *pl2;
        pl2= vrna_plist(structure, 0.95*0.95);
//	Function used to plot the dot_plot graph
        (void) PS_dot_plot_list(toCString(seq), "prova_dot_plot", pl1, pl2, "");
    }
*/
// get size of pl1
    unsigned size;
    for(size = 0, ptr = pl1; ptr->i; size++, ptr++);

    std::cout << "BPPM2 => " << size  << std::endl;
    std::cout << "BPPM2 => " << std::endl;

    std::cout << "size seq = " << length(rnaSeq.seq) << std::endl;
// std::cout << "size graph = " << length(interGraph) << std::endl;
// std::cout << interGraph << std::endl;

    for(unsigned i=0; i<length(rnaSeq.seq);++i)
    {
        addVertex(rnaSeq.graph);
    }
    for(unsigned i=0; i<size;++i)
    {
        if(options.verbose > 2)
            std::cout << i << "_"<< pl1[i].i <<":"<< pl1[i].j <<"|"<< pl1[i].p <<"|"<< pl1[i].type << "\t";
        addEdge(rnaSeq.graph, pl1[i].i, pl1[i].j, pl1[i].p);
    }
    if(options.verbose > 2)
        std::cout << "\n" << rnaSeq.graph << std::endl;


    std::cout << rnaSeq.seq << std::endl;
    std::cout << structure << "\tgibbs = " << gibbs << std::endl;

//	free memory occupied by vrna_fold_compound
    vrna_fold_compound_free(vc);
//	clean up
    free(structure);
}

#endif //_INCLUDE_VIENNA_RNA_H_