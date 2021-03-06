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
// This file contains the SeqStore class.
// ==========================================================================

#ifndef _INCLUDE_STORE_SEQS_H_
#define _INCLUDE_STORE_SEQS_H_

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/seq_io.h>
#include <seqan/rna_io.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function CheckInstanceFormat();
// ----------------------------------------------------------------------------

template <typename TOption, typename TPath>
unsigned CheckInstanceFormat(TOption const & options, TPath & inFilePath) // TODO to be deprecated after all the requested file formats will be available in SeqAn
{
    _V(options, "Check the Instance Format and save the input file code")
    unsigned input_type = UNKNOWN;
    std::ifstream inFile(toCString(inFilePath));
    _V(options, "input file ifstream " << inFilePath);
    std::string cur_line;
    if(inFile.is_open())
    {
        std::getline(inFile, cur_line);
/*
        do
        {
            std::getline(inFile, cur_line);
        } while(cur_line[0] == '#' && !inFile.eof());
        if(cur_line.size() == 0)
            input_type = UNKNOWN;
*/
        if(cur_line[0] == '>' )
        {
            input_type = FASTA;
            do
            {
                std::getline(inFile, cur_line);
                if (cur_line.find("(") != -1 && cur_line[0] != '>') {
                    if (cur_line.find("[") != -1 && cur_line[0] != '>') // this is not exaustive because the EDBN file format should be defined yet (we still do not know if letters or parenthesis should be used for represent the nested loops)
                        input_type = EDBN;
                    else
                        input_type = DBN;
                }
            } while(!inFile.eof());
// if we did not find structure, then return FASTA, otherwis we already returned EXTENDED_FASTA
        } else if(cur_line[0] == '@' )
        {
            input_type = FASTQ;
        } else if(cur_line[0] == '#' && cur_line[1] == '#')
        {
            input_type = EDBN;
        } else if(cur_line[0] == '#' && cur_line[1] != '#')
        {
            input_type = DBN;
        }
    }
    inFile.close();
    return input_type;
}

// ----------------------------------------------------------------------------
// Function readFastaRecords()
// ----------------------------------------------------------------------------
// This function is able to read Fasta, MultiFasta, FASTQ, EMBL or GenBank formats.
template <typename TOption, typename TSeqVect, typename TSeqFileIn>
unsigned readFastaRecords(TOption const & options, TSeqVect & rnaSeqs, TSeqFileIn & seqFileIn)
{
    StringSet<CharString> ids;
    StringSet<Rna5String> seqs;
    StringSet<CharString> quals;
    try
    {
        readRecords(ids, seqs, quals, seqFileIn);
    }
    catch (Exception const &e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    seqan::resize(rnaSeqs, length(ids));
    if(length(quals) == length(ids))
    {
        for (unsigned i = 0; i < length(ids); ++i)
        {
            rnaSeqs[i].name = ids[i];
            rnaSeqs[i].seq = seqs[i];
            rnaSeqs[i].qual = quals[i];
            _VV(options, rnaSeqs[i].name << "\t" << rnaSeqs[i].seq << "\t" << rnaSeqs[i].qual);
        }
    } else {
        for (unsigned i = 0; i < length(ids); ++i)
        {
            rnaSeqs[i].name = ids[i];
            rnaSeqs[i].seq = seqs[i];
            _VV(options, rnaSeqs[i].name << "\t" << rnaSeqs[i].seq);
        }
    }
    return 0;
};

// ----------------------------------------------------------------------------
// Function readDbnRecords()
// ----------------------------------------------------------------------------
// This function is able to read dbn format.
template <typename TOption, typename TSeqVect, typename TSeqFileIn>
unsigned readDbnRecords(TOption const & options, TSeqVect & rnaSeqs, TSeqFileIn & seqFileIn)
{
    std::cout << "function readDbnRecords to be implemented" <<std::endl;
};

// ----------------------------------------------------------------------------
// Function readEdbnRecords()
// ----------------------------------------------------------------------------
// This function is able to read edbn format.
template <typename TOption, typename TSeqVect, typename TSeqFileIn>
unsigned readEdbnRecords(TOption const & options, TSeqVect & rnaSeqs, TSeqFileIn & seqFileIn)
{
    std::cout << "function readEdbnRecords to be implemented" <<std::endl;
};

// ----------------------------------------------------------------------------
// Function readBpseqRecords()
// ----------------------------------------------------------------------------
// This function is able to read edbn format.
template <typename TOption, typename TSeqVect, typename TSeqFileIn>
unsigned readBpseqRecords(TOption const & options, TSeqVect & rnaSeqs, TSeqFileIn & seqFileIn)
{
    std::cout << "function readBpseqRecords to be implemented" <<std::endl;
};

// ----------------------------------------------------------------------------
// Function readEbpseqRecords()
// ----------------------------------------------------------------------------
// This function is able to read edbn format.
template <typename TOption, typename TSeqVect, typename TSeqFileIn>
unsigned readEbpseqRecords(TOption const & options, TSeqVect & rnaSeqs, TSeqFileIn & seqFileIn)
{
    std::cout << "function readEbpseqRecords to be implemented" <<std::endl;
};


// ----------------------------------------------------------------------------
// Function readRnaRecords()
// ----------------------------------------------------------------------------
// This function is able to read Fasta, MultiFasta, FASTQ, EMBL, GenBank, dbn, edbn, bpseq and ebpseq formats.
template <typename TOption, typename TFile, typename TSeqVect>
unsigned readRnaRecords(TOption const & options, TFile file, TSeqVect & rnaSeqs)
{
    unsigned seqanSupportedFiles = UNKNOWN;
    const char * fileName = toCString(file);
    CharString inFilePath = getAbsolutePath(fileName);
    SeqFileIn seqFileIn; // TODO the definition of SeqFileIn, and all the relative functions, should be overloded with the new supported RNA file formats (jump to SeqFileIn and overload this function)
    std::cout << "Open the file, recognize the file format and fill the RNA data structure" << std::endl;
    std::cout << "fasta, fastq, bpseq, ebpseq, dbn and edbn(extended dbn) should be supported" << std::endl;
    seqanSupportedFiles = CheckInstanceFormat(options, inFilePath);  //TODO this function must be modified in order to support the acquisition of the other input file formats
    _V(options,"Reading sequences from file type " << seqanSupportedFiles << " named " << inFilePath );

    if (!open(seqFileIn, toCString(inFilePath))) //FIXME this works for FASTA and FASTQ only, because seqFileIn iso not supported for the other file formats
    {
        std::cerr << "ERROR: Could not open the file " << file << std::endl;
        return 1;
    }
    std::cout << "the file format should be recognized and a readRecord function with the file type flag should be called to acquire the inputs" <<std::endl;
    switch(seqanSupportedFiles)
    {
        case FASTA:
            _V(options, "Input file is in Fasta format");
            readFastaRecords(options, rnaSeqs, seqFileIn);
            break;
        case FASTQ:
            _V(options, "Input file is in Fastq format");
            readFastaRecords(options, rnaSeqs, seqFileIn);
            break;
        case DBN:
            _V(options, "Input file is in DBN format");
            readDbnRecords(options, rnaSeqs, seqFileIn);
            break;
        case EDBN:
            _V(options, "Input file is in EDBN format");
            readEdbnRecords(options, rnaSeqs, seqFileIn);
            break;
        case BPSEQ:
        _V(options, "Input file is in BPSEQ format");
            readBpseqRecords(options, rnaSeqs, seqFileIn);
            break;
        case EBPSEQ:
        _V(options, "Input file is in EBPSEQ format");
            readEbpseqRecords(options, rnaSeqs, seqFileIn);
            break;
        case UNKNOWN:
            _V(options, "Unable to identify the type of your input instance. Please check the allowed file formats");
            _V(options, "Please check, if you supplied a either valid fasta/extended fasta file or a dotplot instance.");
            exit(1);
    }


    return 0;
}

/*
// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clear(SeqStore<TSpec, TConfig> & me)
{
    clear(me.seqs);
    clear(me.names);
}

// ----------------------------------------------------------------------------
// Function shrinkToFit()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void shrinkToFit(SeqStore<TSpec, TConfig> & me)
{
    shrinkToFit(me.seqs);
    shrinkToFit(me.names);
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
inline void reserve(SeqStore<TSpec, TConfig> & me, TSize newCapacity)
{
    reserve(me.seqs, newCapacity, Exact());
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileName>
inline bool open(SeqStore<TSpec, TConfig> & me, TFileName const & fileName, int openMode)
{
    CharString name;

    name = fileName;    append(name, ".txt");
    if (!open(me.seqs, toCString(name), openMode)) return false;

    name = fileName;    append(name, ".rid");
    if (!open(me.names, toCString(name), openMode)) return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileName>
inline bool save(SeqStore<TSpec, TConfig> const & me, TFileName const & fileName)
{
    CharString name;

    name = fileName;    append(name, ".txt");
    if (!save(me.seqs, toCString(name))) return false;

    name = fileName;    append(name, ".rid");
    if (!save(me.names, toCString(name))) return false;

    return true;
}

// --------------------------------------------------------------------------
// Function swap()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void swap(SeqStore<TSpec, TConfig> & a, SeqStore<TSpec, TConfig> & b)
{
    std::swap(a.seqs, b.seqs);
    std::swap(a.names, b.names);
}

// ----------------------------------------------------------------------------
// Function readRecords()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void readRecords(SeqStore<TSpec, TConfig> & me,
                        FormattedFile<Fastq, Input, TFileSpec> & fileIn)
{
    readRecords(me.names, me.seqs, fileIn);
}

template <typename TSpec, typename TConfig, typename TFileSpec, typename TSize>
inline void readRecords(SeqStore<TSpec, TConfig> & me,
                        FormattedFile<Fastq, Input, TFileSpec> & fileIn,
                        TSize maxRecords)
{
    readRecords(me.names, me.seqs, fileIn, maxRecords);
}

// ----------------------------------------------------------------------------
// Function reverse()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reverse(SeqStore<TSpec, TConfig> & me)
{
    for (unsigned seqId = 0; seqId < length(me.seqs); ++seqId)
        reverse(me.seqs[seqId]);
}

// ----------------------------------------------------------------------------
// Function appendReverseComplement()
// ----------------------------------------------------------------------------
// Append reverse complemented sequences.

template <typename TSpec, typename TConfig>
void appendReverseComplement(SeqStore<TSpec, TConfig> & me)
{
    typedef SeqStore<TSpec, TConfig>    TSeqStore;
    typedef typename TSeqStore::TSeqs   TSeqs;
    typedef typename Value<TSeqs>::Type TSeq;
    typedef typename Size<TSeqs>::Type  TSeqId;

    TSeqId seqsCount = length(me.seqs);

    reserve(me.seqs, 2 * seqsCount, Exact());
    reserve(concat(me.seqs), 2 * lengthSum(me.seqs), Exact());

    for (TSeqId seqId = 0; seqId < seqsCount; ++seqId)
    {
        TSeq const & seq = me.seqs[seqId];
        appendValue(me.seqs, seq);
        reverseComplement(back(me.seqs));
    }
}
*/

#endif // #ifndef _INCLUDE_STORE_SEQS_H_
