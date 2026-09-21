/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::jsonWriter

Description
    Minimal JSON object writer for measured results (spec rule 0.4: every
    numerical result goes to JSON). Values are collected in insertion order
    and written by the master rank only. Non-finite numbers are written as
    null (JSON has no NaN/Inf) and counted in the key "_nonFinite".

    Usage:
        jsonWriter j;
        j.add("nIter", 12);
        j.add("residual", 1.2e-9);
        j.add("pass", true);
        j.addList("cellsPerLevel", labelList(...));
        j.addRaw("gamg", other.str());      // nested object
        j.write("results/tests/foo.json");

SourceFiles
    jsonWriter.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_jsonWriter_H
#define coupledFoam_jsonWriter_H

#include "scalar.H"
#include "label.H"
#include "word.H"
#include "fileName.H"
#include "List.H"
#include "DynamicList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class jsonWriter
{
    // Private Data

        //- key and already-serialised value
        DynamicList<Tuple2<word, std::string>> entries_;

        //- Count of non-finite numbers written as null
        label nNonFinite_;


    // Private Member Functions

        static std::string quote(const std::string& s);

        std::string number(const doubleScalar v);


public:

    // Constructors

        jsonWriter();


    // Member Functions

        void add(const word& key, const doubleScalar v);
        void add(const word& key, const floatScalar v);
        void add(const word& key, const label v);
        void add(const word& key, const bool v);
        void add(const word& key, const std::string& v);
        void add(const word& key, const char* v);
        void add(const word& key, const word& v);

        void addList(const word& key, const UList<doubleScalar>& v);
        void addList(const word& key, const UList<label>& v);

        //- Add an already-serialised JSON value (object, array)
        void addRaw(const word& key, const std::string& json);

        //- Serialise
        std::string str() const;

        //- Write to file (master only); creates the parent directory
        void write(const fileName& file) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
