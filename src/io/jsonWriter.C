/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "jsonWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "UPstream.H"
#include <cmath>
#include <cstdio>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::jsonWriter::jsonWriter()
:
    entries_(),
    nNonFinite_(0)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

std::string Foam::jsonWriter::quote(const std::string& s)
{
    std::string q("\"");
    for (const char c : s)
    {
        switch (c)
        {
            case '"':  q += "\\\""; break;
            case '\\': q += "\\\\"; break;
            case '\n': q += "\\n"; break;
            case '\t': q += "\\t"; break;
            default:   q += c; break;
        }
    }
    q += '"';
    return q;
}


std::string Foam::jsonWriter::number(const doubleScalar v)
{
    if (!std::isfinite(v))
    {
        ++nNonFinite_;
        return "null";
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return buf;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::jsonWriter::add(const word& key, const doubleScalar v)
{
    entries_.append(Tuple2<word, std::string>(key, number(v)));
}


void Foam::jsonWriter::add(const word& key, const floatScalar v)
{
    add(key, static_cast<doubleScalar>(v));
}


void Foam::jsonWriter::add(const word& key, const label v)
{
    entries_.append(Tuple2<word, std::string>(key, std::to_string(v)));
}


void Foam::jsonWriter::add(const word& key, const bool v)
{
    entries_.append
    (
        Tuple2<word, std::string>(key, std::string(v ? "true" : "false"))
    );
}


void Foam::jsonWriter::add(const word& key, const std::string& v)
{
    entries_.append(Tuple2<word, std::string>(key, quote(v)));
}


void Foam::jsonWriter::add(const word& key, const char* v)
{
    add(key, std::string(v));
}


void Foam::jsonWriter::add(const word& key, const word& v)
{
    add(key, static_cast<const std::string&>(v));
}


void Foam::jsonWriter::addList(const word& key, const UList<doubleScalar>& v)
{
    std::string s("[");
    forAll(v, i)
    {
        if (i) s += ", ";
        s += number(v[i]);
    }
    s += "]";
    entries_.append(Tuple2<word, std::string>(key, s));
}


void Foam::jsonWriter::addList(const word& key, const UList<label>& v)
{
    std::string s("[");
    forAll(v, i)
    {
        if (i) s += ", ";
        s += std::to_string(v[i]);
    }
    s += "]";
    entries_.append(Tuple2<word, std::string>(key, s));
}


void Foam::jsonWriter::addRaw(const word& key, const std::string& json)
{
    entries_.append(Tuple2<word, std::string>(key, json));
}


std::string Foam::jsonWriter::str() const
{
    std::string s("{\n");
    forAll(entries_, i)
    {
        s += "  " + quote(entries_[i].first()) + ": " + entries_[i].second();
        s += ",\n";
    }
    s += "  \"_nonFinite\": " + std::to_string(nNonFinite_) + "\n}\n";
    return s;
}


void Foam::jsonWriter::write(const fileName& file) const
{
    if (!UPstream::master())
    {
        return;
    }
    if (!file.path().empty())
    {
        Foam::mkDir(file.path());
    }
    OFstream os(file);
    os.writeQuoted(str(), false);
}


// ************************************************************************* //
