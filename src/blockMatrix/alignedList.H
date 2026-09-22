/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::alignedList

Description
    Minimal owning list of trivially copyable T with storage aligned to
    alignedList<T>::alignment bytes (spec 6.5.3): one 4x4 blockScalar block
    or one cache line. Used for the Krylov work vectors and by
    Test-kernelBandwidth.

    It is a UList<T>, so it can be passed wherever the block library takes a
    UList (blockScalarUList&, UList<reduceScalar>&). Size changes only
    through resize_nocopy(); the contents are then unspecified. Not
    copyable; movable (List<alignedList<T>> works through default
    construction + resize_nocopy()).

    Allocation: posix_memalign. The coefficient arrays of blockLduMatrix4
    stay blockScalarList: OpenFOAM v2606 allocates every List of a
    floating-point type with at least 200 entries with 256-byte alignment
    (ListPolicy::use_alignment, default_alignment), which covers the
    64-byte requirement without changing the public matrix interface
    (Test-kernelBandwidth records the alignment it finds).

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_alignedList_H
#define coupledFoam_alignedList_H

#include "UList.H"
#include "error.H"
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
class alignedList
:
    public UList<T>
{
    static_assert
    (
        std::is_trivially_copyable_v<T>,
        "alignedList holds trivially copyable types only"
    );

    // Private Member Functions

        static T* allocate(const label n)
        {
            if (n <= 0)
            {
                return nullptr;
            }
            void* p = nullptr;
            if
            (
                ::posix_memalign
                (
                    &p,
                    alignment,
                    static_cast<std::size_t>(n)*sizeof(T)
                )
            )
            {
                FatalErrorInFunction
                    << "posix_memalign failed for " << n << " entries"
                    << abort(FatalError);
            }
            return static_cast<T*>(p);
        }

        void release() noexcept
        {
            std::free(this->data());
            UList<T>::shallowCopy(nullptr);
        }


public:

    //- Alignment in bytes (one 4x4 float block, one cache line)
    static constexpr std::size_t alignment = 64;


    // Constructors

        //- Empty
        alignedList() noexcept
        :
            UList<T>()
        {}

        //- Size n, contents unspecified
        explicit alignedList(const label n)
        :
            UList<T>(allocate(n), (n > 0 ? n : 0))
        {}

        //- Size n, zero-initialised
        alignedList(const label n, const Foam::zero)
        :
            alignedList(n)
        {
            setZero();
        }

        //- Copy of a list
        explicit alignedList(const UList<T>& list)
        :
            alignedList(list.size())
        {
            if (list.size())
            {
                std::memcpy
                (
                    static_cast<void*>(this->data()),
                    list.cdata(),
                    static_cast<std::size_t>(list.size())*sizeof(T)
                );
            }
        }

        alignedList(const alignedList&) = delete;
        void operator=(const alignedList&) = delete;

        alignedList(alignedList&& rhs) noexcept
        :
            UList<T>(rhs.data(), rhs.size())
        {
            rhs.UList<T>::shallowCopy(nullptr);
        }

        void operator=(alignedList&& rhs) noexcept
        {
            if (this != &rhs)
            {
                release();
                UList<T>::shallowCopy(rhs.data(), rhs.size());
                rhs.UList<T>::shallowCopy(nullptr);
            }
        }


    //- Destructor
    ~alignedList()
    {
        release();
    }


    // Member Functions

        //- Resize without keeping the contents
        void resize_nocopy(const label n)
        {
            if (n != this->size())
            {
                release();
                T* p = allocate(n);
                UList<T>::shallowCopy(p, (n > 0 ? n : 0));
            }
        }

        //- Set every entry to zero
        void setZero() noexcept
        {
            if (this->size())
            {
                std::memset
                (
                    static_cast<void*>(this->data()),
                    0,
                    static_cast<std::size_t>(this->size())*sizeof(T)
                );
            }
        }

        //- Set every entry to zero
        void operator=(const Foam::zero) noexcept
        {
            setZero();
        }

        //- Copy the contents of a list of the same size
        void copyFrom(const UList<T>& list)
        {
            if (list.size() != this->size())
            {
                FatalErrorInFunction
                    << "size mismatch " << list.size() << " != "
                    << this->size() << abort(FatalError);
            }
            if (this->size())
            {
                std::memcpy
                (
                    static_cast<void*>(this->data()),
                    list.cdata(),
                    static_cast<std::size_t>(list.size())*sizeof(T)
                );
            }
        }

        //- True if p is aligned to alignment bytes
        static bool isAligned(const void* p) noexcept
        {
            return reinterpret_cast<std::uintptr_t>(p) % alignment == 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
