#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PROTO_LEVELBOXDATA_H_
#define _PROTO_LEVELBOXDATA_H_
#define LBDT template <typename T, unsigned int C=1,MemType MEMTYPE=HOST, unsigned int CENTERING=DIM>
Enum MemType {HOST,DEVICE};
typedef LBD LevelBoxData<T,C,MEMTYPE,CENTERING>;
//--Forward declarations

template <class T,unsigned int C=1>
class LDOperator
{
public:
  // to keep the compiler happy
  /// PC:: evaluate whether this can go, as well. Test case: reuse in implementing LevelFluxRegister ?
  // PC:: Will at least reimplement the member functions so that they apply only to BoxData.
  virtual ~LDOperator()
  {
  }
  virtual int size(BoxData<T,C>& arg, const Box& b) const
  {
//    return arg.size(b, comps);
    return arg.charsize(b);
  }
  virtual void linearOut(BoxData<T,C>& arg, void* buf, const Box& R,
                         const Interval& comps) const
  {
//    arg.linearOut(buf, R, comps);
    unsigned int start = comps.begin();
    unsigned int ncomp = comps.size();
    arg.linearOut(buf, R, start, ncomp);
  }
  virtual void linearIn(BoxData<T,C>& arg,  void* buf, const Box& R,
                        const Interval& comps)const
  {
//    arg.linearIn(buf, R, comps);
    unsigned int start = comps.begin();
    unsigned int ncomp = comps.size();
    arg.linearIn(buf, R, start, ncomp);
  }

  virtual void op(BoxData<T,C>& dest,
                  const Box& a_bxFr,
                  const Interval& a_Cdest,
                  const Box& a_bxTo,
                  const BoxData<T,C>& a_src,
                  const Interval& a_Csrc) const
  {
    if (LinearizationTest == 0)
    {

      int isrc =  a_Csrc.begin();
      int idst = a_Cdest.begin();
      int inco = a_Csrc.size();
//      dest.copy(a_bxFr, Cdest,a_bxTo, src, Csrc);
      dest.copy(src, a_bxFr, isrc, a_bxTo, idst, inco);
    }
    else
    {
      int sizeSource = size(src,  RegionFrom, Csrc);
      int sizeDest   = size(dest, RegionTo,   Cdest);

      if (T::preAllocatable() == 0 || T::preAllocatable() == 1)
      {
        if (sizeSource != sizeDest)
        {
          MayDay::Abort("LinearizationTest failure: dest and source have different sizes");
        }
      }
      Vector<char> buffer(sizeSource);
      void* b = (void*)&(buffer[0]);
      linearOut(src, b, RegionFrom, Csrc);
      linearIn(dest, b, RegionTo,   Cdest);
    }
  }
};


LBDT
class LBD 
{
public:
  LBD();
  ~LBD();
  LBD(const DisjointBoxLayout& a_dbl);
  void define(const DisjointBoxLayout& a_dbl);
protected: 
  bool            m_isdefined;
  DisjointBoxLayout m_dbl;
  // BVS comment:
  // this used to be std::vector<T>, and I could let vector handle
  // destruction for me, but vector.resize() absolutely demands that
  // I provide a copy constructor for T.  People get uncomfortable when
  // I make them provide me with copy constructors.
  
  Vector<BoxData<T,C>* > m_vector;

  // BVS comment:
  // thinking about making this class self-documenting to a greater degree
  // and having the component names also be kept in the class and carried
  // around through various transformations.....
  //  vector<string> m_componentNames;

  bool m_callDelete = true;
  
  // PC:: Reimplement using BoxData copyTo.
  // PC:: what does this do ?
  //void setVector(const LBD& da,
  //              const Interval& srcComps,
  //             const Interval& destComps);

  // PC:: What does this do? I think that it allocates ghosted data for each box.
  void allocateGhostVector();

  void makeItSo(const Interval&     a_srcComps,
                const LBD&          a_src,
                LBD&                a_dest,
                const Interval&     a_destComps,
                const Copier&       a_copier,
                const LDOP&         a_op =LDOP()) const;
  
  void makeItSoBegin(
                     const Interval&     a_srcComps,
                     const LBD&          a_src,
                     LBD&                a_dest,
                     const Interval&     a_destComps,
                     const Copier&       a_copier,
                     const LDOP&         a_op = LDOP()) const; 
  void makeItSoLocalCopy(
                         const Interval&     a_srcComps,
                         const LBD&          a_src,
                         LBD&                a_dest,
                         const Interval&     a_destComps,
                         const Copier&       a_copier,
                         const LDOP&         a_op = LDOP()) const;
  void makeItSoEnd(
                   const Interval&     a_destComps,
                   const LDOP&         a_op = LDOP());
  
  //========================================================================
  //
  // data structures used by makeItSo when we have some
  // data that needs to be moved (ie. there are entries
  // in the 'FROM' or 'TO' CopyIterators)
  //
  void completePendingSends() const;

  void allocateBuffers(const LBD& a_src,
                       const Interval& a_srcComps,
                       const LBD& a_dest,
                       const Interval& a_destComps,
                       const Copier&   a_copier,
                       const LDOP& a_op) const;

  void writeSendDataFromMeIntoBuffers(const LBD& a_src,
                                      const Interval& a_srcComps,
                                      const LDOP& a_op) const;

  void postSendsFromMe() const ;

  void postReceivesToMe() const ;

  void unpackReceivesToMe(
    const Interval&   a_destComps,
    const LDOP& a_op) ;


  mutable CopierBuffer*  m_buff;
private:
  // disallow copy constructors and assignment operators
  // to avoid v
ery hard-to-find performance problems
 LBD& operator= (const LBD& a_rhs);
 LBD(const LBD& a_rhs);

  // handy usage function, tied to the implementation of this class
  // as a vector<T*>.  Assumes that m_comps and m_boxLayout have already
  // been initialized correctly.  Sets the correct size for the data
  // vector, deletes the extra T objects if m_vector is to shorten, and
  // performs either construction or define on the remaining T objects.
  void allocate();
};

#include "Proto_NamespaceFooter.H"
#include "Proto_NewLevelBoxDataI.H"

#endif //NEWLEVELBOXDATA
