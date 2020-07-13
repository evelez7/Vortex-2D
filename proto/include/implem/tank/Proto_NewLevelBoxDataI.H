#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Chombo_LEVELBOXDATAI_H_
#define _Chombo_LEVELBOXDATAI_H_

#include <cstdlib>
#include <algorithm>
#include <limits.h>
#include <list>
#include "Chombo_CH_OpenMP.H"
#include "Chombo_parstream.H"
#include "Chombo_CH_Timer.H"
#include "Chombo_NamespaceHeader.H"


using std::sort;


LBDT
inline bool LBD::isDefined() const
{
  return m_isdefined;
}

LBDT
inline void LBD::setVector(const LBD& da,
                           const Interval& srcComps,
                           const Interval& destComps)
{
  if(&da != this)
  {
    DataIterator it=this->dataIterator();
    int nbox=it.size();
    // Replace with DataIterator.
    for(int box=0; box<nbox; box++)
    {
//      this->m_vector[it[box].datInd()]->copy( this->box(it[box]), destComps,
//                                              this->box(it[box]), da[it[box]], srcComps);
      int isrc =  srcComps.begin();
      int idst = destComps.begin();
      int inco = destComps.size();
      ::Proto::Box bx = ProtoCh::getProtoBox(this->box(it[box]));
      this->m_vector[it[box].datInd()]->copy(da[it[box]], bx, isrc, bx, idst, inco);

    }
  }
}
LBDT
void LDB::define(const DisjointBoxLayout& a_dbl, 
                 const IntVect& a_ghost)
{
  //CH_TIME("LevelData<T>::define(dbl,comps,ghost,factory)");
  // clear exchange copier if it's already been defined
  if (this->m_isdefined)
  {
    m_exchangeCopier.clear();
  }

  this->m_isdefined = true;
  if (!a_dbl.isClosed())
  {
    MayDay::Error("non-disjoint DisjointBoxLayout: LevelData<T>::define(const DisjointBoxLayout& dp,....)");
  }
  if (comps<=0)
  {
    MayDay::Error("LevelData::LevelData(const BoxLayout& dp, int comps)  comps<=0");
  }
  m_dbl = a_dbl;
  m_ghost = a_ghost;
  //no need for exchange copier if no ghost
  //do not do this if it is periodic because bad things can happen
  //with large number of ghost cells.
  if(m_ghost != IntVect::Zero)
  {
    m_exchangeCopier.define(m_dbl, m_dbl, m_ghost, true);
  }
  allocateGhostVector(ghost);
}
//-----------------------------------------------------------------------

LDBT
inline void LDB::clear()
{
  if (this->m_callDelete == true)
  {
    for (unsigned int i=0; i<this->m_vector.size(); ++i)
    {
      delete this->m_vector[i];
      this->m_vector[i] = NULL;
    }
  }
  m_isdefined = false;
}

LDBT
inline void LDB::allocateGhostVector()
{
  if (this->m_callDelete == true)
  {
    for (unsigned int i=0; i<this->m_vector.size(); ++i)
    {
      delete this->m_vector[i];
      this->m_vector[i] = NULL;
    }
  }

  this->m_callDelete = factory.callDelete();

  DataIterator it(this->dataIterator()); int nbox=it.size();
  this->m_vector.resize(it.size(), NULL);
#pragma omp parallel for if(this->m_threadSafe)
  for(int i=0; i<nbox; i++)
  {
    unsigned int index = it[i].datInd();
    Box abox = this->box(it[i]);
    abox.grow(m_ghost);
    this->m_vector[index] = factory.create(abox, m_comps, it[i]);
    if (this->m_vector[index] == NULL)
    {
      MayDay::Error("OutOfMemory in BoxLayoutData::allocateGhostVector");
    }
  }
}

//======================================================================

template<class T>
void BoxLayoutData<T>::makeItSo(const Interval&   a_srcComps,
                                const BoxLayoutData<T>& a_src,
                                BoxLayoutData<T>& a_dest,
                                const Interval&   a_destComps,
                                const Copier&     a_copier,
                                const LDOperator<T>& a_op) const
{
  if(s_verbosity > 0) 
  {
    pout() << "makeit so copier = " << endl;
    a_copier.print();
  }
  {
    CH_TIME("makeItSoBegin");
    makeItSoBegin(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op);
  }
  {
    CH_TIME("makeItSoLocalCopy");
    makeItSoLocalCopy(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op);
  }
  {
    CH_TIME("makeItSoEnd");
    a_dest.makeItSoEnd(a_destComps, a_op);
  }
}

template<class T>
void BoxLayoutData<T>::makeItSoBegin(const Interval&   a_srcComps,
                                     const BoxLayoutData<T>& a_src,
                                     BoxLayoutData<T>& a_dest,
                                     const Interval&   a_destComps,
                                     const Copier&     a_copier,
                                     const LDOperator<T>& a_op) const
{
  // The following five functions are nullOps in uniprocessor mode

#ifdef CH_MPI

  allocateBuffers(a_src,  a_srcComps,
                  a_dest, a_destComps,
                  a_copier,
                  a_op);  //monkey with buffers, set up 'fromMe' and 'toMe' queues

  writeSendDataFromMeIntoBuffers(a_src, a_srcComps, a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef CH_MPI is for the m_buff->m_toMe and m_buff->m_fromMe
  {
    CH_TIME("post_messages");
    this->m_buff->numReceives = m_buff->m_toMe.size();

    if (this->m_buff->numReceives > 0)
    {
      postReceivesToMe(); // all non-blocking
    }
  

    this->m_buff->numSends = m_buff->m_fromMe.size();
    if (this->m_buff->numSends > 0)
    {
      postSendsFromMe();  // all non-blocking
    }
  }    
#endif 
}

template<class T>
void BoxLayoutData<T>::makeItSoLocalCopy(const Interval&   a_srcComps,
                                         const BoxLayoutData<T>& a_src,
                                         BoxLayoutData<T>& a_dest,
                                         const Interval&   a_destComps,
                                         const Copier&     a_copier,
                                         const LDOperator<T>& a_op) const
{

  CH_TIME("local copying");
  CopyIterator it(a_copier, CopyIterator::LOCAL);  
  int items=it.size();
#ifdef _OPENMP
  bool threadSafe = m_threadSafe && (a_op.threadSafe());
#endif
#pragma omp parallel for if(threadSafe)
  for (int n=0; n<items; n++)
  {
    const otionItem& item = it[n];
//debugging bit to force serial code to run parallel bits
#if 0
    const T & srcFAB = a_src[item.fromIndex];
    T       & dstFAB = a_dest[item.toIndex];
    size_t bufsize_src = a_op.size(srcFAB, item.fromRegion, a_srcComps);
    size_t bufsize_dst = a_op.size(srcFAB, item.fromRegion, a_srcComps);
    if(bufsize_src != bufsize_dst)
    {
      MayDay::Error("buffer size mismatch");
    }
    char* charbuffer = new char[bufsize_src];
    a_op.linearOut(srcFAB, charbuffer, item.fromRegion, a_srcComps);
    a_op.linearIn (dstFAB, charbuffer, item.toRegion,  a_destComps);

    delete[] charbuffer;
#else
    a_op.op(a_dest[item.toIndex], item.fromRegion,
            a_destComps,
            item.toRegion,
            a_src[item.fromIndex],
            a_srcComps);
#endif

  }
}
template<class T>
void BoxLayoutData<T>::makeItSoEnd(
  const Interval&   a_destComps,
  const LDOperator<T>& a_op)
{
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  completePendingSends(); // wait for sends from possible previous operation

  unpackReceivesToMe(a_destComps, a_op); // nullOp in uniprocessor mode

}

#ifndef CH_MPI
// uniprocessor version of all these nullop functions.
template<class T>
void BoxLayoutData<T>::completePendingSends() const
{
}

template<class T>
void BoxLayoutData<T>::allocateBuffers(const BoxLayoutData<T>& a_src,
                                       const Interval& a_srcComps,
                                       const BoxLayoutData<T>& a_dest,
                                       const Interval& a_destComps,
                                       const Copier&   a_copier,
                                       const LDOperator<T>& a_op
  ) const
{
}

template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                                      const Interval&     a_srcComps,
                                                      const LDOperator<T>& a_op) const
{
}

template<class T>
void BoxLayoutData<T>::postSendsFromMe() const
{
}

template<class T>
void BoxLayoutData<T>::postReceivesToMe() const
{
}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe(
  const Interval&   a_destComps,
  const LDOperator<T>& a_op)
{
}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                 const Interval&   a_destComps,                   
                                                 const LDOperator<T>& a_op) const
{
}

#else

// MPI versions of the above codes.

template<class T>
void BoxLayoutData<T>::completePendingSends() const
{
  CH_TIME("completePendingSends");
  if (this->m_buff->numSends > 0)
  {
    CH_TIME("MPI_Waitall");
    m_buff->m_sendStatus.resize(this->m_buff->numSends);
    int result = MPI_Waitall(this->m_buff->numSends, &(m_buff->m_sendRequests[0]), &(m_buff->m_sendStatus[0]));
    if (result != MPI_SUCCESS)
    {
      //hell if I know what to do about failed messaging here
    }
  }
  this->m_buff->numSends = 0;
}

template<class T>
void BoxLayoutData<T>::allocateBuffers(const BoxLayoutData<T>& a_src,
                                       const Interval& a_srcComps,
                                       const BoxLayoutData<T>& a_dest,
                                       const Interval& a_destComps,
                                       const Copier&   a_copier,
                                       const LDOperator<T>& a_op) const
{
  CH_TIME("MPI_allocateBuffers");
  m_buff = &(((Copier&)a_copier).m_buffers);
  a_dest.m_buff = m_buff;
  
  CH_assert(a_srcComps.size() == a_destComps.size());
  if (m_buff->isDefined(a_srcComps.size()) && T::preAllocatable()<2) return;

  if(s_verbosity > 0)
  {
    pout() << " allocate buffers srcComps = " << a_srcComps << ", dest comps = " << a_destComps << endl;
  }
  m_buff->m_ncomps = a_srcComps.size();

  m_buff->m_fromMe.resize(0);
  m_buff->m_toMe.resize(0);
  size_t sendBufferSize = 0;
  size_t recBufferSize  = 0;
  // two versions of code here.  one for preAllocatable T, one not.

  T dummy;
  for (CopyIterator it(a_copier, CopyIterator::FROM); it.ok(); ++it)
  {
    const MotionItem& item = it();
    CopierBuffer::bufEntry b;
    b.item = &item;
    b.size = a_op.size(a_src[item.fromIndex], item.fromRegion, a_srcComps);
    sendBufferSize+=b.size;
    b.procID = item.procID;
    m_buff->m_fromMe.push_back(b);
  }
  sort(m_buff->m_fromMe.begin(), m_buff->m_fromMe.end());
  for (CopyIterator it(a_copier, CopyIterator::TO); it.ok(); ++it)
  {
    const MotionItem& item = it();
    CopierBuffer::bufEntry b;
    b.item = &item;
    if (T::preAllocatable() == 0)
    {
      b.size = a_op.size(dummy, item.fromRegion, a_destComps);
      recBufferSize+=b.size;
    }
    else if (T::preAllocatable() == 1)
    {
      b.size = a_op.size(a_dest[item.toIndex], item.fromRegion, a_destComps);
      recBufferSize+=b.size;
    }
    b.procID = item.procID;
    m_buff->m_toMe.push_back(b);
  }
  sort(m_buff->m_toMe.begin(), m_buff->m_toMe.end());

  // allocate send and receveive buffer space.

  if (sendBufferSize > m_buff->m_sendcapacity)
  {
    free((m_buff->m_sendbuffer));
    if (s_verbosity > 0) pout()<<"malloc send buffer "<<sendBufferSize<<std::endl;
    (m_buff->m_sendbuffer) = malloc(sendBufferSize);
    if ((m_buff->m_sendbuffer) == NULL)
    {
      MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
    }
    m_buff->m_sendcapacity = sendBufferSize;
  }

  if (recBufferSize > m_buff->m_reccapacity)
  {
    free(m_buff->m_recbuffer);
    if (s_verbosity > 0) pout()<<"malloc receive buffer "<<recBufferSize<<std::endl;
    m_buff->m_recbuffer = malloc(recBufferSize);
    if (m_buff->m_recbuffer == NULL)
    {
      MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
    }
    m_buff->m_reccapacity = recBufferSize;
  }

  /*
    pout()<<"\n";
    for (int i=0; i<m_buff->m_fromMe.size(); i++)
    pout()<<m_buff->m_fromMe[i].item->region<<"{"<<m_buff->m_fromMe[i].procID<<"}"<<" ";
    pout() <<"::::";
    for (int i=0; i<m_buff->m_toMe.size(); i++)
    pout()<<m_buff->m_toMe[i].item->region<<"{"<<m_buff->m_toMe[i].procID<<"}"<<" ";
    pout() << endl;
  */

  char* nextFree = (char*)(m_buff->m_sendbuffer);
  if (m_buff->m_fromMe.size() > 0)
  {
    for (unsigned int i=0; i<m_buff->m_fromMe.size(); ++i)
    {
      m_buff->m_fromMe[i].bufPtr = nextFree;
      nextFree += m_buff->m_fromMe[i].size;
    }
  }

  nextFree = (char*)m_buff->m_recbuffer;
  if (m_buff->m_toMe.size() > 0)
  {
    for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
    {
      m_buff->m_toMe[i].bufPtr = nextFree;
      nextFree += m_buff->m_toMe[i].size;
    }
  }

  // since fromMe and toMe are sorted based on procID, messages can now be grouped
  // together on a per-processor basis.

}

template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                                      const Interval&     a_srcComps,
                                                      const LDOperator<T>& a_op) const
{
  CH_TIME("write Data to buffers");
  int isize = m_buff->m_fromMe.size();
#ifdef _OPENMP
  bool threadSafe = m_threadSafe && (a_op.threadSafe());
#endif
#pragma omp parallel for if(threadSafe)
  for (unsigned int i=0; i< isize; ++i)
  {
    const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
    a_op.linearOut(a_src[entry.item->fromIndex], entry.bufPtr,
                   entry.item->fromRegion, a_srcComps);
  }
}

template<class T>
void BoxLayoutData<T>::postSendsFromMe() const
{
  CH_TIME("post_Sends");
  // now we get the magic of message coalescence
  // fromMe has already been sorted in the allocateBuffers() step.

  this->m_buff->numSends = m_buff->m_fromMe.size();

  if (this->m_buff->numSends > 1)
  {
    for (unsigned int i=m_buff->m_fromMe.size()-1; i>0; --i)
    {
      if (m_buff->m_fromMe[i].procID == m_buff->m_fromMe[i-1].procID)
      {
        this->m_buff->numSends--;
        m_buff->m_fromMe[i-1].size = m_buff->m_fromMe[i-1].size + m_buff->m_fromMe[i].size;
        m_buff->m_fromMe[i].size = 0;
      }
    }
  }
  m_buff->m_sendRequests.resize(this->m_buff->numSends);
  std::list<MPI_Request> extraRequests;

  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->m_buff->numSends; ++i)
  {
    const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[next];
    char*  buffer = (char*)entry.bufPtr;
    std::size_t bsize = entry.size;
    int idtag=0;
    while (bsize > CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE)
    {
      extraRequests.push_back(MPI_Request());
      {
        //CH_TIME("MPI_Isend");
        MPI_Isend(buffer, CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                  idtag, CH4_SPMD::Chombo_MPI::comm, &(extraRequests.back()));
      }
      maxSize = CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      bsize -= CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      buffer+=CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      idtag++;
    }
    {
      //CH_TIME("MPI_Isend");
      MPI_Isend(buffer, bsize, MPI_BYTE, entry.procID,
                idtag, CH4_SPMD::Chombo_MPI::comm, &(m_buff->m_sendRequests[i]));
    }
    maxSize = std::max<long long>(bsize, maxSize);
    ++next;
    while (next < m_buff->m_fromMe.size() && m_buff->m_fromMe[next].size == 0) ++next;
  }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
  {
    m_buff->m_sendRequests.push_back(*it);
  }
  this->m_buff->numSends = m_buff->m_sendRequests.size();

  CH4_SPMD::CH_MaxMPISendSize = std::max<long long>(CH4_SPMD::CH_MaxMPISendSize, maxSize);

}

template<class T>
void BoxLayoutData<T>::postReceivesToMe() const
{
  CH_TIME("post_Receives");
  this->m_buff->numReceives = m_buff->m_toMe.size();

  if (this->m_buff->numReceives > 1)
  {
    for (unsigned int i=m_buff->m_toMe.size()-1; i>0; --i)
    {
      if (m_buff->m_toMe[i].procID == m_buff->m_toMe[i-1].procID)
      {
        this->m_buff->numReceives--;
        m_buff->m_toMe[i-1].size += m_buff->m_toMe[i].size;
        m_buff->m_toMe[i].size = 0;
      }

    }
  }
  m_buff->m_receiveRequests.resize(this->m_buff->numReceives);
  std::list<MPI_Request> extraRequests;
  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->m_buff->numReceives; ++i)
  {
    const CopierBuffer::bufEntry& entry = m_buff->m_toMe[next];
    char*  buffer = (char*)entry.bufPtr;
    size_t bsize = entry.size;
    int idtag=0;
    while (bsize > CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE)
    {
      extraRequests.push_back(MPI_Request());
      {
        //CH_TIME("MPI_Irecv");
        MPI_Irecv(buffer, CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                  idtag, CH4_SPMD::Chombo_MPI::comm, &(extraRequests.back()));
      }
      maxSize = CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      bsize -= CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      buffer+=CH4_SPMD::CH_MAX_MPI_MESSAGE_SIZE;
      idtag++;
    }
    {
      //CH_TIME("MPI_Irecv");
      MPI_Irecv(buffer, bsize, MPI_BYTE, entry.procID,
                idtag, CH4_SPMD::Chombo_MPI::comm, &(m_buff->m_receiveRequests[i]));
    }
    ++next;
    maxSize = std::max<long long>(bsize, maxSize);
    while (next < m_buff->m_toMe.size() && m_buff->m_toMe[next].size == 0) ++next;
  }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
  {
    m_buff->m_receiveRequests.push_back(*it);
  }
  this->m_buff->numReceives = m_buff->m_receiveRequests.size();

  CH4_SPMD::CH_MaxMPIRecvSize = std::max<long long>(CH4_SPMD::CH_MaxMPIRecvSize, maxSize);
  //pout()<<"maxSize="<<maxSize<<" posted "<<this->m_buff->numReceives<<" receives\n";

}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe(
  const Interval&   a_destComps,
  const LDOperator<T>& a_op)
{

  CH_TIME("unpack_messages");
  
  if (this->m_buff->numReceives > 0)
  {
    m_buff->m_receiveStatus.resize(this->m_buff->numReceives);
    int result;
    {
      CH_TIME("MPI_Waitall");
      result = MPI_Waitall(this->m_buff->numReceives, &(m_buff->m_receiveRequests[0]),
                           &(m_buff->m_receiveStatus[0]));
    }
    if (result != MPI_SUCCESS)
    {
      //hell if I know what to do about failed messaging here
      //maybe a mayday::warning?
    }

    int isize = m_buff->m_toMe.size();
#ifdef _OPENMP
    bool threadSafe = m_threadSafe && (a_op.threadSafe());
#endif
#pragma omp parallel for if(threadSafe)
    for (unsigned int i=0; i< isize; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
      a_op.linearIn(this->operator[](entry.item->toIndex), entry.bufPtr,
                    entry.item->toRegion, a_destComps);
    }
  }
  this->m_buff->numReceives = 0;
}

LDBT
void LDB::exchange()
{
  CH_TIME("exchange");
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier exchange");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  // PC::need to construct comps, op.
  this->makeItSo(comps, *this, *this, comps, m_exchangeCopier, op);
}

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------

/*template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          BoxLayoutData<T>& dest,
                          const Interval& destComps,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const*/
LDBT
LDB::copyTo(LDB& a_dest, 
            const Copier a_copier)
{
 PR_TIME("copyTo");
 // PC::Need to construct other arguments to makeItSo.
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier copyTo");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  this->makeItSo(srcComps, *this, dest, destComps, a_copier, op);
}

#include "Chombo_NamespaceFooter.H"
#endif
