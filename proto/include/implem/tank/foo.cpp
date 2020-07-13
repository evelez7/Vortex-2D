typedef LBD  LevelBoxData<T,NCOMPS,MEMTYPE,CENTERING>;
#define LBDT template<typename T, unsigned int NCOMPS=1,MemType MEMTYPE=HOST, unsigned int CENTERING=DIM>

LBDT
void LBD::allocateBuffers(const LBD& a_src,
                          const Interval& a_srcComps,
                          const LBD& a_dest,
                          const Interval& a_destComps,
                          const Copier&   a_copier,
                          const LBD& a_op) const
{
  CH_TIME("MPI_allocateBuffers");
  m_buff = &(((Copier&)a_copier).m_buffers);
  a_dest.m_buff = m_buff;
  
  CH_assert(a_srcComps.size() == a_destComps.size());
  if (m_buff->isDefined(a_srcComps.size())) return;

  m_buff->m_ncomps = a_srcComps.size();

  m_buff->m_fromMe.resize(0);
  m_buff->m_toMe.resize(0);
  size_t sendBufferSize = 0;
  size_t recBufferSize  = 0;

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
    m_buff->m_recbuffer = malloc(recBufferSize);
    if (m_buff->m_recbuffer == NULL)
    {
      MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
    }
    m_buff->m_reccapacity = recBufferSize;
  }

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
