//
// Generated file, do not edit! Created by opp_msgc 4.2 from src/base/messages/Transfer.msg.
//

#ifndef _TRANSFER_M_H_
#define _TRANSFER_M_H_

#include <omnetpp.h>

// opp_msgc version check
#define MSGC_VERSION 0x0402
#if (MSGC_VERSION!=OMNETPP_VERSION)
#    error Version mismatch! Probably this file was generated by an earlier version of opp_msgc: 'make clean' should help.
#endif

// cplusplus {{
#include "../base/Manager.h"
	#include "../base/Particle.h"

	typedef Manager *ManagerPtr;
// }}



/**
 * Class generated from <tt>src/base/messages/Transfer.msg</tt> by opp_msgc.
 * <pre>
 * message TransferMessage {
 * 
 * 	int prevSpaceCell;
 * 	int nextSpaceCell;
 * 
 * 	double transferTime;
 * 
 * 	ManagerPtr manager;
 * 
 * }
 * </pre>
 */
class TransferMessage : public ::cMessage
{
  protected:
    int prevSpaceCell_var;
    int nextSpaceCell_var;
    double transferTime_var;
    ManagerPtr manager_var;

  private:
    void copy(const TransferMessage& other);

  protected:
    // protected and unimplemented operator==(), to prevent accidental usage
    bool operator==(const TransferMessage&);

  public:
    TransferMessage(const char *name=NULL, int kind=0);
    TransferMessage(const TransferMessage& other);
    virtual ~TransferMessage();
    TransferMessage& operator=(const TransferMessage& other);
    virtual TransferMessage *dup() const {return new TransferMessage(*this);}
    virtual void parsimPack(cCommBuffer *b);
    virtual void parsimUnpack(cCommBuffer *b);

    // field getter/setter methods
    virtual int getPrevSpaceCell() const;
    virtual void setPrevSpaceCell(int prevSpaceCell);
    virtual int getNextSpaceCell() const;
    virtual void setNextSpaceCell(int nextSpaceCell);
    virtual double getTransferTime() const;
    virtual void setTransferTime(double transferTime);
    virtual ManagerPtr& getManager();
    virtual const ManagerPtr& getManager() const {return const_cast<TransferMessage*>(this)->getManager();}
    virtual void setManager(const ManagerPtr& manager);
};

inline void doPacking(cCommBuffer *b, TransferMessage& obj) {obj.parsimPack(b);}
inline void doUnpacking(cCommBuffer *b, TransferMessage& obj) {obj.parsimUnpack(b);}


#endif // _TRANSFER_M_H_
