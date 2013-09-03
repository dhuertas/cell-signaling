//
// Generated file, do not edit! Created by opp_msgc 4.2 from src/base/messages/Transfer.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "Transfer_m.h"

// Template rule which fires if a struct or class doesn't have operator<<
template<typename T>
std::ostream& operator<<(std::ostream& out,const T&) {return out;}

// Another default rule (prevents compiler from choosing base class' doPacking())
template<typename T>
void doPacking(cCommBuffer *, T& t) {
    throw cRuntimeError("Parsim error: no doPacking() function for type %s or its base class (check .msg and _m.cc/h files!)",opp_typename(typeid(t)));
}

template<typename T>
void doUnpacking(cCommBuffer *, T& t) {
    throw cRuntimeError("Parsim error: no doUnpacking() function for type %s or its base class (check .msg and _m.cc/h files!)",opp_typename(typeid(t)));
}




Register_Class(TransferMessage);

TransferMessage::TransferMessage(const char *name, int kind) : cMessage(name,kind)
{
    this->prevSpaceCell_var = 0;
    this->nextSpaceCell_var = 0;
    this->transferTime_var = 0;
}

TransferMessage::TransferMessage(const TransferMessage& other) : cMessage(other)
{
    copy(other);
}

TransferMessage::~TransferMessage()
{
}

TransferMessage& TransferMessage::operator=(const TransferMessage& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void TransferMessage::copy(const TransferMessage& other)
{
    this->prevSpaceCell_var = other.prevSpaceCell_var;
    this->nextSpaceCell_var = other.nextSpaceCell_var;
    this->transferTime_var = other.transferTime_var;
    this->manager_var = other.manager_var;
}

void TransferMessage::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->prevSpaceCell_var);
    doPacking(b,this->nextSpaceCell_var);
    doPacking(b,this->transferTime_var);
    doPacking(b,this->manager_var);
}

void TransferMessage::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->prevSpaceCell_var);
    doUnpacking(b,this->nextSpaceCell_var);
    doUnpacking(b,this->transferTime_var);
    doUnpacking(b,this->manager_var);
}

int TransferMessage::getPrevSpaceCell() const
{
    return prevSpaceCell_var;
}

void TransferMessage::setPrevSpaceCell(int prevSpaceCell)
{
    this->prevSpaceCell_var = prevSpaceCell;
}

int TransferMessage::getNextSpaceCell() const
{
    return nextSpaceCell_var;
}

void TransferMessage::setNextSpaceCell(int nextSpaceCell)
{
    this->nextSpaceCell_var = nextSpaceCell;
}

double TransferMessage::getTransferTime() const
{
    return transferTime_var;
}

void TransferMessage::setTransferTime(double transferTime)
{
    this->transferTime_var = transferTime;
}

ManagerPtr& TransferMessage::getManager()
{
    return manager_var;
}

void TransferMessage::setManager(const ManagerPtr& manager)
{
    this->manager_var = manager;
}

class TransferMessageDescriptor : public cClassDescriptor
{
  public:
    TransferMessageDescriptor();
    virtual ~TransferMessageDescriptor();

    virtual bool doesSupport(cObject *obj) const;
    virtual const char *getProperty(const char *propertyname) const;
    virtual int getFieldCount(void *object) const;
    virtual const char *getFieldName(void *object, int field) const;
    virtual int findField(void *object, const char *fieldName) const;
    virtual unsigned int getFieldTypeFlags(void *object, int field) const;
    virtual const char *getFieldTypeString(void *object, int field) const;
    virtual const char *getFieldProperty(void *object, int field, const char *propertyname) const;
    virtual int getArraySize(void *object, int field) const;

    virtual std::string getFieldAsString(void *object, int field, int i) const;
    virtual bool setFieldAsString(void *object, int field, int i, const char *value) const;

    virtual const char *getFieldStructName(void *object, int field) const;
    virtual void *getFieldStructPointer(void *object, int field, int i) const;
};

Register_ClassDescriptor(TransferMessageDescriptor);

TransferMessageDescriptor::TransferMessageDescriptor() : cClassDescriptor("TransferMessage", "cMessage")
{
}

TransferMessageDescriptor::~TransferMessageDescriptor()
{
}

bool TransferMessageDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<TransferMessage *>(obj)!=NULL;
}

const char *TransferMessageDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int TransferMessageDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 4+basedesc->getFieldCount(object) : 4;
}

unsigned int TransferMessageDescriptor::getFieldTypeFlags(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeFlags(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISCOMPOUND,
    };
    return (field>=0 && field<4) ? fieldTypeFlags[field] : 0;
}

const char *TransferMessageDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "prevSpaceCell",
        "nextSpaceCell",
        "transferTime",
        "manager",
    };
    return (field>=0 && field<4) ? fieldNames[field] : NULL;
}

int TransferMessageDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='p' && strcmp(fieldName, "prevSpaceCell")==0) return base+0;
    if (fieldName[0]=='n' && strcmp(fieldName, "nextSpaceCell")==0) return base+1;
    if (fieldName[0]=='t' && strcmp(fieldName, "transferTime")==0) return base+2;
    if (fieldName[0]=='m' && strcmp(fieldName, "manager")==0) return base+3;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *TransferMessageDescriptor::getFieldTypeString(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeString(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldTypeStrings[] = {
        "int",
        "int",
        "double",
        "ManagerPtr",
    };
    return (field>=0 && field<4) ? fieldTypeStrings[field] : NULL;
}

const char *TransferMessageDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldProperty(object, field, propertyname);
        field -= basedesc->getFieldCount(object);
    }
    switch (field) {
        default: return NULL;
    }
}

int TransferMessageDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    TransferMessage *pp = (TransferMessage *)object; (void)pp;
    switch (field) {
        default: return 0;
    }
}

std::string TransferMessageDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    TransferMessage *pp = (TransferMessage *)object; (void)pp;
    switch (field) {
        case 0: return long2string(pp->getPrevSpaceCell());
        case 1: return long2string(pp->getNextSpaceCell());
        case 2: return double2string(pp->getTransferTime());
        case 3: {std::stringstream out; out << pp->getManager(); return out.str();}
        default: return "";
    }
}

bool TransferMessageDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    TransferMessage *pp = (TransferMessage *)object; (void)pp;
    switch (field) {
        case 0: pp->setPrevSpaceCell(string2long(value)); return true;
        case 1: pp->setNextSpaceCell(string2long(value)); return true;
        case 2: pp->setTransferTime(string2double(value)); return true;
        default: return false;
    }
}

const char *TransferMessageDescriptor::getFieldStructName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldStructNames[] = {
        NULL,
        NULL,
        NULL,
        "ManagerPtr",
    };
    return (field>=0 && field<4) ? fieldStructNames[field] : NULL;
}

void *TransferMessageDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    TransferMessage *pp = (TransferMessage *)object; (void)pp;
    switch (field) {
        case 3: return (void *)(&pp->getManager()); break;
        default: return NULL;
    }
}


