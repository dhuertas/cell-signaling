//
// Generated file, do not edit! Created by opp_msgc 4.2 from src/base/messages/OutOfNeighborhood.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "OutOfNeighborhood_m.h"

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




Register_Class(OutOfNeighborhoodMessage);

OutOfNeighborhoodMessage::OutOfNeighborhoodMessage(const char *name, int kind) : cMessage(name,kind)
{
    this->outOfNeighborhoodTime_var = 0;
}

OutOfNeighborhoodMessage::OutOfNeighborhoodMessage(const OutOfNeighborhoodMessage& other) : cMessage(other)
{
    copy(other);
}

OutOfNeighborhoodMessage::~OutOfNeighborhoodMessage()
{
}

OutOfNeighborhoodMessage& OutOfNeighborhoodMessage::operator=(const OutOfNeighborhoodMessage& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void OutOfNeighborhoodMessage::copy(const OutOfNeighborhoodMessage& other)
{
    this->outOfNeighborhoodTime_var = other.outOfNeighborhoodTime_var;
}

void OutOfNeighborhoodMessage::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->outOfNeighborhoodTime_var);
}

void OutOfNeighborhoodMessage::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->outOfNeighborhoodTime_var);
}

double OutOfNeighborhoodMessage::getOutOfNeighborhoodTime() const
{
    return outOfNeighborhoodTime_var;
}

void OutOfNeighborhoodMessage::setOutOfNeighborhoodTime(double outOfNeighborhoodTime)
{
    this->outOfNeighborhoodTime_var = outOfNeighborhoodTime;
}

class OutOfNeighborhoodMessageDescriptor : public cClassDescriptor
{
  public:
    OutOfNeighborhoodMessageDescriptor();
    virtual ~OutOfNeighborhoodMessageDescriptor();

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

Register_ClassDescriptor(OutOfNeighborhoodMessageDescriptor);

OutOfNeighborhoodMessageDescriptor::OutOfNeighborhoodMessageDescriptor() : cClassDescriptor("OutOfNeighborhoodMessage", "cMessage")
{
}

OutOfNeighborhoodMessageDescriptor::~OutOfNeighborhoodMessageDescriptor()
{
}

bool OutOfNeighborhoodMessageDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<OutOfNeighborhoodMessage *>(obj)!=NULL;
}

const char *OutOfNeighborhoodMessageDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int OutOfNeighborhoodMessageDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 1+basedesc->getFieldCount(object) : 1;
}

unsigned int OutOfNeighborhoodMessageDescriptor::getFieldTypeFlags(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeFlags(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,
    };
    return (field>=0 && field<1) ? fieldTypeFlags[field] : 0;
}

const char *OutOfNeighborhoodMessageDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "outOfNeighborhoodTime",
    };
    return (field>=0 && field<1) ? fieldNames[field] : NULL;
}

int OutOfNeighborhoodMessageDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='o' && strcmp(fieldName, "outOfNeighborhoodTime")==0) return base+0;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *OutOfNeighborhoodMessageDescriptor::getFieldTypeString(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeString(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldTypeStrings[] = {
        "double",
    };
    return (field>=0 && field<1) ? fieldTypeStrings[field] : NULL;
}

const char *OutOfNeighborhoodMessageDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int OutOfNeighborhoodMessageDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    OutOfNeighborhoodMessage *pp = (OutOfNeighborhoodMessage *)object; (void)pp;
    switch (field) {
        default: return 0;
    }
}

std::string OutOfNeighborhoodMessageDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    OutOfNeighborhoodMessage *pp = (OutOfNeighborhoodMessage *)object; (void)pp;
    switch (field) {
        case 0: return double2string(pp->getOutOfNeighborhoodTime());
        default: return "";
    }
}

bool OutOfNeighborhoodMessageDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    OutOfNeighborhoodMessage *pp = (OutOfNeighborhoodMessage *)object; (void)pp;
    switch (field) {
        case 0: pp->setOutOfNeighborhoodTime(string2double(value)); return true;
        default: return false;
    }
}

const char *OutOfNeighborhoodMessageDescriptor::getFieldStructName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldStructNames[] = {
        NULL,
    };
    return (field>=0 && field<1) ? fieldStructNames[field] : NULL;
}

void *OutOfNeighborhoodMessageDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    OutOfNeighborhoodMessage *pp = (OutOfNeighborhoodMessage *)object; (void)pp;
    switch (field) {
        default: return NULL;
    }
}


