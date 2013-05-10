//
// Generated file, do not edit! Created by opp_msgc 4.2 from Mobility.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "Mobility_m.h"

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




Register_Class(MobilityMessage);

MobilityMessage::MobilityMessage(const char *name, int kind) : cMessage(name,kind)
{
    this->eventType_var = 0;
    this->prevSpaceCell_var = 0;
    this->nextSpaceCell_var = 0;
    this->eventTime_var = 0;
    this->x_var = 0;
    this->y_var = 0;
    this->z_var = 0;
    this->vx_var = 0;
    this->vy_var = 0;
    this->vz_var = 0;
}

MobilityMessage::MobilityMessage(const MobilityMessage& other) : cMessage(other)
{
    copy(other);
}

MobilityMessage::~MobilityMessage()
{
}

MobilityMessage& MobilityMessage::operator=(const MobilityMessage& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void MobilityMessage::copy(const MobilityMessage& other)
{
    this->eventType_var = other.eventType_var;
    this->prevSpaceCell_var = other.prevSpaceCell_var;
    this->nextSpaceCell_var = other.nextSpaceCell_var;
    this->eventTime_var = other.eventTime_var;
    this->x_var = other.x_var;
    this->y_var = other.y_var;
    this->z_var = other.z_var;
    this->vx_var = other.vx_var;
    this->vy_var = other.vy_var;
    this->vz_var = other.vz_var;
    this->partner_var = other.partner_var;
}

void MobilityMessage::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->eventType_var);
    doPacking(b,this->prevSpaceCell_var);
    doPacking(b,this->nextSpaceCell_var);
    doPacking(b,this->eventTime_var);
    doPacking(b,this->x_var);
    doPacking(b,this->y_var);
    doPacking(b,this->z_var);
    doPacking(b,this->vx_var);
    doPacking(b,this->vy_var);
    doPacking(b,this->vz_var);
    doPacking(b,this->partner_var);
}

void MobilityMessage::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->eventType_var);
    doUnpacking(b,this->prevSpaceCell_var);
    doUnpacking(b,this->nextSpaceCell_var);
    doUnpacking(b,this->eventTime_var);
    doUnpacking(b,this->x_var);
    doUnpacking(b,this->y_var);
    doUnpacking(b,this->z_var);
    doUnpacking(b,this->vx_var);
    doUnpacking(b,this->vy_var);
    doUnpacking(b,this->vz_var);
    doUnpacking(b,this->partner_var);
}

int MobilityMessage::getEventType() const
{
    return eventType_var;
}

void MobilityMessage::setEventType(int eventType)
{
    this->eventType_var = eventType;
}

int MobilityMessage::getPrevSpaceCell() const
{
    return prevSpaceCell_var;
}

void MobilityMessage::setPrevSpaceCell(int prevSpaceCell)
{
    this->prevSpaceCell_var = prevSpaceCell;
}

int MobilityMessage::getNextSpaceCell() const
{
    return nextSpaceCell_var;
}

void MobilityMessage::setNextSpaceCell(int nextSpaceCell)
{
    this->nextSpaceCell_var = nextSpaceCell;
}

double MobilityMessage::getEventTime() const
{
    return eventTime_var;
}

void MobilityMessage::setEventTime(double eventTime)
{
    this->eventTime_var = eventTime;
}

double MobilityMessage::getX() const
{
    return x_var;
}

void MobilityMessage::setX(double x)
{
    this->x_var = x;
}

double MobilityMessage::getY() const
{
    return y_var;
}

void MobilityMessage::setY(double y)
{
    this->y_var = y;
}

double MobilityMessage::getZ() const
{
    return z_var;
}

void MobilityMessage::setZ(double z)
{
    this->z_var = z;
}

double MobilityMessage::getVx() const
{
    return vx_var;
}

void MobilityMessage::setVx(double vx)
{
    this->vx_var = vx;
}

double MobilityMessage::getVy() const
{
    return vy_var;
}

void MobilityMessage::setVy(double vy)
{
    this->vy_var = vy;
}

double MobilityMessage::getVz() const
{
    return vz_var;
}

void MobilityMessage::setVz(double vz)
{
    this->vz_var = vz;
}

ParticlePtr& MobilityMessage::getPartner()
{
    return partner_var;
}

void MobilityMessage::setPartner(const ParticlePtr& partner)
{
    this->partner_var = partner;
}

class MobilityMessageDescriptor : public cClassDescriptor
{
  public:
    MobilityMessageDescriptor();
    virtual ~MobilityMessageDescriptor();

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

Register_ClassDescriptor(MobilityMessageDescriptor);

MobilityMessageDescriptor::MobilityMessageDescriptor() : cClassDescriptor("MobilityMessage", "cMessage")
{
}

MobilityMessageDescriptor::~MobilityMessageDescriptor()
{
}

bool MobilityMessageDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<MobilityMessage *>(obj)!=NULL;
}

const char *MobilityMessageDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int MobilityMessageDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 11+basedesc->getFieldCount(object) : 11;
}

unsigned int MobilityMessageDescriptor::getFieldTypeFlags(void *object, int field) const
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
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISEDITABLE,
        FD_ISCOMPOUND,
    };
    return (field>=0 && field<11) ? fieldTypeFlags[field] : 0;
}

const char *MobilityMessageDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "eventType",
        "prevSpaceCell",
        "nextSpaceCell",
        "eventTime",
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "partner",
    };
    return (field>=0 && field<11) ? fieldNames[field] : NULL;
}

int MobilityMessageDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='e' && strcmp(fieldName, "eventType")==0) return base+0;
    if (fieldName[0]=='p' && strcmp(fieldName, "prevSpaceCell")==0) return base+1;
    if (fieldName[0]=='n' && strcmp(fieldName, "nextSpaceCell")==0) return base+2;
    if (fieldName[0]=='e' && strcmp(fieldName, "eventTime")==0) return base+3;
    if (fieldName[0]=='x' && strcmp(fieldName, "x")==0) return base+4;
    if (fieldName[0]=='y' && strcmp(fieldName, "y")==0) return base+5;
    if (fieldName[0]=='z' && strcmp(fieldName, "z")==0) return base+6;
    if (fieldName[0]=='v' && strcmp(fieldName, "vx")==0) return base+7;
    if (fieldName[0]=='v' && strcmp(fieldName, "vy")==0) return base+8;
    if (fieldName[0]=='v' && strcmp(fieldName, "vz")==0) return base+9;
    if (fieldName[0]=='p' && strcmp(fieldName, "partner")==0) return base+10;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *MobilityMessageDescriptor::getFieldTypeString(void *object, int field) const
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
        "int",
        "double",
        "double",
        "double",
        "double",
        "double",
        "double",
        "double",
        "ParticlePtr",
    };
    return (field>=0 && field<11) ? fieldTypeStrings[field] : NULL;
}

const char *MobilityMessageDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int MobilityMessageDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    MobilityMessage *pp = (MobilityMessage *)object; (void)pp;
    switch (field) {
        default: return 0;
    }
}

std::string MobilityMessageDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    MobilityMessage *pp = (MobilityMessage *)object; (void)pp;
    switch (field) {
        case 0: return long2string(pp->getEventType());
        case 1: return long2string(pp->getPrevSpaceCell());
        case 2: return long2string(pp->getNextSpaceCell());
        case 3: return double2string(pp->getEventTime());
        case 4: return double2string(pp->getX());
        case 5: return double2string(pp->getY());
        case 6: return double2string(pp->getZ());
        case 7: return double2string(pp->getVx());
        case 8: return double2string(pp->getVy());
        case 9: return double2string(pp->getVz());
        case 10: {std::stringstream out; out << pp->getPartner(); return out.str();}
        default: return "";
    }
}

bool MobilityMessageDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    MobilityMessage *pp = (MobilityMessage *)object; (void)pp;
    switch (field) {
        case 0: pp->setEventType(string2long(value)); return true;
        case 1: pp->setPrevSpaceCell(string2long(value)); return true;
        case 2: pp->setNextSpaceCell(string2long(value)); return true;
        case 3: pp->setEventTime(string2double(value)); return true;
        case 4: pp->setX(string2double(value)); return true;
        case 5: pp->setY(string2double(value)); return true;
        case 6: pp->setZ(string2double(value)); return true;
        case 7: pp->setVx(string2double(value)); return true;
        case 8: pp->setVy(string2double(value)); return true;
        case 9: pp->setVz(string2double(value)); return true;
        default: return false;
    }
}

const char *MobilityMessageDescriptor::getFieldStructName(void *object, int field) const
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
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        "ParticlePtr",
    };
    return (field>=0 && field<11) ? fieldStructNames[field] : NULL;
}

void *MobilityMessageDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    MobilityMessage *pp = (MobilityMessage *)object; (void)pp;
    switch (field) {
        case 10: return (void *)(&pp->getPartner()); break;
        default: return NULL;
    }
}


