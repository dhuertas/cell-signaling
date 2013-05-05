//
// Generated file, do not edit! Created by opp_msgc 4.2 from WallCollision.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "WallCollision_m.h"

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




Register_Class(WallCollisionMessage);

WallCollisionMessage::WallCollisionMessage(const char *name, int kind) : cMessage(name,kind)
{
    this->x_var = 0;
    this->y_var = 0;
    this->z_var = 0;
    this->vx_var = 0;
    this->vy_var = 0;
    this->vz_var = 0;
    this->collisionTime_var = 0;
}

WallCollisionMessage::WallCollisionMessage(const WallCollisionMessage& other) : cMessage(other)
{
    copy(other);
}

WallCollisionMessage::~WallCollisionMessage()
{
}

WallCollisionMessage& WallCollisionMessage::operator=(const WallCollisionMessage& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void WallCollisionMessage::copy(const WallCollisionMessage& other)
{
    this->x_var = other.x_var;
    this->y_var = other.y_var;
    this->z_var = other.z_var;
    this->vx_var = other.vx_var;
    this->vy_var = other.vy_var;
    this->vz_var = other.vz_var;
    this->collisionTime_var = other.collisionTime_var;
}

void WallCollisionMessage::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->x_var);
    doPacking(b,this->y_var);
    doPacking(b,this->z_var);
    doPacking(b,this->vx_var);
    doPacking(b,this->vy_var);
    doPacking(b,this->vz_var);
    doPacking(b,this->collisionTime_var);
}

void WallCollisionMessage::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->x_var);
    doUnpacking(b,this->y_var);
    doUnpacking(b,this->z_var);
    doUnpacking(b,this->vx_var);
    doUnpacking(b,this->vy_var);
    doUnpacking(b,this->vz_var);
    doUnpacking(b,this->collisionTime_var);
}

double WallCollisionMessage::getX() const
{
    return x_var;
}

void WallCollisionMessage::setX(double x)
{
    this->x_var = x;
}

double WallCollisionMessage::getY() const
{
    return y_var;
}

void WallCollisionMessage::setY(double y)
{
    this->y_var = y;
}

double WallCollisionMessage::getZ() const
{
    return z_var;
}

void WallCollisionMessage::setZ(double z)
{
    this->z_var = z;
}

double WallCollisionMessage::getVx() const
{
    return vx_var;
}

void WallCollisionMessage::setVx(double vx)
{
    this->vx_var = vx;
}

double WallCollisionMessage::getVy() const
{
    return vy_var;
}

void WallCollisionMessage::setVy(double vy)
{
    this->vy_var = vy;
}

double WallCollisionMessage::getVz() const
{
    return vz_var;
}

void WallCollisionMessage::setVz(double vz)
{
    this->vz_var = vz;
}

double WallCollisionMessage::getCollisionTime() const
{
    return collisionTime_var;
}

void WallCollisionMessage::setCollisionTime(double collisionTime)
{
    this->collisionTime_var = collisionTime;
}

class WallCollisionMessageDescriptor : public cClassDescriptor
{
  public:
    WallCollisionMessageDescriptor();
    virtual ~WallCollisionMessageDescriptor();

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

Register_ClassDescriptor(WallCollisionMessageDescriptor);

WallCollisionMessageDescriptor::WallCollisionMessageDescriptor() : cClassDescriptor("WallCollisionMessage", "cMessage")
{
}

WallCollisionMessageDescriptor::~WallCollisionMessageDescriptor()
{
}

bool WallCollisionMessageDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<WallCollisionMessage *>(obj)!=NULL;
}

const char *WallCollisionMessageDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int WallCollisionMessageDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 7+basedesc->getFieldCount(object) : 7;
}

unsigned int WallCollisionMessageDescriptor::getFieldTypeFlags(void *object, int field) const
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
    };
    return (field>=0 && field<7) ? fieldTypeFlags[field] : 0;
}

const char *WallCollisionMessageDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "collisionTime",
    };
    return (field>=0 && field<7) ? fieldNames[field] : NULL;
}

int WallCollisionMessageDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='x' && strcmp(fieldName, "x")==0) return base+0;
    if (fieldName[0]=='y' && strcmp(fieldName, "y")==0) return base+1;
    if (fieldName[0]=='z' && strcmp(fieldName, "z")==0) return base+2;
    if (fieldName[0]=='v' && strcmp(fieldName, "vx")==0) return base+3;
    if (fieldName[0]=='v' && strcmp(fieldName, "vy")==0) return base+4;
    if (fieldName[0]=='v' && strcmp(fieldName, "vz")==0) return base+5;
    if (fieldName[0]=='c' && strcmp(fieldName, "collisionTime")==0) return base+6;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *WallCollisionMessageDescriptor::getFieldTypeString(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldTypeString(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldTypeStrings[] = {
        "double",
        "double",
        "double",
        "double",
        "double",
        "double",
        "double",
    };
    return (field>=0 && field<7) ? fieldTypeStrings[field] : NULL;
}

const char *WallCollisionMessageDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int WallCollisionMessageDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    WallCollisionMessage *pp = (WallCollisionMessage *)object; (void)pp;
    switch (field) {
        default: return 0;
    }
}

std::string WallCollisionMessageDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    WallCollisionMessage *pp = (WallCollisionMessage *)object; (void)pp;
    switch (field) {
        case 0: return double2string(pp->getX());
        case 1: return double2string(pp->getY());
        case 2: return double2string(pp->getZ());
        case 3: return double2string(pp->getVx());
        case 4: return double2string(pp->getVy());
        case 5: return double2string(pp->getVz());
        case 6: return double2string(pp->getCollisionTime());
        default: return "";
    }
}

bool WallCollisionMessageDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    WallCollisionMessage *pp = (WallCollisionMessage *)object; (void)pp;
    switch (field) {
        case 0: pp->setX(string2double(value)); return true;
        case 1: pp->setY(string2double(value)); return true;
        case 2: pp->setZ(string2double(value)); return true;
        case 3: pp->setVx(string2double(value)); return true;
        case 4: pp->setVy(string2double(value)); return true;
        case 5: pp->setVz(string2double(value)); return true;
        case 6: pp->setCollisionTime(string2double(value)); return true;
        default: return false;
    }
}

const char *WallCollisionMessageDescriptor::getFieldStructName(void *object, int field) const
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
    };
    return (field>=0 && field<7) ? fieldStructNames[field] : NULL;
}

void *WallCollisionMessageDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    WallCollisionMessage *pp = (WallCollisionMessage *)object; (void)pp;
    switch (field) {
        default: return NULL;
    }
}


