//
// Generated file, do not edit! Created by opp_msgc 4.2 from src/base/messages/Collision.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#include <iostream>
#include <sstream>
#include "Collision_m.h"

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




Register_Class(CollisionMessage);

CollisionMessage::CollisionMessage(const char *name, int kind) : cMessage(name,kind)
{
    this->collisionTime_var = 0;
    this->x_var = 0;
    this->y_var = 0;
    this->z_var = 0;
    this->vx_var = 0;
    this->vy_var = 0;
    this->vz_var = 0;
}

CollisionMessage::CollisionMessage(const CollisionMessage& other) : cMessage(other)
{
    copy(other);
}

CollisionMessage::~CollisionMessage()
{
}

CollisionMessage& CollisionMessage::operator=(const CollisionMessage& other)
{
    if (this==&other) return *this;
    cMessage::operator=(other);
    copy(other);
    return *this;
}

void CollisionMessage::copy(const CollisionMessage& other)
{
    this->collisionTime_var = other.collisionTime_var;
    this->x_var = other.x_var;
    this->y_var = other.y_var;
    this->z_var = other.z_var;
    this->vx_var = other.vx_var;
    this->vy_var = other.vy_var;
    this->vz_var = other.vz_var;
    this->manager_var = other.manager_var;
    this->partner_var = other.partner_var;
    this->prevPartner_var = other.prevPartner_var;
}

void CollisionMessage::parsimPack(cCommBuffer *b)
{
    cMessage::parsimPack(b);
    doPacking(b,this->collisionTime_var);
    doPacking(b,this->x_var);
    doPacking(b,this->y_var);
    doPacking(b,this->z_var);
    doPacking(b,this->vx_var);
    doPacking(b,this->vy_var);
    doPacking(b,this->vz_var);
    doPacking(b,this->manager_var);
    doPacking(b,this->partner_var);
    doPacking(b,this->prevPartner_var);
}

void CollisionMessage::parsimUnpack(cCommBuffer *b)
{
    cMessage::parsimUnpack(b);
    doUnpacking(b,this->collisionTime_var);
    doUnpacking(b,this->x_var);
    doUnpacking(b,this->y_var);
    doUnpacking(b,this->z_var);
    doUnpacking(b,this->vx_var);
    doUnpacking(b,this->vy_var);
    doUnpacking(b,this->vz_var);
    doUnpacking(b,this->manager_var);
    doUnpacking(b,this->partner_var);
    doUnpacking(b,this->prevPartner_var);
}

double CollisionMessage::getCollisionTime() const
{
    return collisionTime_var;
}

void CollisionMessage::setCollisionTime(double collisionTime)
{
    this->collisionTime_var = collisionTime;
}

double CollisionMessage::getX() const
{
    return x_var;
}

void CollisionMessage::setX(double x)
{
    this->x_var = x;
}

double CollisionMessage::getY() const
{
    return y_var;
}

void CollisionMessage::setY(double y)
{
    this->y_var = y;
}

double CollisionMessage::getZ() const
{
    return z_var;
}

void CollisionMessage::setZ(double z)
{
    this->z_var = z;
}

double CollisionMessage::getVx() const
{
    return vx_var;
}

void CollisionMessage::setVx(double vx)
{
    this->vx_var = vx;
}

double CollisionMessage::getVy() const
{
    return vy_var;
}

void CollisionMessage::setVy(double vy)
{
    this->vy_var = vy;
}

double CollisionMessage::getVz() const
{
    return vz_var;
}

void CollisionMessage::setVz(double vz)
{
    this->vz_var = vz;
}

ManagerPtr& CollisionMessage::getManager()
{
    return manager_var;
}

void CollisionMessage::setManager(const ManagerPtr& manager)
{
    this->manager_var = manager;
}

ParticlePtr& CollisionMessage::getPartner()
{
    return partner_var;
}

void CollisionMessage::setPartner(const ParticlePtr& partner)
{
    this->partner_var = partner;
}

ParticlePtr& CollisionMessage::getPrevPartner()
{
    return prevPartner_var;
}

void CollisionMessage::setPrevPartner(const ParticlePtr& prevPartner)
{
    this->prevPartner_var = prevPartner;
}

class CollisionMessageDescriptor : public cClassDescriptor
{
  public:
    CollisionMessageDescriptor();
    virtual ~CollisionMessageDescriptor();

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

Register_ClassDescriptor(CollisionMessageDescriptor);

CollisionMessageDescriptor::CollisionMessageDescriptor() : cClassDescriptor("CollisionMessage", "cMessage")
{
}

CollisionMessageDescriptor::~CollisionMessageDescriptor()
{
}

bool CollisionMessageDescriptor::doesSupport(cObject *obj) const
{
    return dynamic_cast<CollisionMessage *>(obj)!=NULL;
}

const char *CollisionMessageDescriptor::getProperty(const char *propertyname) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? basedesc->getProperty(propertyname) : NULL;
}

int CollisionMessageDescriptor::getFieldCount(void *object) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    return basedesc ? 10+basedesc->getFieldCount(object) : 10;
}

unsigned int CollisionMessageDescriptor::getFieldTypeFlags(void *object, int field) const
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
        FD_ISCOMPOUND,
        FD_ISCOMPOUND,
        FD_ISCOMPOUND,
    };
    return (field>=0 && field<10) ? fieldTypeFlags[field] : 0;
}

const char *CollisionMessageDescriptor::getFieldName(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldName(object, field);
        field -= basedesc->getFieldCount(object);
    }
    static const char *fieldNames[] = {
        "collisionTime",
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "manager",
        "partner",
        "prevPartner",
    };
    return (field>=0 && field<10) ? fieldNames[field] : NULL;
}

int CollisionMessageDescriptor::findField(void *object, const char *fieldName) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    int base = basedesc ? basedesc->getFieldCount(object) : 0;
    if (fieldName[0]=='c' && strcmp(fieldName, "collisionTime")==0) return base+0;
    if (fieldName[0]=='x' && strcmp(fieldName, "x")==0) return base+1;
    if (fieldName[0]=='y' && strcmp(fieldName, "y")==0) return base+2;
    if (fieldName[0]=='z' && strcmp(fieldName, "z")==0) return base+3;
    if (fieldName[0]=='v' && strcmp(fieldName, "vx")==0) return base+4;
    if (fieldName[0]=='v' && strcmp(fieldName, "vy")==0) return base+5;
    if (fieldName[0]=='v' && strcmp(fieldName, "vz")==0) return base+6;
    if (fieldName[0]=='m' && strcmp(fieldName, "manager")==0) return base+7;
    if (fieldName[0]=='p' && strcmp(fieldName, "partner")==0) return base+8;
    if (fieldName[0]=='p' && strcmp(fieldName, "prevPartner")==0) return base+9;
    return basedesc ? basedesc->findField(object, fieldName) : -1;
}

const char *CollisionMessageDescriptor::getFieldTypeString(void *object, int field) const
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
        "ManagerPtr",
        "ParticlePtr",
        "ParticlePtr",
    };
    return (field>=0 && field<10) ? fieldTypeStrings[field] : NULL;
}

const char *CollisionMessageDescriptor::getFieldProperty(void *object, int field, const char *propertyname) const
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

int CollisionMessageDescriptor::getArraySize(void *object, int field) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getArraySize(object, field);
        field -= basedesc->getFieldCount(object);
    }
    CollisionMessage *pp = (CollisionMessage *)object; (void)pp;
    switch (field) {
        default: return 0;
    }
}

std::string CollisionMessageDescriptor::getFieldAsString(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldAsString(object,field,i);
        field -= basedesc->getFieldCount(object);
    }
    CollisionMessage *pp = (CollisionMessage *)object; (void)pp;
    switch (field) {
        case 0: return double2string(pp->getCollisionTime());
        case 1: return double2string(pp->getX());
        case 2: return double2string(pp->getY());
        case 3: return double2string(pp->getZ());
        case 4: return double2string(pp->getVx());
        case 5: return double2string(pp->getVy());
        case 6: return double2string(pp->getVz());
        case 7: {std::stringstream out; out << pp->getManager(); return out.str();}
        case 8: {std::stringstream out; out << pp->getPartner(); return out.str();}
        case 9: {std::stringstream out; out << pp->getPrevPartner(); return out.str();}
        default: return "";
    }
}

bool CollisionMessageDescriptor::setFieldAsString(void *object, int field, int i, const char *value) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->setFieldAsString(object,field,i,value);
        field -= basedesc->getFieldCount(object);
    }
    CollisionMessage *pp = (CollisionMessage *)object; (void)pp;
    switch (field) {
        case 0: pp->setCollisionTime(string2double(value)); return true;
        case 1: pp->setX(string2double(value)); return true;
        case 2: pp->setY(string2double(value)); return true;
        case 3: pp->setZ(string2double(value)); return true;
        case 4: pp->setVx(string2double(value)); return true;
        case 5: pp->setVy(string2double(value)); return true;
        case 6: pp->setVz(string2double(value)); return true;
        default: return false;
    }
}

const char *CollisionMessageDescriptor::getFieldStructName(void *object, int field) const
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
        "ManagerPtr",
        "ParticlePtr",
        "ParticlePtr",
    };
    return (field>=0 && field<10) ? fieldStructNames[field] : NULL;
}

void *CollisionMessageDescriptor::getFieldStructPointer(void *object, int field, int i) const
{
    cClassDescriptor *basedesc = getBaseClassDescriptor();
    if (basedesc) {
        if (field < basedesc->getFieldCount(object))
            return basedesc->getFieldStructPointer(object, field, i);
        field -= basedesc->getFieldCount(object);
    }
    CollisionMessage *pp = (CollisionMessage *)object; (void)pp;
    switch (field) {
        case 7: return (void *)(&pp->getManager()); break;
        case 8: return (void *)(&pp->getPartner()); break;
        case 9: return (void *)(&pp->getPrevPartner()); break;
        default: return NULL;
    }
}


