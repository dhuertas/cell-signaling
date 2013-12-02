#ifndef RECEIVER_H
#define RECEIVER_H

#include <stdint.h>

#include <csimplemodule.h>
#include <cmessage.h>

#include "../Molecule.h"
#include "../base/Defines.h"

class Receiver {

	private:

	protected:

		uint8_t receiverStatus;

		uint8_t receiverType;

		uint64_t receiverCount;

	public:

		Receiver();

		~Receiver();

		uint64_t getReceiverCount(void) { return receiverCount; };

		uint8_t getReceiverStatus(void) { return receiverStatus; };

		uint8_t getReceiverType(void) { return receiverType; };

		void setReceiverCount(uint64_t c) { receiverCount = c; };

		void setReceiverStatus(uint8_t s) { receiverStatus = s; };

		void setReceiverType(uint8_t t) { receiverType = t; };

};

#endif