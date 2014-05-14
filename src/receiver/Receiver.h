//  Omnet++ project to simulate cell signaling communications 
//  Copyright (C) 2014  Daniel Huertas
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RECEIVER_H
#define RECEIVER_H

#include <stdint.h>

class Receiver {

	private:

	protected:

		bool enabled;

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