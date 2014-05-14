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

#ifndef EMITTER_H
#define EMITTER_H

#include <stdint.h>

#define EM_SQUARE		1
#define EM_GAUSSIAN		2

class Emitter {

	private:

	protected:

		bool enabled;

		uint8_t emitterStatus;

		uint64_t emitterCount;

	public:

		Emitter();

		~Emitter();

		uint64_t getEmitterCount(void) { return emitterCount; };

		uint8_t getEmitterStatus(void) { return emitterStatus; };

		void setEmitterCount(uint64_t c) { emitterCount = c; };

		void setEmitterStatus(uint8_t s) { emitterStatus = s; };

};

#endif
