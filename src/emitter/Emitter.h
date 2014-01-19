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
