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

#ifndef __WEBSERVER_H
#define __WEBSERVER_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <ctime>
#include <list>

// integers
#include <stdint.h>

// connection
#include <sys/select.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

// threading
#include <pthread.h>

// file
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

// others
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <signal.h>

#include "Defines.h"
#include "Particle.h"

#define MAX_LISTEN			30
#define MAX_THREADS			10
#define MAX_DATE_SIZE		64
#define MAX_BUFFER			1024
#define MAX_KEEPALIVE		100

#define TRUE 				1
#define FALSE				0
#define ERROR				-1

#define PORT				8080

#define GET					0
#define HEAD				1
#define POST				2
#define PUT					3
#define DELETE				4
#define TRACE				5
#define CONNECT				6

#define SOCK_RD				0
#define SOCK_WR				1
#define TIME_OUT			10 // Keep Alive timeout 10 seconds
#define RECV_TIME_OUT		10
#define KEEPALIVE_TIME_OUT	1
#define QUIT_TIME_OUT		10

typedef std::vector<std::string> vectstr_t;

typedef struct settings {
	int numberOfParticles;
	vect_t simSpaceSize;
} settings_t;

struct arg_struct {
	int quitFd;							// Quit file descriptor
	settings_t settings;				// Simulation settings structure
	std::list<Particle *> *particles;	// List of particles
};

// utils
vectstr_t split(std::string, char);

// Trim a string left and right
std::string &trim(std::string &);

// Start server
void *startServerThread(void *);

// Stop server
void endServerThread(int);

#endif
