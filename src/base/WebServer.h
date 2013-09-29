#ifndef WEBSERVER_H
#define WEBSERVER_H

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <algorithm>

#include "Defines.h"

struct arg_struct {
	int arg1;
	int arg2;
};

typedef struct header {
	std::string name;
	std::string value;
} header_t;

typedef struct request {
	std::string method;
	std::string uri;
	std::string version;
	std::string body;
	std::vector<header_t> headers;
} request_t;

typedef struct response {
	unsigned short status_code;
	std::string reason_phrase;
	std::string version;
	std::vector<header_t> headers;
} response_t;

typedef std::vector<std::string> vectstr_t;

class WebServer {

	private:

	protected:

		unsigned short listenPort;

		int quitFd;
		int streamFd;

		int serverSockFd;
		int clientSockFd;

		socklen_t clientSize;

		struct sockaddr_in serverAddr, clientAddr;

		std::string documentRoot;

		settings_t settings;

	public:

		WebServer(int);
		WebServer(int, unsigned short);
		WebServer(int, int, unsigned short);
		~WebServer();

		void start();
		void stop(void);

		void acceptConn();

		void handleRequest(void);
		request_t receiveRequest(void);

		void route(request_t *, response_t *);

		void addResponseHeader(response_t *, std::string, std::string);
		void sendResponseHeaders(response_t *);
		void sendResponseBody();

		void sendSettings();
		void sendStream();

		void setSettings(settings_t s) { settings = s; } ;

};

#endif