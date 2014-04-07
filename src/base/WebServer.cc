#include "WebServer.h"

#include <omnetpp.h>
//#include <tcl.h>
//#include <tk.h>
//#include <tkenv.h>

// Splits a string given a delimiter and returns a vector of substrings
vectstr_t split(std::string str, char delim) {

	vectstr_t elems;
	std::stringstream ss(str);
	std::string item;

	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}

	return elems;

}

// Trim a string left and right
std::string &trim(std::string &s) {
	
	while (s.at(0) == ' ') {
		s.erase(0, 1);
	}

	while (s.at(s.length() - 1) == ' ') {
		s.erase(s.length() - 1, 1);
	}

	return s;

}

/*
 * Web Server
 */
namespace WebServer {

	// Mutex
	pthread_cond_t condSockFdFull, condSockFdEmpty;
	pthread_mutex_t mutexSockFd;

	// Threads
	pthread_t thread[MAX_THREADS];
	int threadId[MAX_THREADS];

	// Sockets
	volatile int clientSockFd[MAX_THREADS];
	volatile int clientSockFdRd;
	volatile int clientSockFdWr;
	volatile int clientSockFdCount;

	int serverSockFd, quitServerFd;

	bool quit;

	// Network
	socklen_t clientSize;

	struct sockaddr_in serverAddr;
	struct sockaddr_in clientAddr;

	typedef struct header {
		std::string name;
		std::string value;
	} header_t;

	typedef struct request {
		std::string method;
		std::string uri;
		std::string version;
		std::string resource;
		std::string query;
		std::vector<header_t> headers;
	} request_t;

	typedef struct response {
		uint16_t statusCode;
		std::string reasonPhrase;
		std::string version;
		std::string filePath;
		std::vector<header_t> headers;
	} response_t;

	uint8_t initialize();

	void requestHandler(int, int);

	void receiveRequest(int, request_t *);

	uint8_t handleRequest(int, request_t *);

	void clearRequest(request_t *);

	uint8_t handleResponse(int, request_t *, response_t *);

	void clearResponse(response_t *);

	void addResponseHeader(response_t *, std::string, std::string);

	void sendResponseHeaders(int, response_t *);

	void sendCRLF(int);

	std::string getRequestHeader(request_t *, std::string);

	void route(int, request_t *, response_t *);

	void handleGetRequest(int, request_t *, response_t *);

	void parseQueryParameters(request_t *);

	void * handler(void *);

	void run();

	uint8_t finalize();

	// Simulation related functions and parameters
	namespace Simulation {

		std::list<Particle *> *particles;

		void sendSettings(int);

		void sendStream(int);
	
		void performAction(request_t *);

		// Parameters
		uint16_t rate;

		settings_t settings;

	}

	// Constants
	std::string documentRoot; 

}

/*
 * Initializes the WebServer. Creates the mutex, conditions and threads.
 *
 * @return {uint8_t} 0 if everything went ok
 */
uint8_t WebServer::initialize() {

	int i;

	WebServer::serverSockFd = socket(PF_INET, SOCK_STREAM, 0);

	memset(&WebServer::serverAddr, 0, sizeof(struct sockaddr_in));
	memset(&WebServer::clientAddr, 0, sizeof(struct sockaddr_in));

	WebServer::serverAddr.sin_family = AF_INET;
	WebServer::serverAddr.sin_port = htons(PORT);

	WebServer::serverAddr.sin_addr.s_addr = htonl(INADDR_ANY);

	bind(WebServer::serverSockFd, (struct sockaddr *)&WebServer::serverAddr, sizeof(struct sockaddr_in));

	listen(WebServer::serverSockFd, MAX_LISTEN);

	WebServer::clientSize = sizeof(struct sockaddr_in);

	// Initialize threads mutex and conditions
	pthread_mutex_init(&WebServer::mutexSockFd, NULL);
	pthread_cond_init(&WebServer::condSockFdEmpty, NULL);
	pthread_cond_init(&WebServer::condSockFdFull, NULL);

	WebServer::clientSockFdCount = 0;
	WebServer::clientSockFdRd = 0;
	WebServer::clientSockFdWr = 0;

	WebServer::documentRoot = "./www";

	// Wake up threads
	for (i = 0; i < MAX_THREADS; i++) {

		WebServer::threadId[i] = i;

		if (pthread_create(&WebServer::thread[i], NULL, WebServer::handler, &WebServer::threadId[i]) != 0) {
			// EV << "ERROR: pthread_create" << std::endl;
		}

	}

	return 0;

}

/*
 * Clean up everything.
 *
 * @return {uint8_t} 0 if everything went ok
 */
uint8_t WebServer::finalize() {

	int i;

	for (i = 0; i < MAX_THREADS; i++) {
		pthread_join(WebServer::thread[i], NULL);
	}

	pthread_cond_destroy(&WebServer::condSockFdEmpty);
	pthread_cond_destroy(&WebServer::condSockFdFull);
	pthread_mutex_destroy(&WebServer::mutexSockFd);

	close(WebServer::serverSockFd);

	return 0;

}

/*
 * The WebServer handler. This function is called in a new thread in order
 * to handle the client request.
 *
 * @param {void *} arg: the thread id
 */
void *WebServer::handler(void *arg) {

	int *tid, sockFd;

	tid = (int *) arg;

	sigset_t sigpipe_mask;
	sigemptyset(&sigpipe_mask);
	sigaddset(&sigpipe_mask, SIGPIPE);
	sigset_t saved_mask;

	if (pthread_sigmask(SIG_BLOCK, &sigpipe_mask, &saved_mask) == -1) {
		// ev << "pthread_sigmask error" << endl;
	}

	while ( ! WebServer::quit) {

		// start of mutex area
		pthread_mutex_lock(&WebServer::mutexSockFd);

		while (WebServer::clientSockFdCount == 0) {
			pthread_cond_wait(&WebServer::condSockFdEmpty, &WebServer::mutexSockFd);
		}

		sockFd = WebServer::clientSockFd[WebServer::clientSockFdRd];

		WebServer::clientSockFdRd++;
		WebServer::clientSockFdRd = WebServer::clientSockFdRd % MAX_THREADS;
		WebServer::clientSockFdCount--;

		pthread_cond_broadcast(&WebServer::condSockFdFull);
		pthread_mutex_unlock(&WebServer::mutexSockFd);

		// end of mutex area

		WebServer::requestHandler(*tid, sockFd);

	}

	return 0;
}

/*
 * Return the request header value given a request data structure and a header
 * name.
 *
 * @param {request_t *} req
 * @param {std::string} name
 * @return {std::string}
 */
std::string WebServer::getRequestHeader(request_t *req, std::string name) {

	std::string result;

	std::vector<WebServer::header_t>::iterator it;

	for (it = req->headers.begin(); it != req->headers.end(); ++it) {

		if (it->name.compare(name) == 0) {
			result = it->value;
			break;
		}

	}

	return result;

}

/*
 * Receives the client request and handles the response. Closes the client 
 * socket after handling the response.
 *
 * @param {int} tid: the thread id that handles the request
 * @param {int} clientSockFd: the client socket
 */
void WebServer::requestHandler(int tid, int clientSockFd) {

	short n, m;

	fd_set selectSet;

	struct timeval timeout;

	WebServer::request_t req;
	WebServer::response_t resp;

	timeout.tv_sec = 1;
	timeout.tv_usec = 0;

	std::string conn;

	// Wait for data to be received or connection timeout
	FD_ZERO(&selectSet);
	FD_SET(clientSockFd, &selectSet);

	n = select(clientSockFd + 1, &selectSet, NULL, NULL, &timeout);

	if (n < 0) {

		// Client closed connection
		close(clientSockFd);
		return;

	} else if (n == 0) {

		// timeout
		close(clientSockFd);
		return;

	} else {

		if (FD_ISSET(clientSockFd, &selectSet)) {

			m = WebServer::handleRequest(clientSockFd, &req);

			if (m < 0) {

				WebServer::clearRequest(&req);
				WebServer::clearResponse(&resp);

				close(clientSockFd);
				return;

			}

		}
	
	}

	conn = WebServer::getRequestHeader(&req, "Connection");

	if (conn.length() == 0) {
		// Some http clients (e.g. curl) may not send the Connection header ...
	} else if (strncasecmp(conn.c_str(), "close", strlen("close")) == 0) {

		m = WebServer::handleResponse(clientSockFd, &req, &resp);

	} else {
		
		// Assume keep-alive connection... but keep it simple for now
		timeout.tv_sec = 0;
		timeout.tv_usec = 250;

		WebServer::handleResponse(clientSockFd, &req, &resp);

	}

	WebServer::clearRequest(&req);
	WebServer::clearResponse(&resp);

	FD_CLR(clientSockFd, &selectSet);
	FD_ZERO(&selectSet);

	close(clientSockFd);

}

/*
 * Receives the client request and fills the request_t data structure.
 *
 * @param {int} clientSockFd: the client socket
 * @param {request_t *} req: the request data structure
 * @return {uint8_t} 0 if everything went ok
 */
uint8_t WebServer::handleRequest(int clientSockFd, request_t *req) {

	vectstr_t components;

	WebServer::receiveRequest(clientSockFd, req);

	components = split(req->uri, '?');

	req->resource = components[0];

	if (components.size() == 2) {	
		req->query = components[1];	
	}

	// if Content-length header exists, get the request body part

	return 0;

}

/*
 * Reads the client request and fills the request_t data structure. This
 * function gets called from the WebServer::handleRequest function.
 *
 * @param {int} clientSockFd: the client socket
 * @param {request_t *} req: the request data structure
 */
void WebServer::receiveRequest(int clientSockFd, request_t *req) {

	char chunk;

	int received, i, n;

	std::string line;
	std::stringstream ss;

	WebServer::header_t header;
	vectstr_t splitLine;

	i = 0;
	received = 0;
	n = 0;

	ss.str("");

	// read request byte by byte
	while ((n = read(clientSockFd, &chunk, 1)) > 0) {

		ss.put(chunk);
		received += n;

		if (received >= 4 && ss.str().find("\r\n\r\n") != std::string::npos) {
			break;
		}

	}

	// parse raw data
	while (std::getline(ss, line)) {

		if (i == 0) {

			splitLine = split(line, ' ');

			req->method = trim(splitLine[0]);

			req->uri = trim(splitLine[1]);

			if (splitLine.size() > 2) {
				req->version = trim(splitLine[2]);
			}

		} else {

			splitLine = split(line, ':');

			if (splitLine.size() == 2) {

				header.name = trim(splitLine[0]);
				header.value = trim(splitLine[1]);

				req->headers.push_back(header);

			}

		}

		i++;

	}

}


/*
 * Clears the content from the given request data structure.
 *
 * @param {request_t *} req: the request data structure
 */
void WebServer::clearRequest(request_t *req) {

	req->method.clear();

	req->uri.clear();
	req->version.clear();
	req->resource.clear();
	req->query.clear();

	req->headers.clear();

}

/*
 * The main response handler. Calls the WebServer::route function.
 *
 * @param {int} clientSockFd: the client socket
 * @param {request_t *} req: the request data structure
 * @param {response_t *} resp: the response data structure
 * @return {uint8_t} 0 if everything went ok
 */
uint8_t WebServer::handleResponse(int clientSockFd, WebServer::request_t *req, WebServer::response_t *resp) {

	// get the time
	std::time_t t = std::time(NULL);
	char dateString[100];

	// prepare the response header
	resp->version = "HTTP/1.1";

	std::strftime(dateString, 100, "%c", std::localtime(&t));

	WebServer::addResponseHeader(resp, "Date", std::string(dateString));
	WebServer::addResponseHeader(resp, "Server", "Manager Server");

	// route and send the response accordingly
	route(clientSockFd, req, resp);

	return 0;

}

/*
 * Clears the content from the given response data structure.
 *
 * @param {response_t *} resp: the response data structure
 */
void WebServer::clearResponse(response_t *resp) {

	resp->statusCode = 0;

	resp->reasonPhrase.clear();
	resp->filePath.clear();

	resp->headers.clear();

}

/*
 * Routes the response in depending on the method received. Mainly used to 
 * separate GET from POST requests, though the later is not implemented yet.
 *
 * @param {int} clientSockFd: the client socket
 * @param {request_t *} req: the request data structure
 * @param {response_t *} resp: the response data structure
 */
void WebServer::route(int clientSockFd, request_t *req, response_t *resp) {

	if (strncasecmp(req->method.c_str(), "GET", strlen("GET")) == 0) {

		WebServer::handleGetRequest(clientSockFd, req, resp);

	} else {
		// method not allowed
	}

}

/*
 * Handles the GET requests (both files and/or special resources from the 
 * simulation).
 *
 * @param {int} clientSockFd: the client socket
 * @param {request_t *} req: the request data structure
 * @param {response_t *} resp: the response data structure
 */
void WebServer::handleGetRequest(int clientSockFd, request_t *req, response_t *resp) {

	int size, w;

	std::string filePath, fileStr;
	std::ifstream inFile;
	std::stringstream ss;

	ss << WebServer::documentRoot;

	if (req->resource.at(req->resource.length()-1) == '/') {
		ss << req->resource << "index.html";
	} else {
		ss << req->resource;
	}

	filePath = ss.str();
	inFile.open(filePath.c_str());

	ss.str("");

	w = 0;
	size = 0;

	// Find if client requests a file
	if (inFile) {

		resp->statusCode = 200;
		resp->reasonPhrase = "OK";

		// The file exists, and is open for input
		inFile.seekg(0, std::ifstream::end);
		size = (int) inFile.tellg();

		ss << size;

		WebServer::addResponseHeader(resp, "Length", ss.str());

		if (strncasecmp(WebServer::getRequestHeader(req, "Connection").c_str(), "close", strlen("close")) == 0) {
			WebServer::addResponseHeader(resp, "Connection", "close");
		} else {
			WebServer::addResponseHeader(resp, "Connection", "keep-alive");
		}

		// Send the response headers
		WebServer::sendResponseHeaders(clientSockFd, resp);

		inFile.seekg(0, std::ios::end);   
		fileStr.reserve(inFile.tellg());
		inFile.seekg(0, std::ios::beg);

		fileStr.assign((std::istreambuf_iterator<char>(inFile)),std::istreambuf_iterator<char>());

		if ((w = send(clientSockFd, fileStr.c_str(), fileStr.length(), 0)) != fileStr.length()) {
			// EV << "DEBUG: send" << endl;
		}

		if ((w = send(clientSockFd, "\r\n", 2, 0)) != 2) {
			// EV << "DEBUG: send" << endl;
		}

		// close file
		inFile.close();

	} else if (req->uri.find("/settings") == 0) {

		resp->statusCode = 200;
		resp->reasonPhrase = "OK";

		if (strncasecmp(WebServer::getRequestHeader(req, "Connection").c_str(), "close", strlen("close")) == 0) {
			WebServer::addResponseHeader(resp, "Connection", "close");
		} else {
			WebServer::addResponseHeader(resp, "Connection", "keep-alive");
		}

		WebServer::sendResponseHeaders(clientSockFd, resp);
		WebServer::Simulation::sendSettings(clientSockFd);

		if ((w = send(clientSockFd, "\r\n", 2, 0)) != 2) {
			// EV << "DEBUG: send" << endl;
		}

	} else if (req->uri.find("/simstream") == 0) {

		resp->statusCode = 200;
		resp->reasonPhrase = "OK";

		if (strncasecmp(WebServer::getRequestHeader(req, "Connection").c_str(), "close", strlen("close")) == 0) {
			WebServer::addResponseHeader(resp, "Connection", "close");
		} else {
			WebServer::addResponseHeader(resp, "Connection", "keep-alive");
		}

		WebServer::sendResponseHeaders(clientSockFd, resp);
		WebServer::parseQueryParameters(req);
		WebServer::Simulation::sendStream(clientSockFd);
		WebServer::sendCRLF(clientSockFd);

	} else if (req->uri.find("/action") == 0) {

		resp->statusCode = 200;
		resp->reasonPhrase = "OK";

		WebServer::addResponseHeader(resp, "Connection", "close");
		WebServer::sendResponseHeaders(clientSockFd, resp);

		WebServer::Simulation::performAction(req);

	} else {

		resp->statusCode = 404;
		resp->reasonPhrase = "Not Found";

		WebServer::addResponseHeader(resp, "Connection", "close");
		WebServer::sendResponseHeaders(clientSockFd, resp);
		WebServer::sendCRLF(clientSockFd);

	}

	ss.str("");

}

/*
 * Parse the query parameters.
 *
 * @param {request_t *} req: the request data structure
 */
void WebServer::parseQueryParameters(request_t *req) {

	vectstr_t params, tmp;
	vectstr_t::iterator p;

	if (req->uri.find("?") != std::string::npos) {

		tmp = split(req->uri, '?');
		req->query = tmp[1];

		params = split(req->query, '&');

		tmp.clear();

		for (p = params.begin(); p != params.end(); ++p) {

			tmp = split(*p, '=');

			if (tmp[0].compare("rate") == 0) {

				std::istringstream(tmp[1]) >> WebServer::Simulation::rate;

			}

// Add more if sections here ...

			tmp.clear();

		}

	}

}

/*
 * Sends the response header back to the client.
 * 
 * @param {int} clientSockFd: the client socket
 * @param {response_t *} resp: the response data structure
 */
void WebServer::sendResponseHeaders(int clientSockFd, response_t *resp) {

	unsigned int w;

	std::stringstream out;
	std::vector<header_t>::iterator hi;

	out << resp->version << " ";
	out << resp->statusCode << " ";
	out << resp->reasonPhrase << "\r\n";

	for (hi = resp->headers.begin(); hi != resp->headers.end(); ++hi) {

		out << hi->name << ": " << hi->value << "\r\n";

	}

	out << "\r\n";

	if ((w = send(clientSockFd, out.str().c_str(), out.str().size(), 0)) != out.str().size()) {
		// EV << "DEBUG: unable to send response headers" << endl;
	}

}

/*
 * Sends the settings from the ongoing simulation.
 *
 * @param {int} clientSockFd: the client socket
 */
void WebServer::Simulation::sendSettings(int clientSockFd) {

	unsigned int w;

	std::stringstream out;

	out << "{";

	out << "\"simSpaceSize\":";

	out << "{\"x\":" << WebServer::Simulation::settings.simSpaceSize.x << ",";
	out << " \"y\":" << WebServer::Simulation::settings.simSpaceSize.y << ",";
	out << " \"z\":" << WebServer::Simulation::settings.simSpaceSize.z << "},";

	out << "\"numberOfParticles\":" << WebServer::Simulation::settings.numberOfParticles;

	out << "}";

	if ((w = send(clientSockFd, out.str().c_str(), out.str().size(), 0)) != out.str().size()) {
		// EV << "DEBUG: unable to send settings" << endl;
	}

}

/*
 * Sends a stream of data from the ongoing simulation.
 *
 * @param {int} clientSockFd: the client socket 
 */
void WebServer::Simulation::sendStream(int clientSockFd) {

	bool connClosed = false;
	bool simStopped = false;

	int w, uSleepTime;
	int stoppedCount;
	int statsCount;
	uint64_t t1, t2;

	double st; // simulation time
	double dt; // delta time

	uint64_t uSimTime, uPrevSimTime;

	struct timeval t;

	std::stringstream buffer;
	std::list<Particle *>::iterator p, pt;

	uSleepTime = 1000*1000/rate;

	stoppedCount = rate;
	statsCount = rate;

	while ( ! (connClosed || WebServer::quit)) {

		st = simTime().dbl();
		uSimTime = (uint64_t )SIMTIME_RAW(simTime());

		// if simtime does not change in "rate" times, consider that the simulation is stopped
		if (uPrevSimTime == uSimTime) {

			stoppedCount--;

			if (stoppedCount == 0) {
				simStopped = true;
			}

		} else {

			stoppedCount = rate;
			simStopped = false;

		}

		gettimeofday(&t, NULL);
		t1 = t.tv_sec*1000*1000 + t.tv_usec;

		if ( ! simStopped) {

			// Clear buffer
			buffer.str("");

			// Prepare data to be sent
			buffer << "[";

			for (p = WebServer::Simulation::particles->begin();
				p != WebServer::Simulation::particles->end(); 
				++p) {

				dt = (st - (*p)->getLastCollisionTime());

				buffer << "{";

				buffer << "\"id\":" << (*p)->getParticleId() << ",";

				buffer << "\"radius\":" << (*p)->getRadius() << ",";

				buffer << "\"pos\":{";
				buffer << "\"x\":" << (*p)->getX() + (*p)->getVx()*dt << ",";
				buffer << "\"y\":" << (*p)->getY() + (*p)->getVy()*dt << ",";
				buffer << "\"z\":" << (*p)->getZ() + (*p)->getVz()*dt;
				buffer << "}";

				buffer << "}";

				pt = p;
				pt++;

				if ( ! (p != WebServer::Simulation::particles->end() && 
						pt == WebServer::Simulation::particles->end())) {

					buffer << ",";

				}

			}

			// Add a stats JSON object at the end
			if (statsCount == 0) {

				buffer << ",{";
				buffer << "\"id\": -1,";
				buffer << "\"stats\": true,";
				buffer << "\"st\":" << st;// Simulation time
				buffer << "}";

				statsCount = rate;

			} else {
				statsCount--;
			}

			buffer << "];";

		}

		// Subtract the amount of time it has taken to loop through the list so we truly
		// wait the requested time
		gettimeofday(&t, NULL);
		t2 = t.tv_sec*1000*1000 + t.tv_usec;

		if (uSleepTime - (t2-t1) > 0) {
			usleep(uSleepTime - (t2-t1));
		}

		if ( ! simStopped) {

			w = send(clientSockFd, buffer.str().c_str(), buffer.str().length(), 0);

			if (w == -1) {
				connClosed = true;
			}

		}

		uPrevSimTime = uSimTime;

	}

}

/*
 * Perform the requested action from the client. Mainly used to 
 * start and stop the ongoing simulation.
 *
 * @param {request_t *} req: the request data structure
 */
void WebServer::Simulation::performAction(request_t *req) {

	uint8_t action;

	vectstr_t params, tmp;
	vectstr_t::iterator p;

	action = 0;

	if (req->uri.find("?") != std::string::npos) {

		tmp = split(req->uri, '?');
		req->query = tmp[1];

		params = split(req->query, '&');

		tmp.clear();

		for (p = params.begin(); p != params.end(); ++p) {

			tmp = split(*p, '=');

			if (tmp[0].compare("cmd") == 0) {

				if (tmp[1].compare("start") == 0) {
					action = 1;
				}

				if (tmp[1].compare("stop") == 0) {
					action = 2;
				}

			}

			tmp.clear();

		}
	}

	switch (action) {
/*
		case 1:
			if (Tcl_Eval(getTkenv()->getInterp(), "run_normal") == TCL_ERROR) {
				EV << "Tcl_Eval error" << "\n";
			}
			break;

		case 2:
			if (Tcl_Eval(getTkenv()->getInterp(), "stop_simulation") == TCL_ERROR) {
				EV << "Tcl_Eval error" << "\n";
			}
			break;
*/
		case 0:
		default:
			break;
	}
}

/*
 * Sends a carriage return - line feed tuple to the client.
 *
 * @param {int} clientSockFd: the client socket
 */
void WebServer::sendCRLF(int clientSockFd) {

	uint8_t w;

	if ((w = send(clientSockFd, "\r\n", 2, 0)) != 2) {
		// EV << "DEBUG: send" << endl;
	}

}

/*
 * WebServer main thread. Initializes the web server, and enters the main loop
 * accepting client connections. Sends a broadcast notification to the idle 
 * threads to handle the request at the returned socket from the accept() call.
 *
 * The server loops forever until data is available at the quitServerFd file 
 * descriptor, that is pipe-connected to the Manager Module of the simulator.
 */
void WebServer::run() {

	int sockfd, maxFd;

	fd_set selectSet;

	struct timeval timeout;

	WebServer::quit = false;

	// Initialize() calls pthread_create, which in turn each thread
	// runs the handler() function. This function has a while loop using
	// the WebServer::quit variable.

	if (WebServer::initialize() != 0) {
		// EV << "ERROR: initialize" << std::endl;
	}

	timeout.tv_sec = QUIT_TIME_OUT;
	timeout.tv_usec = 0;

	while ( ! WebServer::quit) {

		sockfd = accept(WebServer::serverSockFd, (struct sockaddr *)&WebServer::clientAddr, &WebServer::clientSize);

		FD_ZERO(&selectSet);
		FD_SET(sockfd, &selectSet);
		FD_SET(WebServer::quitServerFd, &selectSet);

		maxFd = std::max(sockfd, WebServer::quitServerFd);

		select(maxFd + 1, &selectSet, NULL, NULL, &timeout);

		if (FD_ISSET(sockfd, &selectSet)) {

			// start of mutex area
			pthread_mutex_lock(&WebServer::mutexSockFd);

			while (WebServer::clientSockFdCount > MAX_THREADS) {
					pthread_cond_wait(&WebServer::condSockFdFull, &WebServer::mutexSockFd);
			}

			WebServer::clientSockFd[WebServer::clientSockFdWr] = sockfd;

			WebServer::clientSockFdCount++;
			WebServer::clientSockFdWr++;
			WebServer::clientSockFdWr = WebServer::clientSockFdWr % MAX_THREADS;

			pthread_cond_broadcast(&WebServer::condSockFdEmpty);
			pthread_mutex_unlock(&WebServer::mutexSockFd);

			// end of mutex area

		}

		if (FD_ISSET(WebServer::quitServerFd, &selectSet)) {
			WebServer::quit = true;
		}

		FD_CLR(sockfd, &selectSet);
		FD_CLR(WebServer::quitServerFd, &selectSet);

	}

	if (WebServer::finalize() != 0) {
		// EV << "ERROR: finalize" << std::endl;
	}

}

/*
 * Adds a header (key:value pair) to the response data structure.
 *
 * @param {response_t *} resp
 * @param {header_t} header
 */
void WebServer::addResponseHeader(response_t *resp, std::string name, std::string value) {

	header_t header;

	header.name = name;
	header.value = value;

	resp->headers.push_back(header);

}

/*
 * Starts a web server thread. The Manager module spawns a new thread during 
 * the initialization process.
 *
 * @param {void *} arguments: the thread arguments
 */
void *startServerThread(void *arguments) {

	struct arg_struct *args = (struct arg_struct *)arguments;

	WebServer::quitServerFd = (int )args->quitFd;

	WebServer::Simulation::settings = (settings_t )(args->settings);
	WebServer::Simulation::particles = args->particles;

	WebServer::run();

	pthread_exit(NULL);

}

/*
 * Stops the server thread. This functions is called from the Manager module
 * during the finalization process.
 *
 * @param {int} quitFd: the write end of the pipe used to stop the server
 */
void endServerThread(int quitFd) {

	if (write(quitFd, "1", 1) != 1) {
	    EV << "write error\n";
	}

}
