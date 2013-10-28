#include "WebServer.h"

#include <omnetpp.h>

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

		Tkenv *tkenv;

		void sendSettings(int);

		void sendStream(int);
	
		void performAction();

		// Parameters
		uint16_t rate;

		// Actions are:
		// 0 - None
		// 1 - Run
		// 2 - Fast Run
		// 3 - Pause
		// 4 - Stop
		uint16_t action;

		uint16_t lastAction;

		settings_t settings;

	}

	// Constants
	std::string documentRoot; 

}

/*
 *
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
 *
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
 *
 */
void *WebServer::handler(void *arg) {

	int *tid, sockFd;

	tid = (int *) arg;

	while (1) {

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
 *
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
 *
 */
void WebServer::requestHandler(int tid, int clientSockFd) {

	short n, m, reqCount;

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
	reqCount = 1;

	if (conn.length() == 0) {
		// Some http clients (e.g. curl) may not send the Connection header ...
	} else if (strncasecmp(conn.c_str(), "close", strlen("close")) == 0) {

		m = WebServer::handleResponse(clientSockFd, &req, &resp);

	} else {
		
		// Assume keep-alive connection
		timeout.tv_sec = 0;
		timeout.tv_usec = 250;

		WebServer::handleResponse(clientSockFd, &req, &resp);

		while ((n = select(clientSockFd + 1, &selectSet, NULL, NULL, &timeout)) > 0) {

			if (FD_ISSET(clientSockFd, &selectSet)) {

				WebServer::clearRequest(&req);
				WebServer::clearResponse(&resp);

				m = WebServer::handleRequest(clientSockFd, &req);
				
				if (m == 0) {

					WebServer::handleResponse(clientSockFd, &req, &resp);
					reqCount++;

				} else {

					break;

				}

			}

			if (reqCount > MAX_KEEPALIVE) {
				break;
			}

		}

	}

	WebServer::clearRequest(&req);
	WebServer::clearResponse(&resp);

	FD_CLR(clientSockFd, &selectSet);
	FD_ZERO(&selectSet);

	close(clientSockFd);

}

/*
 *
 */
uint8_t WebServer::handleRequest(int clientSockFd, WebServer::request_t *req) {

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
 *
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
 *
 */
void WebServer::clearRequest(WebServer::request_t *req) {

	req->method.clear();

	req->uri.clear();
	req->version.clear();
	req->resource.clear();
	req->query.clear();

	req->headers.clear();

}

/*
 *
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
 *
 */
void WebServer::clearResponse(WebServer::response_t *resp) {

	resp->statusCode = 0;

	resp->reasonPhrase.clear();
	resp->filePath.clear();

	resp->headers.clear();

}

/*
 *
 */
void WebServer::route(int clientSockFd, WebServer::request_t *req, WebServer::response_t *resp) {

	if (strncasecmp(req->method.c_str(), "GET", strlen("GET")) == 0) {

		WebServer::handleGetRequest(clientSockFd, req, resp);

	} else {
		// method not allowed
	}

}

void WebServer::handleGetRequest(int clientSockFd, WebServer::request_t *req, WebServer::response_t *resp) {

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

	} else if (req->uri.find("/action") == 0) {

		resp->statusCode = 200;
		resp->reasonPhrase = "OK";

		if (strncasecmp(WebServer::getRequestHeader(req, "Connection").c_str(), "close", strlen("close")) == 0) {
			WebServer::addResponseHeader(resp, "Connection", "close");
		} else {
			WebServer::addResponseHeader(resp, "Connection", "keep-alive");
		}

		WebServer::parseQueryParameters(req);
		WebServer::Simulation::performAction();
		WebServer::sendResponseHeaders(clientSockFd, resp);
		WebServer::sendCRLF(clientSockFd);

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
 * Parse the query parameters
 *
 * @param {request_t *} req
 */
void WebServer::parseQueryParameters(WebServer::request_t *req) {

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

			if (tmp[0].compare("cmd") == 0) {

				if (tmp[1].compare("start") == 0) {

					WebServer::Simulation::action = 1;

				} else if (tmp[1].compare("fast") == 0) {

					WebServer::Simulation::action = 2;

				} else if (tmp[1].compare("stop") == 0) {

					WebServer::Simulation::action = 3;

				} else {

					WebServer::Simulation::action = 0;

				}

			}
// Add more if sections here ...

			tmp.clear();

		}

	}

}

/*
 *
 */
void WebServer::sendResponseHeaders(int clientSockFd, WebServer::response_t *resp) {

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
 *
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
 *
 */
void WebServer::Simulation::sendStream(int clientSockFd) {

	bool connClosed = false;
	bool simStopped = false;

	int w, uSleepTime;
	unsigned int t1, t2;

	double sTime, dTime;

	struct timeval t;

	std::stringstream buffer;
	std::list<Particle *>::iterator p, pt;

	uSleepTime = 1000*1000/rate;

	simStopped = WebServer::Simulation::tkenv->getStopSimulationFlag();

	while ( ! connClosed) {

		// Clear buffer
		buffer.str("");

		if ( ! simStopped) {

			// Compute the amount of time it takes to generate the JSON data
			gettimeofday(&t, NULL);
			t1 = t.tv_sec*1000*1000 + t.tv_usec;

			// Prepare data to be sent
			buffer << "[";

			sTime = simTime().dbl();

			for (p = WebServer::Simulation::particles->begin();
				p != WebServer::Simulation::particles->end(); 
				++p) {

				dTime = (sTime - (*p)->getLastCollisionTime());

				buffer << "{";

				buffer << "\"id\":" << (*p)->getParticleId() << ",";

				buffer << "\"radius\":" << (*p)->getRadius() << ",";

				buffer << "\"pos\":{";
				buffer << "\"x\":" << (*p)->getX() + (*p)->getVx()*dTime << ",";
				buffer << "\"y\":" << (*p)->getY() + (*p)->getVy()*dTime << ",";
				buffer << "\"z\":" << (*p)->getZ() + (*p)->getVz()*dTime;
				buffer << "}";

				buffer << "}";

				pt = p;
				pt++;

				if ( ! (p != WebServer::Simulation::particles->end() &&
					pt == WebServer::Simulation::particles->end())) {

					buffer << ",";

				}

			}

			buffer << "];";
			// Subtract the amount of time it has taken to loop through the list so we truly
			// wait the requested time
			gettimeofday(&t, NULL);
			t2 = t.tv_sec*1000*1000 + t.tv_usec;

			if (uSleepTime - (t2-t1) > 0) {
				usleep(uSleepTime -t2+t1);
			}

			w = send(clientSockFd, buffer.str().c_str(), buffer.str().length(), 0);

			if (w == -1) connClosed = true;

		}

		simStopped = WebServer::Simulation::tkenv->getStopSimulationFlag();

	}

}

/*
 *
 */
void WebServer::Simulation::performAction() {

	switch (WebServer::Simulation::action) {

		case 1:

			//WebServer::Simulation::tkenv->setSimulationRunMode(Tkenv::RUNMODE_NORMAL);
			//Tcl_GlobalEval(WebServer::Simulation::tkenv->getInterp(), "update");

			if (Tcl_Eval(WebServer::Simulation::tkenv->getInterp(), "opp_run") == TCL_ERROR)
            	ev << "Tkenv: " << Tcl_GetStringResult(WebServer::Simulation::tkenv->getInterp());

			break;

		case 2:

			WebServer::Simulation::tkenv->setSimulationRunMode(Tkenv::RUNMODE_FAST);

			break;

		case 3:

			WebServer::Simulation::tkenv->setStopSimulationFlag();

			break;

		case 0:
		default:
			break;
	}

	WebServer::Simulation::lastAction = WebServer::Simulation::action;
	WebServer::Simulation::action = 0;

}

/*
 *
 */
void WebServer::sendCRLF(int clientSockFd) {

    uint8_t w;

	if ((w = send(clientSockFd, "\r\n", 2, 0)) != 2) {
		// EV << "DEBUG: send" << endl;
	}

}

/*
 * WebServer main thread
 */
void WebServer::run() {

	int sockfd, maxFd;

	fd_set selectSet;

	struct timeval timeout;

	if (WebServer::initialize() != 0) {
		// EV << "ERROR: initialize" << std::endl;
	}

	WebServer::quit = false;

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
 * Adds a header (key:value pair) to the response data structure
 *
 * @param {response_t*} resp
 * @param {header_t} header
 */
void WebServer::addResponseHeader(response_t *resp, std::string name, std::string value) {

	header_t header;

	header.name = name;
	header.value = value;

	resp->headers.push_back(header);

}

// Starts the server thread
void *startServerThread(void *arguments) {

	struct arg_struct *args = (struct arg_struct *)arguments;

	WebServer::quitServerFd = (int )args->quitFd;

	WebServer::Simulation::settings = (settings_t )(args->settings);
	WebServer::Simulation::particles = args->particles;
	WebServer::Simulation::tkenv = args->tkenv;

	WebServer::run();

	pthread_exit(NULL);

}

// Stops the server thread
void endServerThread(int quitFd) {

	write(quitFd, "1", 1);

}
