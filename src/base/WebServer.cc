#include "WebServer.h"

using namespace std;

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

// Starts the server thread
void *serverThread(void *arguments) {

	int quitFd, streamFd;

	struct arg_struct *args = (struct arg_struct *)arguments;

	quitFd = (int )args->arg1;
	streamFd = (int )args->arg2;

	// Web server thread
	WebServer ws(quitFd, streamFd, LISTEN_PORT);

	ws.start();

	ws.stop();

	pthread_exit(NULL);

}

// Stops the server thread
void endServerThread(int quitFd) {

	write(quitFd, "1", 1);

}

/*
 *
 */
WebServer::WebServer(int fd) {

	memset(&serverAddr, 0, sizeof(serverAddr));
	memset(&clientAddr, 0, sizeof(clientAddr));

	listenPort = 80;

	quitFd = fd;

	documentRoot = ".";

}

/*
 *
 */
WebServer::WebServer(int quitFd, unsigned short port) {

	memset(&serverAddr, 0, sizeof(serverAddr));
	memset(&clientAddr, 0, sizeof(clientAddr));

	listenPort = port;

	this->quitFd = quitFd;

	documentRoot = ".";

}

WebServer::WebServer(int quitFd, int streamFd, unsigned short port) {

	memset(&serverAddr, 0, sizeof(serverAddr));
	memset(&clientAddr, 0, sizeof(clientAddr));

	listenPort = port;

	this->quitFd = quitFd;
	this->streamFd = streamFd;

	documentRoot = ".";

}

/*
 *
 */
WebServer::~WebServer() {

	memset(&serverAddr, 0, sizeof(serverAddr));
	memset(&clientAddr, 0, sizeof(clientAddr));

}

/*
 *
 */
void WebServer::addResponseHeader(response_t *resp, std::string name, std::string value) {

	header_t header;

	header.name = name;
	header.value = value;

	resp->headers.push_back(header);

}

/*
 *
 */
void WebServer::start() {

	bool quit = false;

	int in;

	fd_set selectSet;
 	struct timeval timeout;

	timeout.tv_sec = 10;
	timeout.tv_usec = 0;

	serverSockFd = socket(PF_INET, SOCK_STREAM, 0);

	serverAddr.sin_family = AF_INET;
	serverAddr.sin_port = htons(listenPort);
	serverAddr.sin_addr.s_addr = htonl(INADDR_ANY);

	bind(serverSockFd, (struct sockaddr *) &serverAddr, sizeof(serverAddr));

	listen(serverSockFd, MAX_LISTEN);

	clientSize = sizeof(clientAddr);

	// Start the main loop
	while ( ! quit) {

		FD_ZERO(&selectSet);
		FD_SET(serverSockFd, &selectSet);
		FD_SET(quitFd, &selectSet);

		in = select(max(serverSockFd, quitFd) + 1, &selectSet, NULL, NULL, &timeout);

		if (in < 0) {
			// error
		} else if (in == 0) {
			// timeout
		} else {

			if (FD_ISSET(serverSockFd, &selectSet)) {
				acceptConn();
			}

			if (FD_ISSET(quitFd, &selectSet)) {
				quit = true;
			}

		}

	}

}

/*
 *
 */
void WebServer::acceptConn() {

	int in;

	fd_set selectSet;
 	struct timeval timeout;

	timeout.tv_sec = 1;
	timeout.tv_usec = 0;

	clientSockFd = accept(serverSockFd, (struct sockaddr *) &clientAddr, &clientSize);

	FD_ZERO(&selectSet);
	FD_SET(clientSockFd, &selectSet);

	in = select(clientSockFd + 1, &selectSet, NULL, NULL, &timeout);

	if (in < 0) {
		// error
	} else if (in == 0) {
		// time out
	} else {

		if (FD_ISSET(clientSockFd, &selectSet)) {
			handleRequest();
		}

	}

	close(clientSockFd);

}

/*
 *
 */
void WebServer::stop() {

	close(serverSockFd);

}

/*
 *
 */
void WebServer::handleRequest() {

	header_t header;
	request_t req;
	response_t resp;

	// get the time
	std::time_t t = std::time(NULL);
	char date_string[100];

	// get the request
	req = receiveRequest();

	// prepare the response header
	resp.status_code = 200;
	resp.reason_phrase = "OK";
	resp.version = "HTTP/1.1";

	std::strftime(date_string, 100, "%c", std::localtime(&t));

	addResponseHeader(&resp, "Date", std::string(date_string));
	addResponseHeader(&resp, "Server", "Manager Server");

	// route and send the response accordingly
	route(&req, &resp);

}

/*
 *
 */
request_t WebServer::receiveRequest() {

	char chunk;

	int received, i, n;

	std::string line;
	std::stringstream ss;

	header_t header;
	request_t req;
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
	while(std::getline(ss, line)) {

		if (i == 0) {

			splitLine = split(line, ' ');

			req.method = trim(splitLine[0]);
			req.uri = trim(splitLine[1]);

			if (splitLine.size() > 2) {
				req.version = trim(splitLine[2]);
			}

		} else {

			splitLine = split(line, ':');

			if (splitLine.size() == 2) {

				header.name = trim(splitLine[0]);
				header.value = trim(splitLine[1]);

				req.headers.push_back(header);

			}

		}

		i++;

	}

	return req;

}

/*
 *
 */
void WebServer::route(request_t *req, response_t *resp) {

	int size, w;

	std::string filePath, fileStr;
	std::ifstream inFile;
	std::stringstream ss;

	vectstr_t components, resourcePath;

	components = split(req->uri, '?');

	ss << documentRoot;

	if (components[0].at(components[0].length()-1) == '/') {

		ss << components[0] << "index.html";

	} else {

		ss << components[0];

	}

	filePath = ss.str();
	inFile.open(filePath.c_str());

	ss.str("");

	w = 0;
	size = 0;

	cout << req->uri << endl;

	// Find if client requests a file
	if (inFile) {

		// The file exists, and is open for input
		inFile.seekg(0, std::ifstream::end);
		size = (int) inFile.tellg();

		ss << size;

		addResponseHeader(resp, "Connection", "close");
		addResponseHeader(resp, "Length", ss.str());

		// Send the response headers
		sendResponseHeaders(resp);

		inFile.seekg(0, std::ios::end);   
		fileStr.reserve(inFile.tellg());
		inFile.seekg(0, std::ios::beg);

		fileStr.assign((std::istreambuf_iterator<char>(inFile)),
			std::istreambuf_iterator<char>());

		if ((w = send(clientSockFd, fileStr.c_str(), fileStr.length(), 0)) != fileStr.length()) {
			cout << "DEBUG: send" << endl;
		}

		if ((w = send(clientSockFd, "\r\n", 2, 0)) != 2) {
			cout << "DEBUG: send" << endl;
		}

		// close file
		inFile.close();

	} else if (req->uri.compare("/settings") == 0) {
		// send settings as a JSON object
		sendSettings();

	} else if (req->uri.compare("/simstream") == 0) {
		// perform requested action
		addResponseHeader(resp, "Connection", "keep-alive");

		sendResponseHeaders(resp);

		// Read query options if necessary

		// Send data stream
		sendStream();

	} else {

		resp->status_code = 404;
		resp->reason_phrase = "Not Found";

		addResponseHeader(resp, "Connection", "close");

		sendResponseHeaders(resp);

		if ((w = send(clientSockFd, "\r\n", 2, 0)) != 2) {
			cout << "DEBUG: send" << endl;
		}

	}

	ss.str("");

}

/*
 *
 */
void WebServer::sendResponseHeaders(response_t *resp) {

	unsigned int w;

	std::stringstream out;
	std::vector<header_t>::iterator hi;

	out << resp->version << " ";
	out << resp->status_code << " ";
	out << resp->reason_phrase << "\r\n";

	for (hi = resp->headers.begin(); hi != resp->headers.end(); ++hi) {

		out << hi->name << ": " << hi->value << "\r\n";

	}

	out << "\r\n";

	if ((w = send(clientSockFd, out.str().c_str(), out.str().size(), 0)) != out.str().size()) {
		cout << "DEBUG: unable to send response headers" << endl;
	}

}

/*
 *
 */
void WebServer::sendSettings() {
	
}

/*
 *
 */
void WebServer::sendStream() {

	bool quit = false;

	char buffer[1024];

	int in, r, w;

	fd_set selectSet;
 	struct timeval timeout;

	timeout.tv_sec = 10;
	timeout.tv_usec = 0;

	while ( ! quit) {

		FD_ZERO(&selectSet);
		FD_SET(quitFd, &selectSet);
		FD_SET(streamFd, &selectSet);

		in = select(max(quitFd, streamFd) + 1, &selectSet, NULL, NULL, &timeout);

		if (in < 0) {
			// error
		} else if (in == 0) {
			// timeout
		} else {

			if (FD_ISSET(quitFd, &selectSet)) {
				quit = true;
			}

			if (FD_ISSET(streamFd, &selectSet)) {

				r = read(streamFd, buffer, 1024);

				if (r > 0) {

					w = send(clientSockFd, buffer, r, 0);

					if (w != r) {
						// error
					}

				} else if (r == 0) {

				} else {
					//error
				}

				r = 0;
				memset(buffer, 0, 1024);

			}

		}

	}

}