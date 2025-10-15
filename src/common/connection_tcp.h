/*
 * Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
 *
 * This program is free software : you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <stdlib.h>

#ifdef _WIN32

//#      include <iostream>
//#      include <winsock2.h>
//#      include <ws2tcpip.h>
//#	include <windows.h>
#	include <winsock2.h>
#	include <mswsock.h>
#	include <ws2tcpip.h>

#	pragma comment(lib, "Ws2_32.lib")
#	pragma comment(lib, "Mswsock.lib")
#	pragma comment(lib, "AdvApi32.lib")

#else
#	include <arpa/inet.h>
#	include <netdb.h>
#	include <netinet/in.h>
#	include <netinet/tcp.h>
#	include <sys/socket.h>
#	include <unistd.h>
#endif

#define TCP_OPTIMIZATION
//#define TCP_FLOAT
#define MAX_CONNECTIONS 100


class TcpConnection {
private:
	int g_port_offset = -1;

	int g_server_id_cam[MAX_CONNECTIONS];
	int g_client_id_cam[MAX_CONNECTIONS];

	int g_server_id_data[MAX_CONNECTIONS];
	int g_client_id_data[MAX_CONNECTIONS];

	int g_timeval_sec = 60;
	int g_connection_error = 0;

	sockaddr_in g_client_sockaddr_cam[MAX_CONNECTIONS];
	sockaddr_in g_server_sockaddr_cam[MAX_CONNECTIONS];

	sockaddr_in g_client_sockaddr_data[MAX_CONNECTIONS];
	sockaddr_in g_server_sockaddr_data[MAX_CONNECTIONS];

	bool g_is_server = true;

public:
	void write_data_kernelglobal(void* data, size_t size);
	bool read_data_kernelglobal(void* data, size_t size);
	void close_kernelglobal();

	bool is_error();

	void init_sockets_cam(const char* server = NULL, int port_cam = 0, int port_data = 0, bool is_server = true);
	void init_sockets_data(const char* server = NULL, int port = 0, bool is_server = true);

	bool client_check();
	bool server_check();

	void client_close();
	void server_close();

	void send_data_cam(char* data, size_t size, bool ack = false);
	void recv_data_cam(char* data, size_t size, bool ack = false);

	void send_data_data(char* data, size_t size, bool ack = false);
	void recv_data_data(char* data, size_t size, bool ack = false);

	void set_port_offset(int offset);

private:
	int setsock_tcp_windowsize(int inSock, int inTCPWin, int inSend);
	bool init_wsa();
	void init_port();
	void close_wsa();
	bool server_create(int port,
		int& server_id,
		int& client_id,
		sockaddr_in& server_sock,
		sockaddr_in& client_sock,
		bool only_accept);

	bool client_create(const char* server_name, int port, int& client_id, sockaddr_in& client_sock);
	void close_tcp(int id);

	//void send_data(char* data, size_t size);
	//void recv_data(char* data, size_t size);
};
