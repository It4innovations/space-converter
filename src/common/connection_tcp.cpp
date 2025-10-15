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

#include "connection_tcp.h"

#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>


//#include <omp.h>

// RGB
#  define TCP_WIN_SIZE_SEND (32L * 1024L * 1024L)
#  define TCP_WIN_SIZE_RECV (32L * 1024L * 1024L)

#  define TCP_BLK_SIZE (1L * 1024L * 1024L * 1024L
#  define TCP_MAX_SIZE (128L * 1024L * 1024L)


#ifdef _WIN32
#  define KERNEL_SOCKET_SEND(s, buf, len) send(s, buf, (int)len, 0)
#  define KERNEL_SOCKET_RECV(s, buf, len) recv(s, buf, (int)len, 0)
#  define KERNEL_SOCKET_SEND_IGNORE_RC(s, buf, len) send(s, buf, (int)len, 0)
#  define KERNEL_SOCKET_RECV_IGNORE_RC(s, buf, len) recv(s, buf, (int)len, 0)
#  define socklen_t int
#else
#  define KERNEL_SOCKET_SEND(s, buf, len) write(s, buf, len); 
#  define KERNEL_SOCKET_RECV(s, buf, len) read(s, buf, len); 
#  define KERNEL_SOCKET_SEND_IGNORE_RC(s, buf, len) { auto rc = write(s, buf, len); assert(rc == len); }
#  define KERNEL_SOCKET_RECV_IGNORE_RC(s, buf, len) { auto rc = read(s, buf, len); assert(rc == len); }
#endif

int TcpConnection::setsock_tcp_windowsize(int inSock, int inTCPWin, int inSend)
{
#  ifdef SO_SNDBUF
	int rc;
	int newTCPWin;

	// assert( inSock >= 0 );

	if (inTCPWin > 0) {

#    ifdef TCP_WINSHIFT

		/* UNICOS requires setting the winshift explicitly */
		if (inTCPWin > 65535) {
			int winShift = 0;
			int scaledWin = inTCPWin >> 16;
			while (scaledWin > 0) {
				scaledWin >>= 1;
				winShift++;
			}

			/* set TCP window shift */
			rc = setsockopt(inSock, IPPROTO_TCP, TCP_WINSHIFT, (char*)&winShift, sizeof(winShift));
			if (rc < 0) {
				return rc;
			}

			/* Note: you cannot verify TCP window shift, since it returns
			 * a structure and not the same integer we use to set it. (ugh) */
		}
#    endif /* TCP_WINSHIFT  */

#    ifdef TCP_RFC1323
		/* On AIX, RFC 1323 extensions can be set system-wide,
		 * using the 'no' network options command. But we can also set them
		 * per-socket, so let's try just in case. */
		if (inTCPWin > 65535) {
			/* enable RFC 1323 */
			int on = 1;
			rc = setsockopt(inSock, IPPROTO_TCP, TCP_RFC1323, (char*)&on, sizeof(on));
			if (rc < 0) {
				return rc;
			}
		}
#    endif /* TCP_RFC1323 */

		if (!inSend) {
			/* receive buffer -- set
			 * note: results are verified after connect() or listen(),
			 * since some OS's don't show the corrected value until then. */
			newTCPWin = inTCPWin;
			rc = setsockopt(inSock, SOL_SOCKET, SO_RCVBUF, (char*)&newTCPWin, sizeof(newTCPWin));
		}
		else {
			/* send buffer -- set
			 * note: results are verified after connect() or listen(),
			 * since some OS's don't show the corrected value until then. */
			newTCPWin = inTCPWin;
			rc = setsockopt(inSock, SOL_SOCKET, SO_SNDBUF, (char*)&newTCPWin, sizeof(newTCPWin));
		}
		if (rc < 0) {
			return rc;
		}
	}
#  endif /* SO_SNDBUF */

	return 0;
} /* end setsock_tcp_windowsize */

bool TcpConnection::client_check()
{
	return (g_client_id_cam[g_port_offset] != -1 && g_client_id_data[g_port_offset] != -1);
	// check_socket(g_client_id_cam) || check_socket(g_client_id_data);
}

bool TcpConnection::server_check()
{
	return (g_server_id_cam[g_port_offset] != -1 && g_server_id_data[g_port_offset] != -1);
	// check_socket(g_server_id_cam) || check_socket(g_server_id_data);
}

bool TcpConnection::is_error()
{
	return g_connection_error != 0;
}

bool TcpConnection::init_wsa()
{
#  ifdef WIN32
	WSADATA wsaData;
	// Request Winsock version 2.2
	if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
		WSACleanup();
		return false;
	}

#  endif
	return true;
}

void TcpConnection::init_port()
{
	if (g_port_offset == -1) {
		for (int tid = 0; tid < MAX_CONNECTIONS; tid++) {
			g_server_id_cam[tid] = -1;
			g_client_id_cam[tid] = -1;
			g_server_id_data[tid] = -1;
			g_client_id_data[tid] = -1;
		}

		g_connection_error = 0;
		g_port_offset = 0;
	}
}

void TcpConnection::close_wsa()
{
#  ifdef WIN32
	WSACleanup();
#  endif
}

bool TcpConnection::server_create(int port,
	int& server_id,
	int& client_id,
	sockaddr_in& server_sock,
	sockaddr_in& client_sock,
	bool only_accept)
{
	init_port();

	if (!only_accept) {
		if (!init_wsa()) {
			return false;
		}

		int type = SOCK_STREAM;
		int protocol = IPPROTO_TCP;

		//#  ifdef WITH_SOCKET_UDP
		//		type = SOCK_DGRAM;
		//		protocol = IPPROTO_UDP;
		//#  endif

		server_id = socket(AF_INET, type, protocol);

		if (server_id == -1) {
			printf("server_id == -1\n");
			fflush(0);
			return false;
		}

#  if !defined(__MIC__) && !defined(WIN32)
		int enable = 1;
		setsockopt(server_id, SOL_SOCKET, SO_REUSEPORT, &enable, sizeof(int));
#  endif

		// timeval tv;
		// tv.tv_sec = g_timeval_sec;
		// tv.tv_usec = 0;
		// if (setsockopt(server_id, SOL_SOCKET, SO_RCVTIMEO, (char *)&tv, sizeof(tv)) < 0) {
		//  printf("setsockopt == -1\n");
		//  fflush(0);
		//  return false;
		//}

		// sockaddr_in sock_name;
		memset(&server_sock, 0, sizeof(server_sock));
		memset(&client_sock, 0, sizeof(client_sock));
		server_sock.sin_family = AF_INET;
		server_sock.sin_port = htons(port);
		server_sock.sin_addr.s_addr = INADDR_ANY;

		int err_bind = bind(server_id, (sockaddr*)&server_sock, sizeof(server_sock));
		if (err_bind == -1) {
			printf("err_bind == -1\n");
			fflush(0);
			return false;
		}

		//#  ifdef WITH_SOCKET_UDP
		//		client_id = server_id;
		//#  else

		int err_listen = listen(server_id, 1);
		if (err_listen == -1) {
			printf("err_listen == -1\n");
			fflush(0);
			return false;
		}
		//#    if defined(WITH_SOCKET_ONLY_DATA)
		//		return true;
		//#    endif
	}

	sockaddr_in client_info;
	socklen_t addr_len = sizeof(client_info);

	printf("listen on %d\n", port);

	client_id = accept(server_id, (sockaddr*)&client_info, &addr_len);
	if (client_id == -1) {
		printf("client_id == -1\n");
		fflush(0);
		return false;
	}
	//#  endif

		// printf("accept\n");
	printf("accept on %d <-> %d\n", port, client_info.sin_port);

	fflush(0);

	g_connection_error = 0;

	return true;
}

bool TcpConnection::client_create(const char* server_name, int port, int& client_id, sockaddr_in& client_sock)
{
	// printf("connect to %s:%d\n", server_name, port);
	init_port();

	if (!init_wsa()) {
		return false;
	}

	hostent* host = gethostbyname(server_name);
	if (host == NULL) {
		printf("host == NULL\n");
		fflush(0);
		return false;
	}

	int type = SOCK_STREAM;
	int protocol = IPPROTO_TCP;

	//#  ifdef WITH_SOCKET_UDP
	//	type = SOCK_DGRAM;
	//	protocol = IPPROTO_UDP;
	//#  endif

	client_id = socket(AF_INET, type, protocol);
	if (client_id == -1) {
		printf("client_id == -1\n");
		fflush(0);
		return false;
	}

	// timeval tv;
	// tv.tv_sec = g_timeval_sec;
	// tv.tv_usec = 0;
	// if (setsockopt(6client_id, SOL_SOCKET, SO_RCVTIMEO, (char *)&tv, sizeof(tv)) < 0) {
	//  printf("setsockopt == -1\n");
	//  fflush(0);
	//  return false;
	//}
	// netsh int tcp set global autotuninglevel=normal

#  ifdef TCP_OPTIMIZATION
#    ifdef _WIN32
#      define SIO_TCP_SET_ACK_FREQUENCY _WSAIOW(IOC_VENDOR, 23)
	int freq = 1;
	unsigned long bytes = 0;
	int result = WSAIoctl(
		client_id, SIO_TCP_SET_ACK_FREQUENCY, &freq, sizeof(freq), NULL, 0, &bytes, NULL, NULL);
	int i = 1;
	setsockopt(client_id, IPPROTO_TCP, TCP_NODELAY, (char*)&i, sizeof(i));
#    else
	int i = 1;
	setsockopt(client_id, IPPROTO_TCP, TCP_NODELAY, (char*)&i, sizeof(i));
	i = 1;
	//setsockopt(client_id, IPPROTO_TCP, TCP_QUICKACK, (char*)&i, sizeof(i));
#    endif

#    if 1 
	setsock_tcp_windowsize(client_id, TCP_WIN_SIZE_SEND, 1);
	setsock_tcp_windowsize(client_id, TCP_WIN_SIZE_RECV, 0);
#    endif
#  endif

	// sockaddr_in client_sock;
	memset(&client_sock, 0, sizeof(client_sock));
	client_sock.sin_family = AF_INET;
	client_sock.sin_port = htons(port);
	memcpy(&(client_sock.sin_addr), host->h_addr, host->h_length);

	//#  ifndef WITH_SOCKET_UDP

	int connect_count = 0;
	g_connection_error = 0;

	while (true) {
#    ifdef _WIN32
		Sleep(2);
#    else
		usleep(2000000);
#    endif

		int err_connect = connect(client_id, (sockaddr*)&client_sock, sizeof(client_sock));
		if (client_id == -1) {
			printf("disconnect\n");
			return false;
		}

		connect_count++;

		if (connect_count < 2 && err_connect == -1) {
			printf("%d: wait on server %s:%d\n", connect_count, server_name, port);

			//#      ifdef _WIN32
			//      Sleep(2);
			//#      else
			//      usleep(2000000);
			//#      endif
			continue;
		}
		g_connection_error = err_connect;
		break;
	}
	//#  endif

		// printf("connect\n");
	printf("connect to %s:%d\n", server_name, port);
	fflush(0);

	return true;
}

void TcpConnection::close_tcp(int id)
{
#  ifdef WIN32
	closesocket(id);
#  else
	close(id);
#  endif
}

void TcpConnection::client_close()
{
	//#  if 0  // ndef _WIN32
	//#    pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
	//#  endif
	for (int tid = 0; tid < MAX_CONNECTIONS; tid++) {
		// int tid = omp_get_thread_num();
		close_tcp(g_client_id_cam[tid]);
		close_tcp(g_client_id_data[tid]);

		//g_server_id_cam[tid] = -1;
		g_client_id_cam[tid] = -1;

		//g_server_id_data[tid] = -1;
		g_client_id_data[tid] = -1;
	}

	g_connection_error = 0;
	g_port_offset = -1;
}

void TcpConnection::server_close()
{
	//#  if 0  // ndef _WIN32
	//#    pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
	//#  endif
	for (int tid = 0; tid < MAX_CONNECTIONS; tid++) {
		// int tid = omp_get_thread_num();
		close_tcp(g_server_id_cam[tid]);
		close_tcp(g_server_id_data[tid]);

		g_server_id_cam[tid] = -1;
		g_server_id_data[tid] = -1;
	}

	g_connection_error = 0;
	g_port_offset = -1;

	//close_wsa();
}

void TcpConnection::init_sockets_cam(const char* server, int port_cam, int port_data, bool is_server)
{
	g_is_server = is_server;
	init_port();

	if (g_client_id_cam[g_port_offset] == -1) {
		init_wsa();
		//#  if defined(BLENDER_CLIENT)
		if (g_is_server) {

			const char* env_p_port_cam = std::getenv("SOCKET_SERVER_PORT_CAM");
			if (port_cam == 0) {
				port_cam = (env_p_port_cam) ? atoi(env_p_port_cam) : 7000;
			}

			//#    if 0  // ndef _WIN32
			//#      pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
			//#    endif
					//for (int tid = 0; tid < SOCKET_CONNECTIONS; tid++) {
						// int tid = omp_get_thread_num();
			server_create(port_cam + g_port_offset,
				g_server_id_cam[g_port_offset],
				g_client_id_cam[g_port_offset],
				g_server_sockaddr_cam[g_port_offset],
				g_client_sockaddr_cam[g_port_offset], false);
			//}

	//#    ifdef WITH_SOCKET_UDP
	//		char ack = -1;
	//		recv_data_cam(&ack, sizeof(ack), false);
	//#    endif

			init_sockets_data(server, port_data);

			//#  else
		}
		else {

			const char* env_p_port_cam = std::getenv("SOCKET_SERVER_PORT_CAM");
			if (port_cam == 0) {
				// port_cam = atoi(env_p_port_cam);
				port_cam = (env_p_port_cam) ? atoi(env_p_port_cam) : 7000;
			}

			const char* env_p_name_cam = std::getenv("SOCKET_SERVER_NAME_CAM");
			char server_temp[1024];
			strcpy(server_temp, "localhost");

			if (env_p_name_cam != NULL) {
				strcpy(server_temp, env_p_name_cam);
			}

			if (server != NULL) {
				strcpy(server_temp, server);
			}

			//#    if 0  // ndef _WIN32
			//#      pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
			//#    endif
					//for (int tid = 0; tid < SOCKET_CONNECTIONS; tid++) {
						// int tid = omp_get_thread_num();
			client_create(server_temp, port_cam + g_port_offset, g_client_id_cam[g_port_offset], g_client_sockaddr_cam[g_port_offset]);
			//}

#    ifndef WITH_CLIENT_RENDERENGINE_SENDER
			init_sockets_data(server, port_data);
#    endif

			//#    ifdef WITH_SOCKET_UDP
			//		char ack = -1;
			//		send_data_cam(&ack, sizeof(ack), false);
			//#    endif

	//#  endif
		}
	}
}

void TcpConnection::init_sockets_data(const char* server, int port, bool is_server)
{
	g_is_server = is_server;
	init_port();

	if (g_client_id_data[g_port_offset] == -1) {
		init_wsa();

		//#  if (!defined(WITH_SOCKET_ONLY_DATA) && !defined(BLENDER_CLIENT) && \
		//       !defined(WITH_CLIENT_MPI_VRCLIENT)) || \
		//      (defined(WITH_SOCKET_ONLY_DATA) && defined(BLENDER_CLIENT))
		if (!g_is_server) {

			const char* env_p_port_data = std::getenv("SOCKET_SERVER_PORT_DATA");
			if (port == 0) {
				// port = atoi(env_p_port_data);
				port = (env_p_port_data) ? atoi(env_p_port_data) : 7001;
			}

			const char* env_p_name_data = std::getenv("SOCKET_SERVER_NAME_DATA");
			char server_temp[1024];
			strcpy(server_temp, "localhost");

			if (env_p_name_data != NULL) {
				strcpy(server_temp, env_p_name_data);
			}

			if (server != NULL) {
				strcpy(server_temp, server);
			}

			//#    ifdef WITH_SOCKET_ONLY_DATA
			//		//#ifndef _WIN32
			//		//#        pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
			//		//#endif
			//		//for (int tid = 0; tid < SOCKET_CONNECTIONS; tid++) {
			//			// int tid = i;//omp_get_thread_num();
			//			//#        pragma omp critical
			//			client_create(server_temp, port + g_port_offset, g_client_id_data[g_port_offset], g_client_sockaddr_data[g_port_offset]);
			//		//}
			//#    else
			//#      if 0  // ndef _WIN32
			//#        pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
			//#      endif
					//for (int tid = 0; tid < SOCKET_CONNECTIONS; tid++) {
						// int tid = omp_get_thread_num();
			client_create(server_temp, port + g_port_offset, g_client_id_data[g_port_offset], g_client_sockaddr_data[g_port_offset]);
			//}
	//#    endif
			// char ack = -1;
			// send_data_data(&ack, sizeof(ack));

	//#  else
		}
		else {

			const char* env_p_port_data = std::getenv("SOCKET_SERVER_PORT_DATA");
			if (port == 0) {
				// port = atoi(env_p_port_data);
				port = (env_p_port_data) ? atoi(env_p_port_data) : 7001;
			}

			//#    if 0// defined(WITH_SOCKET_ONLY_DATA)
			//		server_create(port,
			//			g_server_id_data[0],
			//			g_client_id_data[0],
			//			g_server_sockaddr_data[0],
			//			g_client_sockaddr_data[0],
			//			false);
			//
			//		for (int tid = 1; tid < SOCKET_CONNECTIONS; tid++) {
			//			g_server_id_data[tid] = g_server_id_data[0];
			//			g_server_sockaddr_data[tid] = g_server_sockaddr_data[0];
			//		}
			//
			//		//#        pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
			//		for (int i = 0; i < SOCKET_CONNECTIONS; i++) {
			//			int tid = i;  // omp_get_thread_num();
			//			//#pragma omp critical
			//			server_create(port,
			//				g_server_id_data[tid],
			//				g_client_id_data[tid],
			//				g_server_sockaddr_data[tid],
			//				g_client_sockaddr_data[tid],
			//				true);
			//		}
			//
			//#    else
			//#      if 0  // ndef _WIN32
			//#        pragma omp parallel for num_threads(SOCKET_CONNECTIONS)
			//#      endif
					//for (int tid = 0; tid < SOCKET_CONNECTIONS; tid++) {
						// int tid = omp_get_thread_num();
			server_create(port + g_port_offset,
				g_server_id_data[g_port_offset],
				g_client_id_data[g_port_offset],
				g_server_sockaddr_data[g_port_offset],
				g_client_sockaddr_data[g_port_offset], false);
			//}
			// char ack = -1;
			// recv_data_data(&ack, sizeof(ack));
	//#    endif
	//#  endif
		}
	}
}

#  define DEBUG_PRINT(size) //printf("%s: %lld\n", __FUNCTION__, size);

void TcpConnection::send_data_cam(char* data, size_t size, bool ack_enabled)
{
	DEBUG_PRINT(size);

	init_sockets_cam();

	if (is_error())
		return;

	size_t sended_size = 0;

	while (sended_size != size) {
		size_t size_to_send = size - sended_size;
		if (size_to_send > TCP_MAX_SIZE) {
			size_to_send = TCP_MAX_SIZE;
		}

		int temp = KERNEL_SOCKET_SEND(g_client_id_cam[g_port_offset], (char*)data + sended_size, size_to_send);

		if (temp < 1) {
			g_connection_error = 1;
			break;
		}

		sended_size += temp;
	}

	if (ack_enabled) {
		char ack = 0;
		KERNEL_SOCKET_RECV_IGNORE_RC(g_client_id_cam[g_port_offset], &ack, 1);
		if (ack != 0) {
			printf("error in send_data_cam\n");
			g_connection_error = 1;
		}
	}
}

void TcpConnection::send_data_data(char* data, size_t size, bool ack_enabled)
{
	DEBUG_PRINT(size);

	init_sockets_data();

	if (is_error())
		return;

	size_t sended_size = 0;

	while (sended_size != size) {
		size_t size_to_send = size - sended_size;
		if (size_to_send > TCP_MAX_SIZE) {
			size_to_send = TCP_MAX_SIZE;
		}

		int temp = KERNEL_SOCKET_SEND(g_client_id_data[g_port_offset], (char*)data + sended_size, size_to_send);

		if (temp < 1) {
			g_connection_error = 1;
			break;
		}

		sended_size += temp;
	}

	if (ack_enabled) {
		char ack = 0;
		KERNEL_SOCKET_RECV_IGNORE_RC(g_client_id_data[g_port_offset], &ack, 1);
		if (ack != 0) {
			printf("error in g_client_id_data\n");
			g_connection_error = 1;
		}
	}
}

void TcpConnection::recv_data_cam(char* data, size_t size, bool ack_enabled)
{
	DEBUG_PRINT(size);

	init_sockets_cam();

	if (is_error())
		return;

	size_t sended_size = 0;

	while (sended_size != size) {
		size_t size_to_send = size - sended_size;
		if (size_to_send > TCP_MAX_SIZE) {
			size_to_send = TCP_MAX_SIZE;
		}

		int temp = KERNEL_SOCKET_RECV(g_client_id_cam[g_port_offset], (char*)data + sended_size, size_to_send);

		if (temp < 1) {
			g_connection_error = 1;
			break;
		}

		sended_size += temp;
	}

	if (ack_enabled) {
		char ack = 0;
		KERNEL_SOCKET_SEND_IGNORE_RC(g_client_id_cam[g_port_offset], &ack, 1);
		if (ack != 0) {
			printf("error in g_client_id_cam\n");
			g_connection_error = 1;
		}
	}
}

void TcpConnection::recv_data_data(char* data, size_t size, bool ack_enabled)
{
	DEBUG_PRINT(size);

	init_sockets_data();

	if (is_error())
		return;

	size_t sended_size = 0;

	while (sended_size != size) {
		size_t size_to_send = size - sended_size;
		if (size_to_send > TCP_MAX_SIZE) {
			size_to_send = TCP_MAX_SIZE;
		}

		int temp = KERNEL_SOCKET_RECV(g_client_id_data[g_port_offset], (char*)data + sended_size, size_to_send);

		if (temp < 1) {
			g_connection_error = 1;
			break;
		}

		sended_size += temp;
	}

	if (ack_enabled) {
		char ack = 0;
		KERNEL_SOCKET_SEND_IGNORE_RC(g_client_id_data[g_port_offset], &ack, 1);
		if (ack != 0) {
			printf("error in g_client_id_data\n");
			g_connection_error = 1;
		}
	}
}

// limit UDP 65,507 bytes

//void TcpConnection::send_data(char* data, size_t size)
//{
//	send_data_data(data, size);
//}
//
//void TcpConnection::recv_data(char* data, size_t size)
//{
//	recv_data_data(data, size);
//}

void TcpConnection::close_kernelglobal()
{
	client_close();
}

void TcpConnection::write_data_kernelglobal(void* data, size_t size)
{
	send_data_data((char*)data, size);
}

bool TcpConnection::read_data_kernelglobal(void* data, size_t size)
{
	recv_data_cam((char*)data, size);

	return true;
}

void TcpConnection::set_port_offset(int offset)
{
	init_port();
	g_port_offset = offset;
}

