#####################################################################################################################
# Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
#
# This program is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import socket
import select
import socketserver as SocketServer
from . import bspace_pref

def connect_to_server(context):
    """connect_to_server"""

    pref = bspace_pref.preferences()

    port = pref.bspace_port
    server_hostname = pref.bspace_server_name

    try:
        context.scene.bspace.connection_channel = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        context.scene.bspace.connection_channel.connect((server_hostname, port))

    except Exception as ex:
        print(ex)
        context.scene.bspace.disconnect(context)