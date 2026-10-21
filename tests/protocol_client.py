#!/usr/bin/env python3
"""Minimal TCP client for the space-converter remote protocol (smoke test T8).

Mirrors the wire format used by the bspace addon (see the protocol comment in
src/data_communication.cpp `recv_requested_data` and bspace_panel.py). Connects
to a running server, requests dataset info, then one sparse extraction, and
checks that a parseable payload with sane metadata comes back.

Usage: protocol_client.py <host> <port> [particle_type] [block_id]
Exit code 0 on success, 1 on failure.
"""

import socket
import struct
import sys

MSG_EXIT, MSG_EMPTY, MSG_INFO, MSG_DATA, MSG_BBOX = -1, 0, 1, 2, 3


def recvall(sock, n):
    buf = b""
    while len(buf) < n:
        part = sock.recv(n - len(buf))
        if not part:
            raise ConnectionError(f"connection closed while reading {n} bytes, got {len(buf)}")
        buf += part
    return buf


def send(sock, fmt, *vals):
    sock.sendall(struct.pack(fmt, *vals))


def main():
    host = sys.argv[1] if len(sys.argv) > 1 else "localhost"
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
    ptype = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    block = int(sys.argv[4]) if len(sys.argv) > 4 else 0

    sock = socket.create_connection((host, port), timeout=120)
    sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)

    # ── INFO request ──────────────────────────────────────────────────────────
    send(sock, "<i", MSG_INFO)
    anim_type, anim_start, anim_end = struct.unpack("<3i", recvall(sock, 12))
    (name_len,) = struct.unpack("<i", recvall(sock, 4))
    names = recvall(sock, name_len).decode(errors="replace")
    send(sock, "<i", 1)  # ack
    print(f"[client] info: anim_type={anim_type} range=[{anim_start},{anim_end}]")
    print(f"[client] particle data types: {names[:200]}{'...' if len(names) > 200 else ''}")
    if name_len <= 0:
        print("[client] FAIL: empty type/block list")
        return 1

    # ── DATA request (sparse extraction, full box) ────────────────────────────
    send(sock, "<i", MSG_DATA)
    send(sock, "<3f", 0.0, 0.0, 0.0)          # bbox_min (object space)
    send(sock, "<3f", 1000.0, 1000.0, 1000.0)  # bbox_max
    send(sock, "<i", 64)                       # bbox_dim
    send(sock, "<f", 1.0)                      # grid_transform
    send(sock, "<i", ptype)                    # particle_type
    send(sock, "<i", block)                    # block_name_id
    send(sock, "<i", 0)                        # extracted_type = sparse
    send(sock, "<i", 0)                        # dense_type = none
    send(sock, "<i", 0)                        # dense_norm = none
    send(sock, "<f", 1000.0)                   # object_size
    send(sock, "<f", 0.0)                      # particle_radius_multiplier
    send(sock, "<f", -3.4e38)                  # filter_min
    send(sock, "<f", 3.4e38)                   # filter_max
    send(sock, "<i", 0)                        # frame
    send(sock, "<i", 0)                        # anim_type
    send(sock, "<i", 0)                        # anim_task_counter

    (file_type,) = struct.unpack("<i", recvall(sock, 4))
    (size,) = struct.unpack("<q", recvall(sock, 8))
    print(f"[client] receiving file_type={file_type}, {size} bytes")
    if not (0 < size < (1 << 31)):
        print("[client] FAIL: implausible payload size")
        return 1
    payload = recvall(sock, size)
    min_v, max_v, min_r, max_r = struct.unpack("<4f", recvall(sock, 16))
    (frames,) = struct.unpack("<i", recvall(sock, 4))
    send(sock, "<i", 1)  # ack
    print(f"[client] min/max={min_v:g}/{max_v:g} reduced={min_r:g}/{max_r:g} frames={frames}")

    ok = True
    if file_type == 1 and not payload[:8].startswith(b"\x20\x42\x44\x56"):  # ' BDV' magic of .vdb
        print(f"[client] FAIL: OpenVDB payload lacks VDB magic (got {payload[:8]!r})")
        ok = False
    if not (max_r >= min_r):
        print("[client] FAIL: reduced max < min")
        ok = False
    if frames != 1:
        print(f"[client] FAIL: expected 1 frame, got {frames}")
        ok = False

    # ── EXIT ──────────────────────────────────────────────────────────────────
    send(sock, "<i", MSG_EXIT)
    sock.close()
    print("[client] OK" if ok else "[client] FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
