import socket

def _blocked(*args, **kwargs):
    raise RuntimeError('OFFLINE CHECK: network access blocked')

socket.socket.connect = _blocked
socket.socket.connect_ex = _blocked
socket.create_connection = _blocked
